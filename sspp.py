import os, sys
import numpy as np
import pandas as pd
import logging
import argparse
import configparser
import time
import tracemalloc
import sqlite3 as sql

from modules import PPTranslateMagnitude
from modules import PPAddUncertainties
from modules import PPRandomizeMeasurements
from modules import PPTrailingLoss
from modules import PPFootprintFilter
from modules import PPOutWriteCSV
from modules import PPOutWriteSqlite3
from modules import PPVignetting

# def loadPandas(path, delim_whitespace=False):
#     suffix = path.split(',')[-1]
#     if suffix == '.h5':
#         return pd.read_hdf(path).reset_index(drop=True)
#     elif suffix == '.csv':
#         return pd.read_csv(path, delim_whitespace=delim_whitespace)
#     else:
#         print('input path suffix could not be parsed. File format may be incorrect')

def checkFootprint(df, cadence, fp, recordDetector=False):
    onSensor, detectorIDs = fp.applyFootprint(df, cadence)
    out = df.iloc[onSensor].reset_index(drop=True, inplace=True)
    if recordDetector:
        out['detector'] = detectorIDs
    return out

def run(parser):
    args=parser.parse_args()
    path2cadence = args.cadence
    path2input = args.input
    path2colors = args.colors

    if args.outputstem=='input':
        outstem = path2input.split('.')[0]
    else:
        outstem = args.outputstem

    if args.outpath == 'None':
        outpath = outstem
    else:
        outpath = os.path.join(args.outpath, outstem)

    rng = np.random.default_rng(2021) # should randomize this

    # pplogger.info('Reading pointing database')
     # print(cadence.columns)
    # oif = loadPandas(path2input)
    if args.format=='hdf':
        if args.h5table==None: # assumes pandas formatted hdf5 file
            oif = pd.read_hdf(args.input).reset_index(drop=True)
        else:
            oif = pd.read_hdf(args.input, key=args.h5table).reset_index(drop=True)
            # note that while we can load table from such an h5 file, they must be written to seperate output. 
            # adds a bit of time, but not complicated to have another post-post-processing script
    elif args.format=='csv':
        oif = pd.read_csv( args.inout )
    else: # for now assume csv
        oif = pd.read_csv( args.inout )

    cadencestart = oif['FieldID'].min()-1
    cadenceend= oif['FieldID'].max()+1

    cadencesuffix = path2cadence.split('.')[-1]
    querycolumns= ['observationId', 'observationStartMJD', 'filter', 'seeingFwhmGeom', 'seeingFwhmEff', 'fiveSigmaDepth', 'fieldRA', 'fieldDec', 'rotSkyPos']
    if cadencesuffix=='db':
        query = 'SELECT ' + ', '.join(querycolumns) + f' FROM observations ORDER BY observationId LIMIT %i, %i' %(cadencestart, cadenceend-cadencestart+1)
        con=sql.connect(path2cadence)
        cadence=pd.read_sql_query(query, con)
    elif cadencesuffix=='csv':
        skiprows=lambda x: x > 0 and x < cadencestart and x > cadenceend
        # nrows=cadenceend-cadencestart
        cadence = pd.read_csv(path2cadence, usecols=querycolumns, skiprows=skiprows, header=0)
    else:
        sys.exit('unrecognized cadence file suffix. terminating...')
    colors=pd.read_csv(path2colors, delim_whitespace=True).reset_index(drop=True)

    fpPath = os.path.join(os.path.dirname(__file__), "data/detectors_corners.csv")
    fovPath = os.path.join(os.path.dirname(__file__), "data/cameraFootprintNoGaps.csv")
    footprint = PPFootprintFilter.Footprint(fpPath)
    fov = PPFootprintFilter.Footprint(fovPath)

    # surveydb_join= pd.merge(oif["FieldID"], cadence, left_on="FieldID", right_on="observationId", how="left")
    # for name in ["fiveSigmaDepth", 'filter', 'seeingFwhmEff', 'seeingFwhmGeom']:
    #     oif[name] = surveydb_join[name]
    oif = pd.merge(oif, cadence, left_on="FieldID", right_on="observationId", how="left")

    #==========================================================================
    # do all the processing here

    # drop observations not inside the camera outline
    oif = checkFootprint(oif, cadence, fov)

    # calculate apparent magnitude
    oif["Mag"]=PPTranslateMagnitude.PPTranslateMagnitude(oif, cadence, colors)

    # calculate uncertainties
    oif['AstrometricSigma(mas)'], oif['PhotometricSigma(mag)'], oif["SNR"] = PPAddUncertainties.addUncertainties(oif, cadence, obsIdNameEph='FieldID', obsIdName='observationId', filterMagName='Mag')
    oif["AstrometricSigma(deg)"] = oif['AstrometricSigma(mas)'] / 3600 / 1000

    oif['dmagDetect']=PPTrailingLoss.PPTrailingLoss(oif, cadence)
    oif['dmagVignet']=PPVignetting.vignettingLosses(oif, cadence)

    oif = checkFootprint(oif, cadence, footprint, recordDetector=True)

    oif.drop(columns=["AstrometricSigma(mas)"], inplace=True)
    oif['Mag_fading'] = PPRandomizeMeasurements.randomizePhotometry(oif, magName='Mag', sigName="PhotometricSigma(mag)", rng=rng)

    if args.save_footprint:
        if args.h5table==None:
            fp_out_path=os.path.join(args.outpath, outstem+'_footprint'+'.h5')
        else:
            fp_out_path=os.path.join(args.outpath, outstem+'_'+args.h5table+'_footprint'+'.h5')
        oif.drop(columns=cadence.columns).to_hdf(
            fp_out_path,
            key='data', 
            complevel=3,
            format='table')
        
    oif.drop( np.where(oif["SNR"] <= 2.)[0], inplace=True)
    oif.reset_index(drop=True, inplace=True)

    oif["AstRATrue(deg)"] = oif["AstRA(deg)"]
    oif["AstDecTrue(deg)"] = oif["AstDec(deg)"]
    oif["AstRA(deg)"], oif["AstDec(deg)"] = PPRandomizeMeasurements.randomizeAstrometry(oif, sigName='AstrometricSigma(deg)', rng=rng)

    oif.drop( np.where(oif["Mag_fading"] + oif["dmagDetect"] + oif['dmagVignet'] >= oif["fiveSigmaDepth"])[0], inplace=True)
    oif.reset_index(drop=True, inplace=True)

    #==========================================================================

    if args.h5table==None:
        out_path_final=outpath + '.h5'
    else:
        out_path_final=outpath + '_' + args.h5table + '.h5'
        # fp_out_path=os.path.join(args.outpath, outstem+'_'+args.h5table+'_footprint'+'.h5')
    out_path_final = outpath + '.h5'

    columns2drop = []
    for column in cadence.columns:
        if np.isin(column, oif.columns):
            columns2drop.append(column)

    oif.drop(columns=columns2drop).to_hdf(
            out_path_final,
            key='data', format='table', index=False, mode='w',
            complevel=3,)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='path to input file', type=str)
    parser.add_argument('--cadence', help='path to survey database', type=str)
    parser.add_argument('--colors', help='path to color file', type=str)
    parser.add_argument('--save-footprint', dest='save_footprint', help='saves a copy after the footprint is run', action=argparse.BooleanOptionalAction)
    parser.add_argument('--outputstem', help='output file name', type=str, default='input')
    parser.add_argument('--format', help='output file format', type=str, default='csv')
    parser.add_argument('--outpath', type=str, default='None')
    parser.add_argument('--h5table', type=str, default=None, help='If input file is hdf5 with multiple tables, which table to read from')

    t0=time.time()
    tracemalloc.start()
    run(parser)
    _, peak = tracemalloc.get_traced_memory()
    t1=time.time()
    print(f"Peak memory use was %3.3f GB, runtime was %3.3f seconds" %(peak / 10**9, (t1-t0) ))
    # parser =
    # parser = argparse.ArgumentParser('--no-footprint', help='whether to run camera footprint', action=argparse.BooleanOptionalAction)
    # parser = argparse.ArgumentParser('--no-cutoff', help='disables limiting magnitude cutoff', action=argparse.BooleanOptionalAction)
    # parser = argparse.ArgumentParser('--no-cutoff', help='disables limiting magnitude cutoff', action=argparse.BooleanOptionalAction)

main()
    
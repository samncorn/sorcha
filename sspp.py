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

def run(parser):
    args=parser.parse_args()
    path2opsim = args.opsim
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
    con=sql.connect(path2opsim)
    opsim=pd.read_sql_query('SELECT observationId, observationStartMJD, filter, seeingFwhmGeom, seeingFwhmEff, fiveSigmaDepth, fieldRA, fieldDec, rotSkyPos FROM observations order by observationId', con)
    # print(opsim.columns)
    # oif = loadPandas(path2input)
    oif = pd.read_hdf(args.input).reset_index(drop=True)
    colors=pd.read_csv(path2colors, delim_whitespace=True).reset_index(drop=True)

    # surveydb_join= pd.merge(oif["FieldID"], opsim, left_on="FieldID", right_on="observationId", how="left")
    # for name in ["fiveSigmaDepth", 'filter', 'seeingFwhmEff', 'seeingFwhmGeom']:
    #     oif[name] = surveydb_join[name]
    oif = pd.merge(oif, opsim, left_on="FieldID", right_on="observationId", how="left")

    # calculate 
    oif["Mag"]=PPTranslateMagnitude.PPTranslateMagnitude(oif, opsim, colors)

    # print(oif.columns)
    oif['AstrometricSigma(mas)'], oif['PhotometricSigma(mag)'], oif["SNR"] = PPAddUncertainties.addUncertainties(oif, opsim)
    # print(len(uncert[0]))
    # print(len(oif))
    oif["AstrometricSigma(deg)"] = oif['AstrometricSigma(mas)'] / 3600 / 1000

    oif['dmagDetect']=PPTrailingLoss.PPTrailingLoss(oif, opsim)
    oif['dmagVignet']=PPVignetting.vignettingLosses(oif, opsim)

    footprint = PPFootprintFilter.Footprint("detectors_corners.csv")
    onSensor, detectorIDs = footprint.applyFootprint(oif, opsim)
    oif=oif.iloc[onSensor]
    oif["detectorID"] = detectorIDs
    oif.reset_index(drop=True, inplace=True)

    if args.save_footprint:
        oif.drop(columns=opsim.columns).to_hdf(
            os.path.join(args.outpath, outstem, '.h5'),
            key='data', 
            complevel=3)
        
    oif.drop(columns=["AstrometricSigma(mas)"], inplace=True)
    oif.drop( np.where(oif["SNR"] <= 2.)[0], inplace=True)
    oif.reset_index(drop=True, inplace=True)

    oif["AstRATrue(deg)"] = oif["AstRA(deg)"]
    oif["AstDecTrue(deg)"] = oif["AstDec(deg)"]
    oif["AstRA(deg)"], oif["AstDec(deg)"] = PPRandomizeMeasurements.randomizeAstrometry(oif, sigName='AstrometricSigma(deg)', rng=rng)

    oif.drop( np.where(oif["Mag"] + oif["dmagDetect"] + oif['dmagVignet'] >= oif["fiveSigmaDepth"])[0], inplace=True)
    oif.reset_index(drop=True, inplace=True)

    out_path_final = outpath + '.h5'
    oif.drop(columns=opsim.columns).to_hdf(
            out_path_final,
            key='data', format='table', index=False, mode='w',
            complevel=3,)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='path to input file', type=str)
    parser.add_argument('--opsim', help='path to survey database', type=str)
    parser.add_argument('--colors', help='path to color file', type=str)
    parser.add_argument('--save-footprint', dest='save_footprint', help='saves a copy after the footprint is run', action=argparse.BooleanOptionalAction)
    parser.add_argument('--outputstem', help='output file name', type=str, default='input')
    # parser.add_argument('--outputformat', help='output file format', type=str)
    parser.add_argument('--outpath', type=str, default='None')

    run(parser)
    # parser =
    # parser = argparse.ArgumentParser('--no-footprint', help='whether to run camera footprint', action=argparse.BooleanOptionalAction)
    # parser = argparse.ArgumentParser('--no-cutoff', help='disables limiting magnitude cutoff', action=argparse.BooleanOptionalAction)
    # parser = argparse.ArgumentParser('--no-cutoff', help='disables limiting magnitude cutoff', action=argparse.BooleanOptionalAction)

main()
    
# A basic test case for the lightcurve model

import cProfile
import pstats

from pstats import SortKey
from sorcha.sorcha import runLSSTSimulation  # noqa: F401
from sorcha.modules.PPConfigParser import PPConfigFileParser
import argparse

if __name__ == "__main__":  # pragma: no cover
    # Parse the command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("--object_type", default="mba", help="The type of objects to test (mba or tno).")
    args = parser.parse_args()

    cmd_args_dict = {
        "paramsinput": f"./{args.object_type}_sample_1000_physical_lca.csv",
        "orbinfile": f"./{args.object_type}_sample_1000_orbit.csv",
        "oifoutput": f"./{args.object_type}_sample_1000_eph.csv",
        "configfile": "./OIFconfig_benchmark.ini",
        "pointing_database": "./baseline_v2.0_1yr.db",
        "outpath": "../tests/out",
        "makeTemporaryEphemerisDatabase": False,
        "readTemporaryEphemerisDatabase": False,
        "deleteTemporaryEphemerisDatabase": False,
        "surveyname": "LSST",
        "outfilestem": f"out_{args.object_type}",
        "verbose": False,
    }

    configs = PPConfigFileParser("./OIFconfig_benchmark.ini", "LSST")

    cProfile.run("runLSSTSimulation(cmd_args_dict, configs)", "../tests/out/restats")

    p = pstats.Stats("../tests/out/restats")
    p.strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats()

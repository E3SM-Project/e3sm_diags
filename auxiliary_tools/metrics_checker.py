#!usr/bin/env python
from __future__ import division, print_function

import argparse
import os

"""
Usage: metrics_checker.py [options]

Options:
  -t FILE, Path to test e3sm_diags results directory
  -r FILE, Path to reference e3sm_diags results directory

About:
	This script is used to compare seasonal mean tables between a rest and reference e3sm_diags run,
        and to print out lines of variables being changed in test.
"""


parser = argparse.ArgumentParser()
parser.add_argument(
    "--ref_path",
    "-r",
    dest="ref",
    help="Path to reference e3sm_diags output",
    metavar="FILE",
)
parser.add_argument(
    "--test_path",
    "-t",
    dest="test",
    help="Path to test e3sm_diags output",
    metavar="FILE",
)

args = parser.parse_args()

seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]

ref_path = args.ref
test_path = args.test
# ref_path = '/global/cfs/cdirs/e3sm/www/chengzhu/e3sm_diags_for_cmip/CMIP6_20220324_v2_paper_linear/E3SM-1-0/historical/r1i1p1f1'
# test_path = '/global/cfs/cdirs/e3sm/www/chengzhu/e3sm_diags_for_cmip/CMIP6_20220324_v2_paper_linear/E3SM-1-1/historical/r1i1p1f1'


def compare_metrics(ref_path, test_path, season):
    fref = os.path.join(ref_path, "viewer/table-data", f"{season}_metrics_table.csv")
    ftest = os.path.join(test_path, "viewer/table-data", f"{season}_metrics_table.csv")
    try:
        with open(fref, "r") as ref, open(ftest, "r") as test:
            file_ref = ref.readlines()
            file_test = test.readlines()
            for line in file_ref:
                if line not in file_test:
                    print(f"Found difference in {season}", line)
    except Exception as e: 
        print('Failed to open file:'+ str(e))



for season in seasons:
    compare_metrics(ref_path, test_path, season)

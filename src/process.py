#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  process.py
 * @Package: 
 * @Description: Control the deNSF pipeline
 * @author: tjiang
 * @date: June 11 2018
 * @version V1.0.1     
'''

import argparse
import sys
import pysam
# import extract, Map, call_TE
# from process import *

# STAGES = {"extract": extract.run, \
#           "map": Map.run, \
#           "call": call_TE.run}

STAGES = dict()

VERSION="1.0.1"

USAGE = """\
         _           __     _   ______   ______
        | |         |  \   | | |  ____| |  ____|
        | |         |   \  | | | |____  | |____
     ___| |   ____  | |\ \ | | |_____ | |  ____|
    |  _  |  / __ \ | | \ \| |  _   | | | |
    | |_| | |  ___/ | |  \   | | |__| | | |
    |_____|  \____| |_|   \__| |______| |_|

    deNSF - deBruijn Graph-based Novel Sequence Insertion Finder

  STAGE is one of
    
  See README.md for documentation or --help for details
  
  deNSF V%s
"""%(VERSION)

def parseArgs():

	parser = argparse.ArgumentParser(prog="process.py", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)

	# parser.add_argument("-h", "--help", action="store_true")
	parser.add_argument("stage", metavar="STAGE", choices=STAGES.keys(), type=str, help="Stage to execute")
	parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER, help="Options to pass to the stage")

	args = parser.parse_args()

	STAGES[args.stage](args.options)

if __name__ == '__main__':
	parseArgs()
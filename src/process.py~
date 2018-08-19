#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  process.py
 * @Package: 
 * @Description: Control rCANID pipeline
 * @author: tjiang
 * @date: June 11 2018
 * @version V1.0.1     
'''

import argparse
import sys
# import pysam
# import extract, Map, call_TE
# from process import *
import cluster
import assembly
import detection

STAGES = {"Cluster": cluster.run, \
          "Assemble": assembly.run, \
          "Detect": detection.run}

# STAGES = dict()

VERSION="1.0.1"

USAGE = """\
             _____    ____    __    _   _____   _____
      _ __  / __  \  / __ \  |  \  | | |_   _| |  __ \\
     | ^__| | | |_| / /  \ \ |   \ | |   | |   | |  \ \\
     | |    | |  _  | |__| | | |\ \| |   | |   | |  | |
     | |    | |_| | |  __  | | | \   |  _| |_  | |__/ /
     |_|    \_____/ |_|  |_| |_|  \__| |_____| |_____/

  rCANID - read Clustering and Assembliy-based Novel insertion Detection tool

  STAGE is one 
    Cluster  cluster all of signal reads and unmapped reads respectively
    Assemble generate high-quality contigs for each cluster
    Detect   detect novel sequence insertions
    
  See README.md for documentation or --help for details
  
  rCANID V%s
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

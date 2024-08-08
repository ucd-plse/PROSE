#!/usr/bin/env python3

import argparse
import os, sys; sys.path.append(os.path.join(os.path.dirname(__file__)))
from proselib import ProseProjectTransformer

# Initialize CLI argument parser
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-s', metavar='setup_file_path', type=str, required=True,
                    help='setup.ini file path that describes how to build, run,'
                         ' and evaluate the codebase to be tuned.')
args = parser.parse_args()

myProseProjectTransformer = ProseProjectTransformer(args.s)
myProseProjectTransformer.preliminary_analysis()
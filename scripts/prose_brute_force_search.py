#!/usr/bin/env python3

import argparse
import os
from proselib import BruteForceSearch, ProseProjectTransformer

# Initialize and execute CLI argument parser
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-s', metavar='setup_file_path', type=str, required=True,
                    help='setup.ini file path that describes how to build, run,'
                        ' and evaluate the codebase to be tuned.')
args = parser.parse_args()

prose = ProseProjectTransformer(args.s)
prose.preliminary_analysis()
search_algorithm = BruteForceSearch(prose.generate_search_space(args.s))
prose.search(search_algorithm)
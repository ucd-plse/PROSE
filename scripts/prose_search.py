#!/usr/bin/env python3

import argparse
import os
from proselib import PrecimoniousSearch, ProseProjectTransformer

# Initialize and execute CLI argument parser
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-s', metavar='setup_file_path', type=str, required=True,
                    help='setup.ini file path that describes how to build, run,'
                        ' and evaluate the codebase to be tuned.')
args = parser.parse_args()

# if there are saved transformers, load them
try:
    saved_transformer = [x for x in os.listdir("prose_workspace/") if x.endswith("ProseProjectTransformer.pckl")][0]
    prose = ProseProjectTransformer.load(path_to_transformer=os.path.join("prose_workspace", saved_transformer), resume=True)

    # if there is a saved search algorithm object, load it -- we are resuming an interrupted search
    if os.path.exists("prose_workspace/__PrecimoniousSearch.pckl"):
        search_algorithm = PrecimoniousSearch.load()
        
    # otherwise, we are assuming that this script is being called after prose_preliminary.py
    else:
        search_algorithm = PrecimoniousSearch(prose.generate_search_space(args.s))

# otherwise, run the whole pipeline
except (FileNotFoundError, IndexError):
    prose = ProseProjectTransformer(args.s)
    prose.preliminary_analysis()
    search_algorithm = PrecimoniousSearch(prose.generate_search_space(args.s))

prose.search(search_algorithm)
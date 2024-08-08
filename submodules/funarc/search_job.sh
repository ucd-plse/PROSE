#!/bin/bash

# reset experiment directory
make reset

# commence search
python3 ${PROSE_REPO_PATH}/scripts/prose_search.py -s setup.ini
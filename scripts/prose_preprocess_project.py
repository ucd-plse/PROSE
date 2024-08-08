#!/usr/bin/env python3
import sys
from proselib import preprocess_project

if len(sys.argv) > 1:
    preprocess_project(sys.argv[1])
else:
    preprocess_project()
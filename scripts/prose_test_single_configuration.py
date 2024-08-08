#!/usr/bin/env python3

import os
import sys
from proselib import ProseProjectTransformer

prose = ProseProjectTransformer.load(path_to_transformer=os.path.join("prose_workspace", [x for x in os.listdir("prose_workspace/") if x.endswith("ProseProjectTransformer.pckl")][0]))
if len(sys.argv) > 2:
    prose.test_configuration(sys.argv[1], sys.argv[2])
else:
    prose.test_configuration(sys.argv[1])

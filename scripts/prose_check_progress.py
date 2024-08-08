#!/usr/bin/env python3

import os
from proselib import ProseProjectTransformer

prose = ProseProjectTransformer.load(path_to_transformer=os.path.join("prose_workspace", [x for x in os.listdir("prose_workspace/") if x.endswith("ProseProjectTransformer.pckl")][0]))
prose.report()

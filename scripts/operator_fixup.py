#!/usr/bin/env python3
import sys
import os
import re

# not accounting for line continuation
operator_overload_re = re.compile(r"^\s*public\s+(assignment|operator, (//|[<>/][=]?|\.\w+\.|==), \2|operator, (==|\*|\+|-|/))", re.IGNORECASE)

input_file_name = sys.argv[1]
if input_file_name.endswith(".rmod"):
    input_file_name = input_file_name.lower()

if not os.path.isfile(input_file_name):
    print("** module name doesn't match file name. rmod file not fixed with respect to operator overloading.")
else:
    with open(input_file_name, "r") as f:
        lines = f.readlines()

    i = -1
    while i + 1 < len(lines):
        i += 1
        if operator_overload_re.match(lines[i]):

            j = i
            while lines[i].rstrip().endswith("&"):
                j += 1
                lines[i] = lines[i].rstrip()[:-1] + lines[j].lstrip()
                lines[j] = ""

            lines[i] = re.sub(r"assignment", "assignment(=)", lines[i], flags=re.IGNORECASE)
            lines[i] = re.sub(r"operator, (//|[<>/][=]?|\.\w+\.|==), \1", r"operator(\1)", lines[i], flags=re.IGNORECASE)
            lines[i] = re.sub(r"operator, (\*|\+|-|/)", r"operator(\1)", lines[i], flags=re.IGNORECASE)

    with open(input_file_name, "w") as f:
        f.writelines(lines)
#!/usr/bin/env python3

from copy import deepcopy
import os
import re
import subprocess
import shutil
from glob import glob
from parsing import preprocess, find_valid_fortran_names

FORTRAN_DECLARATION_MODIFIERS = [
    "real",
    "integer",
    "common",
    "character",
    "save",
    "parameter",
    "logical",
    "intent",
    "optional",
    "pointer",
    "target",
    "allocatable",
    "dimension",
]

FORTRAN_NON_SCOPE_END_TYPES = [
    "associate",
    "block",
    "critical",
    "do",
    "enum",
    "forall",
    "if",
    "select",
    "team", 
    "where"
]

INTRINSIC_OR_OMITTED_PROC = "intrinsic or omitted"

FIXED_FORM_FORTRAN = False

PROCESSED_SRC_PATHS = {}
NAME_TO_SRC_PATH_MAP = {}
TARGET_SRC_CODE_SLICES = {}

fp_literal_re1 = re.compile(r"(?<![a-z0-9_])-?([0-9]+\.[0-9]*|[0-9]+\.[0-9]+)([ed]-?[0-9]*)?(_[a-z0-9_]+)?", re.IGNORECASE)
fp_literal_re2 = re.compile(r"(?<![a-z0-9_])-?([0-9]+)([ed]-?[0-9]+)(_[a-z0-9_]+)?", re.IGNORECASE)
fp_literal_nokind_re1 = re.compile(r"(?<![a-z0-9_])-?(([0-9]+\.[0-9]*)|([0-9]*\.[0-9]+))([ed]-?[0-9]*)?", re.IGNORECASE)
fp_literal_nokind_re2 = re.compile(r"(?<![a-z0-9_])-?([0-9]+)([ed]-?[0-9]+)", re.IGNORECASE)
integer_literal_nokind_re = re.compile(r"(?<![a-z0-9_])-?[0-9]+(?![\.de])")
derived_type_begin_re = re.compile(r"^\s*type(?![a-z0-9_])(.*::)?\s*([a-z][a-z0-9_]*)", re.IGNORECASE)
derived_type_var_decl_re = re.compile(r"^\s*type\s*\(\s*([a-z][a-z0-9_]*)", re.IGNORECASE)
use_import_re = re.compile(r"^\s*use(\s+|.*::\s*)([a-z][a-z0-9_]*)", re.IGNORECASE)
procedure_begin_re = re.compile(r"(^|\s)(subroutine|function|interface)\s+([a-z0-9_]+)", re.IGNORECASE)
program_begin_re = re.compile(r"^\s*program\s+([a-z][a-z0-9_]*)", re.IGNORECASE)
module_begin_re = re.compile(r"^\s*module\s+([a-z][a-z0-9_]*)", re.IGNORECASE)
line_continuation_free_form_re = re.compile(r"^.+&\s*(!.*)?$")
line_continuation_fixed_form_re = re.compile(r"^\s\s\s\s\s\S(?!0)")
comment_line_fixed_form_re = re.compile(r"^[cd!*]", re.IGNORECASE) 
comment_line_free_form_re = re.compile(r"^\s*!")
variable_declaration_re = re.compile(r"((^\s*((real[\s(,\*])|(integer[\s(,\*])|(character[\*\s(,])|(common[\s(,])|(save[\s(,](?!\s*$))|(parameter[\s(,])|(type[\s(,])|(data[\s(,])|(external[\s(,])|(double\s+precision[\s(,])|(logical[\s,])))|(::))", re.IGNORECASE)
real_variable_declaration_re = re.compile(r"^\s*real[\s(,\*]", re.IGNORECASE)
possible_proc_call_re = re.compile(r"(call\s+)?([a-z][a-z0-9_]*)\s*\(", re.IGNORECASE)
possible_string_literal_re = re.compile(r"[\"']")
open_paren_re = re.compile(r"\(")
valid_fortran_names_re = re.compile(r"(?<![a-z0-9_])([a-z][a-z0-9_]*(\s*%\s*[a-z][a-z0-9_]*)*)", re.IGNORECASE)
keep_line_re = re.compile(r"^\s*(?<!!)\s*(program|module|submodule|function|subroutine|implicit|contains|public|private|interface)", re.IGNORECASE)
cpp_directive_re = re.compile(r"^\s*#")
include_directive_re = re.compile(r"^\s*(#\s*)?include\s*[\'\"<]([a-z0-9_\./-]+)", re.IGNORECASE)
cpp_define_directive_re = re.compile(r"^\s*#\s*define\s*", re.IGNORECASE)
end_re = re.compile(r"^\s*end((\s*$)|\s+([a-z][a-z0-9_]*))", re.IGNORECASE)
prose_wrapper_name_re = re.compile(r"[a-zA-Z][a-zA-Z0-9_]*_(wrapper_(id[0-9]+_)?[048x]+_to_[048x]+|wrap_[0-9]+)")
sign_intrinsic_re = re.compile(r"[,\s]sign[\s/(]", re.IGNORECASE)
public_statement_re = re.compile(r"^\s*public\s+[a-z][a-z0-9_]*", re.IGNORECASE)
implicit_statement_re = re.compile(r"^\s*implicit\s", re.IGNORECASE)
contains_statement_re = re.compile(r"^\s*contains\s", re.IGNORECASE)
fixed_sign_call_arg_string_re = re.compile(r"\(REAL\((.*),8\)\),\(REAL\((.*),8\)\)\)")

# using some ugly parsing here to avoid splitting matching string literal delimiters which have shown up in many if statements
stmt_comment_removal = lambda stmt : stmt[:stmt.find("!")] + "\n" if "!" in stmt and stmt[:stmt.find("!")].count("'")%2 == 0 and stmt[:stmt.find("!")].count('"')%2 == 0 else stmt

comment_line = lambda line : comment_line_fixed_form_re.search(line) if (FIXED_FORM_FORTRAN and not comment_line_free_form_re.search(line)) else comment_line_free_form_re.search(line)
possibly_contains_string_literal = lambda line : possible_string_literal_re.search(stmt_comment_removal(line))
contains_end_statement = lambda line: end_re.search(stmt_comment_removal(line))
procedure_begin = lambda line : procedure_begin_re.search(line) if procedure_begin_re.search(stmt_comment_removal(line)) else None

def get_scoped_name(src_lines, start_idx, var_name, containing_scope_name=""):

    if containing_scope_name == "":
        if "%" in var_name:
            containing_derived_type_specific_name = var_name.split("%")[-2].strip().lower()
            var_name = var_name.split("%")[-1].strip().lower()
            i = start_idx + 1
            while i - 1 >= 0:
                i -= 1

                # check for declaration of this variable to get the name of the derived type
                mmatch = derived_type_var_decl_re.search(src_lines[i])
                if mmatch:
                    if containing_derived_type_specific_name in src_lines[i][mmatch.end():].lower():
                        containing_scope_name = mmatch.group(1).lower()
                        break

        else:
            # read previous lines
            i = start_idx + 1

            # first check to see if the variable is declared in this src file
            while i - 1 >= 0:
                i -= 1
                if variable_declaration_re.search(src_lines[i]):
                    non_paren_string = ""
                    j = -1
                    while j + 1 < len(src_lines[i]):
                        j += 1
                        if src_lines[i][j] == "(":
                            while j < len(src_lines[i]) and src_lines[i][j] != ")":
                                j += 1
                        else:
                            non_paren_string += src_lines[i][j]

                    valid_fortran_names = [mmatch.group(0).lower().replace(" ", "") for mmatch in valid_fortran_names_re.finditer(non_paren_string) if mmatch]
                    if var_name in valid_fortran_names:
                        break
            
            # once we know the variable was declared in this src file, look for containing scope information
            while i - 1 >= 0:
                i -= 1

                if comment_line(src_lines[i]):
                    continue

                # check for procedure
                mmatch = procedure_begin(src_lines[i])
                if mmatch and not contains_end_statement(src_lines[i]) and not possibly_contains_string_literal(src_lines[i][:mmatch.end(0)]):
                    containing_scope_name = mmatch.group(3).lower()
                    break
                        
                # check for derived type
                mmatch =  derived_type_begin_re.search(src_lines[i])
                if mmatch and not derived_type_var_decl_re.search(src_lines[i]):
                    containing_scope_name = mmatch.group(2).lower()
                    break

                # check for module
                mmatch =  module_begin_re.search(src_lines[i])
                if mmatch:
                    containing_scope_name = mmatch.group(1).lower()
                    break

    return containing_scope_name + "::" + var_name


def semicolons_to_newlines(src_lines):
    
    temp_lines = []
    for i in range(len(src_lines)):

        if not comment_line(src_lines[i]):

            # if there is a semicolon in the line...
            if ";" in src_lines[i]:

                # ... and if there is a comment character before it, skip it
                if "!" in src_lines[i] and src_lines[i].find("!") < src_lines[i].rfind(";"):
                    pass

                # ... and if it is in a string literal, skip it
                elif src_lines[i][:src_lines[i].find(";")].count('"')%2 == 1 or src_lines[i][:src_lines[i].find(";")].count("'")%2 == 1:
                    pass

                # otherwise, split it up and put each statement on a new line
                else:
                    for l in src_lines[i].split(";"):
                        if l.strip():
                            
                            # accounting for fixed-form continuation lines
                            if len(l) - len(l.lstrip()) == 5 and l.lstrip()[0] in ["*", "&", "$"]:
                                whitespace = "     "
                            else:
                                whitespace = "      "
                            temp_lines.append(whitespace + l.strip() + "\n")
                    continue
            
        temp_lines.append(src_lines[i])

    return temp_lines


def comment_line_length_fix(src_lines):
    for i in range(len(src_lines)):
        if len(src_lines[i]) > 132 and (comment_line(src_lines[i])):
            src_lines[i] = ""
    return src_lines


def gather_statement_text(src_lines, init_idx):

    gathered_line = ""

    # fixed form
    if FIXED_FORM_FORTRAN:
        if comment_line_fixed_form_re.search(src_lines[init_idx]) or not line_continuation_fixed_form_re.search(src_lines[init_idx]):
            return init_idx, init_idx, src_lines[init_idx]
        
        # find potential statement end
        # conditional tests for any line continuation characters in column 5, blank lines, cpp directives, or comment lines
        last_idx = init_idx
        while last_idx < len(src_lines) and (line_continuation_fixed_form_re.search(src_lines[last_idx]) or src_lines[last_idx].strip() == "" or comment_line_fixed_form_re.search(src_lines[last_idx]) or src_lines[last_idx].strip().startswith("#")):
            last_idx += 1
        
        # find actual statement end backtracking from potential statement end, 
        # this will be the first statement we encounter that starts with a line continuation character in column 5
        while last_idx > 0 and not line_continuation_fixed_form_re.search(src_lines[last_idx]):
            last_idx -= 1
        first_idx = last_idx + 1

        # find actual statement begin backtracking from last_idx
        # this will be the first statement we encounter -- ignoring comment lines, cpp directives, and
        # blank lines -- that doesn't start with a line continuation character in column 5 
        # also gather the line along the way
        while first_idx - 1 > 0 and line_continuation_fixed_form_re.search(src_lines[first_idx - 1]):
            first_idx -= 1
            gathered_line = stmt_comment_removal(src_lines[first_idx].strip()[1:]).strip() + gathered_line

            while first_idx - 1 > 0 and (src_lines[first_idx - 1].strip() == "" or comment_line_fixed_form_re.search(src_lines[first_idx - 1]) or src_lines[first_idx - 1].strip().startswith("#")):
                first_idx -= 1
                if comment_line_fixed_form_re.search(src_lines[first_idx]) or src_lines[first_idx].strip() == "":
                    continue
                elif src_lines[first_idx].strip().startswith("#"):
                    gathered_line = "\n" + src_lines[first_idx].strip() + "\n     & " + gathered_line
        
        first_idx -= 1
        gathered_line = stmt_comment_removal(src_lines[first_idx].rstrip()).rstrip() + gathered_line

    # free form
    else:          
        if comment_line_free_form_re.search(src_lines[init_idx]) or not line_continuation_free_form_re.search(src_lines[init_idx]):
            return init_idx, init_idx, src_lines[init_idx]

        # find actual statement end tracking forward from first_idx (which is init_idx)
        # this will be the first statement we encounter -- ignoring comment lines, cpp directives, and
        # blank lines -- that doesn't end with a line continuation character 
        # also gather the line along the way
        first_idx = init_idx
        last_idx = first_idx - 1
        while last_idx + 1 < len(src_lines) and line_continuation_free_form_re.search(src_lines[last_idx + 1]):
            last_idx += 1
            trimmed_line = stmt_comment_removal(src_lines[last_idx].strip()).strip()
            if trimmed_line.endswith("&"):
                trimmed_line = trimmed_line[:trimmed_line.rfind("&")]
            if trimmed_line.startswith("&"):
                trimmed_line = trimmed_line[trimmed_line.rfind("&") + 1:]
            gathered_line = gathered_line + trimmed_line

            while last_idx + 1 < len(src_lines) and (src_lines[last_idx + 1].strip() == "" or comment_line_free_form_re.search(src_lines[last_idx + 1]) or src_lines[last_idx + 1].strip().startswith("#")):
                last_idx += 1
                if comment_line_free_form_re.search(src_lines[last_idx]) or src_lines[last_idx].strip() == "":
                    continue
                elif src_lines[last_idx].strip().startswith("#"):
                    gathered_line = gathered_line + "&\n" + src_lines[last_idx].strip() + "\n"
                
        last_idx += 1
        trimmed_line = stmt_comment_removal(src_lines[last_idx].strip()).strip()
        if trimmed_line.endswith("&"):
            trimmed_line = trimmed_line[:trimmed_line.rfind("&")]
        if trimmed_line.startswith("&"):
            trimmed_line = trimmed_line[trimmed_line.rfind("&") + 1:]
        gathered_line = gathered_line + trimmed_line

    return first_idx, last_idx, gathered_line


def old_preprocess(src_lines):
    src_lines = comment_line_length_fix(src_lines)

    i = -1
    while i + 1 < len(src_lines):
        i += 1

        # exclude any program blocks
        match = program_begin_re.search(src_lines[i])
        if match:
            src_lines[i] = ""
            unmatched_begin_count = 1
            while i + 1 < len(src_lines) and unmatched_begin_count > 0:
                i += 1

                if not possibly_contains_string_literal(src_lines[i]) and not comment_line(src_lines[i]):

                    match_begin = procedure_begin(src_lines[i])
                    match_end = contains_end_statement(src_lines[i])

                    # second clause is a workaround for nameless interfaces
                    if (match_begin and not match_end) or (src_lines[i].lower().strip().startswith("interface") and not match_end):
                        unmatched_begin_count += 1
                    elif match_end and (match_end.group(3) == None or match_end.group(3).lower() not in FORTRAN_NON_SCOPE_END_TYPES):
                        unmatched_begin_count -= 1

                src_lines[i] = ""
            continue

        # don't touch cpp directives
        if cpp_directive_re.search(src_lines[i]):
            continue

        i, end_idx, gathered_line = gather_statement_text(src_lines, i)
        if gathered_line != src_lines[i]:
            gathered_lines = gathered_line.split("\n")
            while i <= end_idx:
                if len(gathered_lines) != 0:
                    src_lines[i] = gathered_lines.pop(0) + "\n"
                else:
                    src_lines[i] = ""
                i += 1
            i = end_idx

    src_lines = semicolons_to_newlines(src_lines)

    return src_lines


# if you provide multiple src paths, they should be in the order of dependence, if any
def slice_program_and_build_graphs(target_src_paths, src_search_paths, additional_plugin_flags):

    global FIXED_FORM_FORTRAN

    # clean up any files that may be lingering from previously failed slices
    for src_search_path in src_search_paths:
        subprocess.run(
            f"find {src_search_path} -type f -iname '*_postprocessed.f*' | xargs rm -f",
            shell=True,
            executable="/bin/bash",
        )    

    target_src_paths = [os.path.relpath(x) for x in target_src_paths]
    upstream_src_paths = []

    src_queue = deepcopy(target_src_paths)
    i = -1
    print(f"** PASS 1/6 gathering upstream procedures")
    while i + 1 < len(src_queue):
        i += 1
        this_src_path = src_queue[i]
        FIXED_FORM_FORTRAN = this_src_path.lower().endswith(".f") 

        with open(this_src_path, "r") as f:
            src_lines = preprocess(f.readlines())
            
        print(f"\t processing file {i+1: >4}/{len(src_queue): <3} {this_src_path: <120}", end='\r', flush=True)

        module_name = ""
        j = -1
        while j + 1 < len(src_lines) and module_name == "":
            j += 1

            # check for module statements
            match = module_begin_re.search(src_lines[j])
            if match:
                module_name = match.group(1).lower()
                NAME_TO_SRC_PATH_MAP[module_name] = this_src_path
                break

        for src_search_path in src_search_paths:
            result = subprocess.run(
                f"find {src_search_path} -type f -iname '*.f90' -o -type f -iname '*.f' | xargs -r -n 128 -P 16 grep -lirE '^\s*use\s+{module_name}(,|\s+|$)'",
                shell=True,
                stdout=subprocess.PIPE,
                text=True,
                executable="/bin/bash",
            )
            upstream_src_paths = list(set(upstream_src_paths + [os.path.relpath(x) for x in result.stdout.split("\n") if x]))

    src_queue = list(set(deepcopy(upstream_src_paths) + deepcopy(target_src_paths)))
    i = -1
    excluded = set()
    print(f"\n** PASS 2/6 gathering symbols from all procedures")
    while i + 1 < len(src_queue):
        i += 1
        this_src_path = src_queue[i]
        FIXED_FORM_FORTRAN = this_src_path.lower().endswith(".f") 
        scope_stack = []

        PROCESSED_SRC_PATHS[this_src_path] = {
            'all_proc_symbols'         : set(),
            'real_vars'                : set(),
            'nonreal_vars'             : set(),
            'derived_types'            : set(),
            'all_depends'              : set(),
            'minimal_depends'          : set(),
            'required_proc_symbols'    : set(),
            'required_derived_type_symbols' : set(),
        }
    
        with open(this_src_path, "r") as f:
            src_lines = preprocess(f.readlines())
            
        print(f"\t processing file {i+1: >4}/{len(src_queue): <3} {this_src_path: <120}", end='\r', flush=True)

        j = -1
        while j + 1 < len(src_lines):
            j += 1

            if comment_line(src_lines[j]):
                continue

            # check for module statements
            match = module_begin_re.search(src_lines[j])
            if match:
                module_name = match.group(1).lower()
                NAME_TO_SRC_PATH_MAP[module_name] = this_src_path
                scope_stack.append(module_name)
                continue

            # check for USE statements
            match = use_import_re.search(src_lines[j])
            if match:
                module_name = match.group(2).lower()

                if module_name in NAME_TO_SRC_PATH_MAP.keys():
                    depend_src_path = NAME_TO_SRC_PATH_MAP[module_name]
                elif module_name in excluded:
                    continue
                else:
                    unique_found = False
                    multiple_found = []
                    for src_search_path in src_search_paths:
                        result = subprocess.run(
                            f"find {src_search_path} -type f -iname '*.f90' -o -type f -iname '*.f' | xargs -r -n 128 -P 16 grep -lirE '^\s*module\s+{module_name}(\s+|$)'",
                            shell=True,
                            stdout=subprocess.PIPE,
                            text=True,
                            executable="/bin/bash",
                        )
                        grep_results = [x for x in result.stdout.split("\n") if x]
                        if len(grep_results) == 1:
                            unique_found = True
                            depend_src_path = os.path.relpath(grep_results[0])
                            NAME_TO_SRC_PATH_MAP[module_name] = depend_src_path
                            break
                        elif len(grep_results) > 1:
                            multiple_found = deepcopy(grep_results)
                    if not unique_found and multiple_found == []:
                        excluded.add(module_name)
                        continue
                    elif not unique_found and multiple_found != []:
                        print(f"multiple matches for {module_name} found: {multiple_found}")
                        assert(False)

                PROCESSED_SRC_PATHS[this_src_path]['all_depends'].add(depend_src_path)
                if depend_src_path not in PROCESSED_SRC_PATHS.keys() and depend_src_path not in src_queue:
                    src_queue.append(depend_src_path)
                continue

            # check for procedure declarations
            match = procedure_begin(src_lines[j])
            if match and not contains_end_statement(src_lines[j]) and not possibly_contains_string_literal(src_lines[j][:match.end(0)]):
                proc_name = match.group(3).lower()
                PROCESSED_SRC_PATHS[this_src_path]['all_proc_symbols'].add(proc_name)
                if this_src_path in target_src_paths:
                    PROCESSED_SRC_PATHS[this_src_path]['required_proc_symbols'].add(proc_name)

                if match.group(2).lower() == "interface":
                    while j + 1 < len(src_lines):
                        j += 1

                        if comment_line(src_lines[j]):    
                            continue

                        match_end = contains_end_statement(src_lines[j])
                        if match_end and not possibly_contains_string_literal(src_lines[j][:match_end.end(0)]) and " interface" in match_end.group(0).lower():
                            break
                else:
                    scope_stack.append(proc_name)
                continue

            # check for end statements
            match = contains_end_statement(src_lines[j])
            if match and not possibly_contains_string_literal(src_lines[j][:match.end(0)]):
                if match.group(3) == None or match.group(3).lower() not in ["interface", "type"] + FORTRAN_NON_SCOPE_END_TYPES:
                    scope_stack.pop(-1)
                continue

            # check for include directives and add the src lines
            match = include_directive_re.search(src_lines[j])
            if match:
                include_file_name = match.group(2)
                if include_file_name in NAME_TO_SRC_PATH_MAP.keys():
                    depend_src_path = NAME_TO_SRC_PATH_MAP[include_file_name]
                elif include_file_name in excluded:
                    src_lines[j] = ""
                    continue
                else:
                    unique_found = False
                    multiple_found = []
                    for src_search_path in src_search_paths:
                        if "/" not in include_file_name:
                            result = subprocess.run(
                                f"find {src_search_path} -type f -iname '{include_file_name}'",
                                shell=True,
                                stdout=subprocess.PIPE,
                                text=True,
                                executable="/bin/bash",
                            )
                        else:
                            result = subprocess.run(
                                f"find {src_search_path} -type f -wholename '{include_file_name}'",
                                shell=True,
                                stdout=subprocess.PIPE,
                                text=True,
                                executable="/bin/bash",
                            )

                        find_results = [x for x in result.stdout.split("\n") if x]
                        if len(find_results) == 1:
                            unique_found = True
                            depend_src_path = os.path.relpath(find_results[0])
                            NAME_TO_SRC_PATH_MAP[include_file_name] = depend_src_path
                            with open(depend_src_path, "r") as f:
                                include_lines = preprocess(f.readlines())
                            src_lines = src_lines[:j] + ["!     !PROSE!" + src_lines[j]] + include_lines + ["!     !PROSE!" + src_lines[j]] + src_lines[j+1:]
                            break
                        elif len(find_results) > 1:
                            multiple_found = deepcopy(find_results)
                    if not unique_found and multiple_found == []:
                        excluded.add(include_file_name)
                        src_lines[j] = ""
                        continue
                    elif not unique_found and multiple_found != []:
                        print(f"multiple matches for {include_file_name} found: {multiple_found}")
                        assert(False)
                    continue

            # check for var declarations
            match = variable_declaration_re.search(src_lines[j])
            if match and not use_import_re.search(src_lines[j]): # second conditional to handle use statements of the form "USE, INTRINSIC :: IEEE_ARITHMETIC"
                
                match_derived_type = derived_type_begin_re.search(src_lines[j])
                if match_derived_type and not derived_type_var_decl_re.search(src_lines[j]):
                    derived_type_name = match_derived_type.group(2).lower()
                    PROCESSED_SRC_PATHS[this_src_path]['derived_types'].add(derived_type_name)
                    if this_src_path in target_src_paths:
                        PROCESSED_SRC_PATHS[this_src_path]['required_derived_type_symbols'].add(derived_type_name)

                else:

                    # extract variable names
                    # make this easier by first removing anything within parentheses from the line
                    non_paren_string = ""
                    k = -1
                    while k + 1 < len(src_lines[j]):
                        k += 1
                        if src_lines[j][k] == "(":
                            while k < len(src_lines[j]) and src_lines[j][k] != ")":
                                k += 1
                        else:
                            non_paren_string += src_lines[j][k]

                    valid_fortran_names = [mmatch.group(0).lower().replace(" ", "") for mmatch in valid_fortran_names_re.finditer(non_paren_string) if mmatch]
                    if scope_stack != []:
                        containing_scope = scope_stack[-1]
                    else:
                        containing_scope = ""
                    scoped_names = set([get_scoped_name(src_lines, j, name, containing_scope) for name in valid_fortran_names if name not in FORTRAN_DECLARATION_MODIFIERS])

                    # save each variable name as either "real" or "nonreal"
                    match_real = real_variable_declaration_re.search(src_lines[j])
                    if match_real:
                        PROCESSED_SRC_PATHS[this_src_path]['real_vars'] = PROCESSED_SRC_PATHS[this_src_path]['real_vars'].union(scoped_names)
                    else:
                        PROCESSED_SRC_PATHS[this_src_path]['nonreal_vars'] = PROCESSED_SRC_PATHS[this_src_path]['nonreal_vars'].union(scoped_names)
                continue

    print("\n\n\t\t include-files/modules excluded from analysis:")
    for x in excluded:
        print(f"\t\t\t{x}")
    print()

    print("\n\t\t called from targeted source:")
    called = set()
    for target_src in target_src_paths:
        called = called.union(PROCESSED_SRC_PATHS[target_src]["all_depends"])
    for x in called:
        print(f"\t\t\t{x}")
    print()

    print(f"** PASS 3/6 propagating FP flow upstream")
    i = -1
    while i + 1 < len(upstream_src_paths):
        i += 1
        this_src_path = upstream_src_paths[i]
        FIXED_FORM_FORTRAN = this_src_path.lower().endswith(".f") 

        print(f"\t processing file {i+1: >4}/{len(upstream_src_paths): <3} {this_src_path: <120}", end='\r', flush=True)       

        with open(this_src_path, "r") as f:
            src_lines = preprocess(f.readlines())

        # pass over src code
        src_contains_call_to_target_procs = False
        j = -1
        while j + 1 < len(src_lines):
            j += 1

            # skip blank and commented lines
            if src_lines[j] == "" or comment_line(src_lines[j]):
                continue

            # check for #include directives and add the src lines
            match = include_directive_re.search(src_lines[j])
            if match:
                include_file_name = match.group(2)
                if include_file_name in NAME_TO_SRC_PATH_MAP.keys():
                    with open(NAME_TO_SRC_PATH_MAP[include_file_name], "r") as f:
                        include_lines = preprocess(f.readlines())
                    src_lines = src_lines[:j] + ["!     !PROSE!" + src_lines[j]] + include_lines + ["!     !PROSE!" + src_lines[j]] + src_lines[j+1:]
                elif include_file_name in excluded:
                    src_lines[j] = ""
                else:
                    assert("Should have processed this include statement earlier" == False)
                continue

            # check each procedure declaration
            match = procedure_begin(src_lines[j])
            if match and not contains_end_statement(src_lines[j]) and not possibly_contains_string_literal(src_lines[j][:match.end(0)]):
                proc_name = match.group(3).lower()

                # if symbol is an interface, remove it since there shouldn't be any procedure calls here
                if match.group(2) == "interface":

                    unmatched_begin_count = 1
                    while unmatched_begin_count > 0:
                        j += 1
    
                        if comment_line(src_lines[j]):
                            continue
                        
                        match_begin = procedure_begin(src_lines[j])
                        match_end = contains_end_statement(src_lines[j])

                        # second clause is a workaround for nameless interfaces
                        if ((match_begin and not match_end) and not possibly_contains_string_literal(src_lines[j][:match_begin.end(0)])) or (src_lines[j].lower().strip().startswith("interface") and not match_end):
                            unmatched_begin_count += 1
                        elif match_end and not possibly_contains_string_literal(src_lines[j][:match_end.end(0)]) and (match_end.group(3) == None or match_end.group(3).lower() not in FORTRAN_NON_SCOPE_END_TYPES):
                            unmatched_begin_count -= 1

                        src_lines[j] = ""

                    continue
                
                else:
                    # iterate over the procedure
                    while j + 1 < len(src_lines):
                        j += 1

                        if comment_line(src_lines[j]) or contains_statement_re.search(src_lines[j]):
                            continue

                        # stop at nested procedures or end of the current procedure
                        match_begin = procedure_begin(src_lines[j])
                        match_end = contains_end_statement(src_lines[j])
                        if match_end and not possibly_contains_string_literal(src_lines[j][:match_end.end(0)]) and (match_end.group(3) == None or match_end.group(3).lower() not in FORTRAN_NON_SCOPE_END_TYPES):
                            break
                        elif match_begin and not possibly_contains_string_literal(src_lines[j][:match_begin.end(0)]):
                            if len(match_begin.groups()) > 2 and match_begin.group(2) != "interface":
                                j -= 1
                                break
                            else: # remove nested interface blocks
                                while j + 1 < len(src_lines):
                                    j += 1

                                    if comment_line(src_lines[j]):
                                        continue

                                    match_end = contains_end_statement(src_lines[j])
                                    src_lines[j] = ""
                                    if match_end and not possibly_contains_string_literal(src_lines[j][:match_end.end(0)]) and (match_end.group(3) != None and match_end.group(3).lower() == "interface"):
                                        break

                        # check each procedure call to see if it calls a proc symbol from the targeted scopes
                        mmatches = possible_proc_call_re.finditer(src_lines[j])
                        line_contains_call_to_target_procs = None
                        for mmatch in mmatches:
                            referenced_proc_symbol = mmatch.group(2).lower()

                            proc_symbol_origin = ""
                            for src_path in set([this_src_path]).union(PROCESSED_SRC_PATHS[this_src_path]['all_depends']):
                                if referenced_proc_symbol in PROCESSED_SRC_PATHS[src_path]['all_proc_symbols']:
                                    proc_symbol_origin = src_path
                                    break

                            # if the referenced symbol is a procedure from either this src or one of its depends, we have found a procedure call! (over-approximates, assumes public accessibility of all declared procedures)
                            if proc_symbol_origin != "":

                                if line_contains_call_to_target_procs == None:
                                    line_contains_call_to_target_procs = False

                                for target_src_path in target_src_paths:
                                    if referenced_proc_symbol in PROCESSED_SRC_PATHS[target_src_path]['all_proc_symbols']:
                                        PROCESSED_SRC_PATHS[this_src_path]['required_proc_symbols'].add(proc_name)
                                        PROCESSED_SRC_PATHS[this_src_path]['minimal_depends'].add(target_src_path)
                                        line_contains_call_to_target_procs = True
                                        src_contains_call_to_target_procs = True
                                        break

                        # remove procedure calls not to the target scope
                        if line_contains_call_to_target_procs == False:
                            src_lines[j] = ""
                        else:
                            continue
        
        # save src code if there is actual fp flow to the targeted source code
        # otherwise, remove it
        if src_contains_call_to_target_procs:
            TARGET_SRC_CODE_SLICES[this_src_path] = src_lines
        else:
            upstream_src_paths.pop(i)        
            i -= 1


    print(f"** PASS 4/6 propagating FP flow downstream")
    src_queue = list(set(deepcopy(upstream_src_paths) + deepcopy(target_src_paths)))
    required_proc_symbol_count = -1
    while sum([len(PROCESSED_SRC_PATHS[src_file_path]["required_proc_symbols"]) for src_file_path in src_queue]) > required_proc_symbol_count:
        required_proc_symbol_count = sum([len(PROCESSED_SRC_PATHS[src_file_path]["required_proc_symbols"]) for src_file_path in src_queue])

        i = -1
        while i + 1 < len(src_queue):
            i += 1
            this_src_path = src_queue[i]
            FIXED_FORM_FORTRAN = this_src_path.lower().endswith(".f") 

            print(f"\t processing file {i+1: >4}/{len(src_queue): <3} {this_src_path: <120}", end='\r', flush=True)

            # load src code if alread read; otherwise read it in and preprocess
            if this_src_path in TARGET_SRC_CODE_SLICES:
                src_lines = TARGET_SRC_CODE_SLICES[this_src_path]
            else:
                with open(this_src_path, "r") as f:
                    src_lines = preprocess(f.readlines())

            # pass over src code
            scope_stack = []
            j = -1
            while j + 1 < len(src_lines):
                j += 1

                # skip blank lines and comment lines
                if src_lines[j] == "" or comment_line(src_lines[j]):
                    continue

                # check for module statements
                match = module_begin_re.search(src_lines[j])
                if match:
                    module_name = match.group(1).lower()
                    scope_stack.append(module_name)
                    continue

                # check for #include directives and add the src lines
                match = include_directive_re.search(src_lines[j])
                if match:
                    include_file_name = match.group(2)
                    if include_file_name in NAME_TO_SRC_PATH_MAP.keys():
                        with open(NAME_TO_SRC_PATH_MAP[include_file_name], "r") as f:
                            include_lines = preprocess(f.readlines())
                        src_lines = src_lines[:j] + ["!     !PROSE!" + src_lines[j]] + include_lines + ["!     !PROSE!" + src_lines[j]] + src_lines[j+1:]
                    elif include_file_name in excluded:
                        src_lines[j] = ""
                    else:
                        assert("Should have processed this include statement earlier" == False)
                    continue

                # check each procedure declaration
                match = procedure_begin(src_lines[j])
                if match and not contains_end_statement(src_lines[j]) and not possibly_contains_string_literal(src_lines[j][:match.end(0)]):
                    proc_name = match.group(3).lower()

                    # if this procedure is in one of the original target src code paths OR
                    # if FP data has been propagated here from those original target scopes, process the procedure
                    if proc_name in PROCESSED_SRC_PATHS[this_src_path]['required_proc_symbols']:

                        # get_bounds of the procedure's statements
                        statement_end_idx = j
                        while statement_end_idx < len(src_lines):
                            statement_end_idx += 1

                            if comment_line(src_lines[statement_end_idx]):
                                continue

                            match_begin = procedure_begin(src_lines[statement_end_idx])
                            match_end = contains_end_statement(src_lines[statement_end_idx])

                            # second clause is a workaround for nameless interfaces
                            if ((match_begin and not match_end) and not possibly_contains_string_literal(src_lines[statement_end_idx][:match_begin.end(0)])) or (src_lines[statement_end_idx].lower().strip().startswith("interface") and not match_end):
                                break
                            elif match_end and not possibly_contains_string_literal(src_lines[statement_end_idx][:match_end.end(0)]) and (match_end.group(3) == None or match_end.group(3).lower() not in FORTRAN_NON_SCOPE_END_TYPES):
                                break

                        # if the required symbol is an interface, we mark everything within it as required (over-approximation)
                        if match.group(2) == "interface":
                            # get bounds of the entire interface block
                            scope_end_idx = j
                            unmatched_begin_count = 1
                            while scope_end_idx < len(src_lines) and unmatched_begin_count > 0:
                                scope_end_idx += 1

                                if comment_line(src_lines[scope_end_idx]):
                                    continue

                                match_begin = procedure_begin(src_lines[scope_end_idx])
                                match_end = contains_end_statement(src_lines[scope_end_idx])

                                # second clause is a workaround for nameless interfaces
                                if ((match_begin and not match_end) and not possibly_contains_string_literal(src_lines[scope_end_idx][:match_begin.end(0)])) or (src_lines[scope_end_idx].lower().strip().startswith("interface") and not match_end):
                                    unmatched_begin_count += 1
                                elif match_end and not possibly_contains_string_literal(src_lines[scope_end_idx][:match_end.end(0)]) and (match_end.group(3) == None or match_end.group(3).lower() not in FORTRAN_NON_SCOPE_END_TYPES):
                                    unmatched_begin_count -= 1

                            potential_proc_symbols = set()
                            while j <= scope_end_idx:
                                potential_proc_symbols = potential_proc_symbols.union(set(src_lines[j].lower().replace(",", " ").split()))
                                j += 1
                            PROCESSED_SRC_PATHS[this_src_path]['required_proc_symbols'] = PROCESSED_SRC_PATHS[this_src_path]['required_proc_symbols'].union(potential_proc_symbols.intersection(PROCESSED_SRC_PATHS[this_src_path]['all_proc_symbols']))
                            continue
                        
                        # if it is not an interface, iterate over the statements within it and
                        # check for anything that adds to the minimal depends, that is:
                        #   1. derived type dependencies in the declared variables of a required procedure
                        #   2. real-valued arguments to procedure calls
                        # remove anything that is not a variable declaration, a call with real-valued arguments,
                        # comments, a contains statement, a use statement, or a cpp directive. Also parse in "include" files.
                        else:
                            scope_stack.append(proc_name)

                            # iterate over the procedure
                            while j + 1 < statement_end_idx:
                                j += 1
                                
                                # keep comments, contains
                                if comment_line(src_lines[j]) or contains_statement_re.search(src_lines[j]):
                                    continue

                                # keep use statements for now; we will prune these as a last step
                                # once the minimal depends have been identified
                                if use_import_re.search(src_lines[j]):
                                    continue

                                # check for #include directives and add the src lines
                                mmatch = include_directive_re.search(src_lines[j])
                                if mmatch:
                                    include_file_name = mmatch.group(2)
                                    if include_file_name in NAME_TO_SRC_PATH_MAP.keys():
                                        with open(NAME_TO_SRC_PATH_MAP[include_file_name], "r") as f:
                                            include_lines = preprocess(f.readlines())
                                        src_lines = src_lines[:j] + ["!     !PROSE!" + src_lines[j]] + include_lines + ["!     !PROSE!" + src_lines[j]] + src_lines[j+1:]
                                        statement_end_idx += len(include_lines)
                                    elif include_file_name in excluded:
                                        src_lines[j] = ""
                                    else:
                                        assert("Should have processed this include statement earlier" == False)
                                    continue

                                if cpp_directive_re.search(src_lines[j]):
                                    continue

                                # check variable declarations
                                mmatch = variable_declaration_re.search(src_lines[j])
                                if mmatch and not use_import_re.search(src_lines[j]): # second conditional to handle use statements of the form "USE, INTRINSIC :: IEEE_ARITHMETIC"

                                    # special handling for possible symbol imports involved in variable declarations
                                    for mmmatch in open_paren_re.finditer(src_lines[j]):

                                        # extract the argument list
                                        arg_string = ""
                                        k = mmmatch.end() - 1
                                        unmatched_paren_count = 1
                                        while k+1 < len(src_lines[j]):
                                            k += 1
                                            arg_string += src_lines[j][k]
                                            
                                            if src_lines[j][k] == "(":
                                                unmatched_paren_count += 1
                                            elif src_lines[j][k] == ")":
                                                unmatched_paren_count -= 1
                                            
                                            if unmatched_paren_count == 0:
                                                break

                                        # get the referenced arg symbols from the argument list
                                        arg_string = arg_string.replace(" ", "")
                                        referenced_var_symbols = [mmmmatch.group(0).lower() for mmmmatch in valid_fortran_names_re.finditer(arg_string) if mmmmatch]
                                        for referenced_var_symbol in referenced_var_symbols:

                                            if "%" in referenced_var_symbol:
                                                scoped_name = get_scoped_name(src_lines, j, referenced_var_symbol, "")
                                            else:
                                                scoped_name = get_scoped_name(src_lines, j, referenced_var_symbol, scope_stack[0])

                                            # if the referenced symbol is not from this src file...
                                            if referenced_var_symbol in PROCESSED_SRC_PATHS[this_src_path]['derived_types']:
                                                PROCESSED_SRC_PATHS[this_src_path]['required_derived_type_symbols'].add(referenced_var_symbol)                                    
                                            elif not scoped_name in PROCESSED_SRC_PATHS[this_src_path]['real_vars'].union(PROCESSED_SRC_PATHS[this_src_path]['nonreal_vars']):

                                                # ...then save the depend src path from which this proc symbol was imported so we can be sure to include it in the minimal program
                                                for depend_src_path in PROCESSED_SRC_PATHS[this_src_path]['all_depends']:

                                                    if "%" in referenced_var_symbol:
                                                        scoped_name = get_scoped_name(src_lines, j, referenced_var_symbol, "")
                                                    else:
                                                        module_name = [x for x in NAME_TO_SRC_PATH_MAP.keys() if NAME_TO_SRC_PATH_MAP[x] == depend_src_path][0]
                                                        scoped_name = get_scoped_name(src_lines, j, referenced_var_symbol, module_name)

                                                    if referenced_var_symbol in PROCESSED_SRC_PATHS[depend_src_path]['derived_types']:
                                                        PROCESSED_SRC_PATHS[this_src_path]['minimal_depends'].add(depend_src_path)
                                                        PROCESSED_SRC_PATHS[depend_src_path]['required_derived_type_symbols'].add(referenced_var_symbol)                                    
                                                    elif scoped_name in PROCESSED_SRC_PATHS[depend_src_path]['real_vars'].union(PROCESSED_SRC_PATHS[depend_src_path]['nonreal_vars']):
                                                        PROCESSED_SRC_PATHS[this_src_path]['minimal_depends'].add(depend_src_path)

                                    continue

                                # check each procedure call to see if it might contain reals
                                mmatches = possible_proc_call_re.finditer(src_lines[j])
                                proc_calls_with_real_args = set()
                                for mmatch in mmatches:

                                    referenced_proc_symbol = mmatch.group(2).lower()

                                    proc_call_with_real_args = False
                                    proc_symbol_origin = INTRINSIC_OR_OMITTED_PROC
                                    for src_path in set([this_src_path]).union(PROCESSED_SRC_PATHS[this_src_path]['all_depends']):
                                        if referenced_proc_symbol in PROCESSED_SRC_PATHS[src_path]['all_proc_symbols']:
                                            proc_symbol_origin = src_path
                                            break

                                    # if the referenced symbol is a procedure from either this src or one of its depends, we have found a procedure call! (over-approximates, assumes public accessibility of all declared procedures)
                                    # also, calls to the sign intrinsic are included so we can fix the signs of them in the plugin pass
                                    if proc_symbol_origin != INTRINSIC_OR_OMITTED_PROC or referenced_proc_symbol == "sign":

                                        # extract the argument list
                                        # start by extracting the text between the following matching set of parens
                                        arg_string = ""
                                        k = mmatch.end() - 1
                                        unmatched_paren_count = 1
                                        while k + 1 < len(src_lines[j]):
                                            k += 1
                                            arg_string += src_lines[j][k]
                                            
                                            if src_lines[j][k] == "(":
                                                unmatched_paren_count += 1
                                            elif src_lines[j][k] == ")":
                                                unmatched_paren_count -= 1
                                            
                                            if unmatched_paren_count == 0:
                                                break

                                        arg_string_without_string_literals = ""
                                        k = -1
                                        while k + 1 < len(arg_string):
                                            k += 1
                                            if arg_string[k] == '"':
                                                while k + 1 < len(arg_string):
                                                    k += 1 
                                                    if arg_string[k] == '"':
                                                        if arg_string[k] == arg_string[k + 1] == '"':
                                                            k += 1
                                                        elif arg_string[k - 1] == "\\":
                                                            continue
                                                        else:
                                                            break
                                            elif arg_string[k] == "'":
                                                while k + 1 < len(arg_string):
                                                    k += 1 
                                                    if arg_string[k] == "'":
                                                        if arg_string[k] == arg_string[k + 1] == "'":
                                                            k += 1
                                                        elif arg_string[k - 1] == "\\":
                                                            continue
                                                        else:
                                                            break
                                            else:
                                                arg_string_without_string_literals += arg_string[k]

                                        # first, check the arg string for the presence of FP symbols
                                        referenced_arg_symbols = [mmmatch.group(0).lower().replace(" ", "") for mmmatch in valid_fortran_names_re.finditer(arg_string_without_string_literals) if mmmatch]
                                        imported_arg_symbols = set()
                                        for referenced_arg_symbol in referenced_arg_symbols:

                                            scoped_arg_symbol = get_scoped_name(src_lines, j, referenced_arg_symbol, "")

                                            # local nonreals
                                            if scoped_arg_symbol in PROCESSED_SRC_PATHS[this_src_path]['nonreal_vars']:
                                                continue

                                            # local reals
                                            elif scoped_arg_symbol in PROCESSED_SRC_PATHS[this_src_path]['real_vars']:
                                                proc_call_with_real_args = True
                                            
                                            else:
                                                imported_arg_symbols.add(scoped_arg_symbol)

                                                for arg_symbol in imported_arg_symbols:
                                                    for depend_src_path in PROCESSED_SRC_PATHS[this_src_path]['all_depends'].difference(set([this_src_path])):
                                                        
                                                        module_name = [x for x in NAME_TO_SRC_PATH_MAP.keys() if NAME_TO_SRC_PATH_MAP[x] == depend_src_path][0]
                                                        imported_scoped_arg_symbol = module_name + arg_symbol

                                                        if imported_scoped_arg_symbol in PROCESSED_SRC_PATHS[depend_src_path]['nonreal_vars']:
                                                            break

                                                        elif imported_scoped_arg_symbol in PROCESSED_SRC_PATHS[depend_src_path]['real_vars']:
                                                            proc_call_with_real_args = True
                                                            break
                                                
                                        # if none of the extracted arg symbols were identified as real, check the arg string for the presence of FP literals
                                        if not proc_call_with_real_args:
                                            arg_string_without_string_literals_and_parens = ""
                                            k = -1
                                            while k + 1 < len(arg_string_without_string_literals):
                                                k += 1
                                                if arg_string_without_string_literals[k] == "(" and arg_string_without_string_literals[k + 1] != "/":
                                                    while k < len(arg_string_without_string_literals) and arg_string_without_string_literals[k] != ")":
                                                        k += 1
                                                else:
                                                    arg_string_without_string_literals_and_parens += arg_string_without_string_literals[k]
                                            if fp_literal_re1.search(arg_string_without_string_literals_and_parens) or fp_literal_re2.search(arg_string_without_string_literals_and_parens):
                                                proc_call_with_real_args = True

                                        # if we ultimately identified any real arguments in the arg string...
                                        if proc_call_with_real_args:

                                            # ...save the proc call; either just the call for subroutines or the call in context for functions which rose sometimes needs for some reason?
                                            if mmatch.group(1) != None:
                                                proc_calls_with_real_args.add(mmatch.group(0) + arg_string)
                                            else: 
                                                call_in_context = src_lines[j]
                                                if call_in_context.strip().lower().startswith("if"):
                                                    conditional_end_idx = call_in_context.find("(")
                                                    unmatched_paren_count = 1
                                                    while conditional_end_idx + 1 < len(call_in_context):
                                                        conditional_end_idx += 1

                                                        if call_in_context[conditional_end_idx] == "(":
                                                            unmatched_paren_count += 1
                                                        elif call_in_context[conditional_end_idx] == ")":
                                                            unmatched_paren_count -= 1
                                                        
                                                        if unmatched_paren_count == 0:
                                                            break

                                                    if mmatch.start() > conditional_end_idx:
                                                        call_in_context = call_in_context[conditional_end_idx + 1:]
                                                    else:
                                                        call_in_context = call_in_context[:call_in_context.rfind(")") + 1] + " then; endif"
                                                proc_calls_with_real_args.add(call_in_context)                                    

                                            # remember any required imported proc_symbols
                                            if proc_symbol_origin != this_src_path and proc_symbol_origin != INTRINSIC_OR_OMITTED_PROC:
                                                PROCESSED_SRC_PATHS[this_src_path]['minimal_depends'].add(proc_symbol_origin)
                                                PROCESSED_SRC_PATHS[proc_symbol_origin]['required_proc_symbols'].add(referenced_proc_symbol)

                                            # remember any required imported arg symbols
                                            for arg_symbol in imported_arg_symbols:
                                                for depend_src_path in PROCESSED_SRC_PATHS[this_src_path]['all_depends'].difference(set([this_src_path])):
                                                    
                                                    module_name = [x for x in NAME_TO_SRC_PATH_MAP.keys() if NAME_TO_SRC_PATH_MAP[x] == depend_src_path][0]
                                                    scoped_arg_symbol = module_name + arg_symbol

                                                    if scoped_arg_symbol in PROCESSED_SRC_PATHS[depend_src_path]['real_vars'].union(PROCESSED_SRC_PATHS[depend_src_path]['nonreal_vars']):
                                                        PROCESSED_SRC_PATHS[this_src_path]['minimal_depends'].add(depend_src_path)
                                                        break

                                # include procedure calls with real arguments
                                if len(proc_calls_with_real_args) > 0:
                                    new_text = ""
                                    for proc_call in proc_calls_with_real_args:
                                        new_text = new_text + "      " + proc_call + "\n"
                                    src_lines[j] = new_text
                                    continue

                                # any line that has made it to this point will now be removed
                                src_lines[j] = ""

                            scope_stack.pop(-1)

            # save src code
            TARGET_SRC_CODE_SLICES[this_src_path] = src_lines

            # add any new dependencies to the end of the queue
            src_queue += [src_path for src_path in PROCESSED_SRC_PATHS[this_src_path]['minimal_depends'] if src_path not in src_queue]

    print(f"\n\n** PASS 5/6 performing minimal slice")
    src_queue = list(set(src_queue))
    i = -1
    while i + 1 < len(src_queue):
        i += 1
        this_src_path = src_queue[i]
        FIXED_FORM_FORTRAN = this_src_path.lower().endswith(".f") 

        print(f"\t processing file {i+1: >4}/{len(src_queue): <3} {this_src_path: <120}", end='\r', flush=True)

        src_lines = TARGET_SRC_CODE_SLICES[this_src_path]

        # pass over src code
        j = -1
        while j + 1 < len(src_lines):
            j += 1

            # skip blank lines
            if src_lines[j] == "" or comment_line(src_lines[j]):
                continue

            # keep required use statements
            match = use_import_re.search(src_lines[j])
            if match:
                module_name = match.group(2).lower()

                if module_name in NAME_TO_SRC_PATH_MAP.keys():
                    depend_src_path = NAME_TO_SRC_PATH_MAP[module_name]
                    if depend_src_path in PROCESSED_SRC_PATHS[this_src_path]['minimal_depends']:
    
                        continue

            # keep required procedures
            match = procedure_begin(src_lines[j])
            if match and not contains_end_statement(src_lines[j]) and not possibly_contains_string_literal(src_lines[j][:match.end(0)]):
                proc_name = match.group(3).lower()

                # get bounds of the block that defines this procedure
                scope_end_idx = j
                statement_end_idx = len(src_lines)
                unmatched_begin_count = 1
                while scope_end_idx < len(src_lines) and unmatched_begin_count > 0:
    
                    scope_end_idx += 1

                    if comment_line(src_lines[scope_end_idx]):
                        continue

                    match = use_import_re.search(src_lines[scope_end_idx])
                    if match:
                        module_name = match.group(2).lower()
                        if module_name in NAME_TO_SRC_PATH_MAP.keys():
                            depend_src_path = NAME_TO_SRC_PATH_MAP[module_name]
                            if depend_src_path in PROCESSED_SRC_PATHS[this_src_path]['minimal_depends']:
                                continue

                        src_lines[scope_end_idx] = ""
                        continue

                    if not possibly_contains_string_literal(src_lines[scope_end_idx]):

                        match_begin = procedure_begin(src_lines[scope_end_idx])
                        match_end = contains_end_statement(src_lines[scope_end_idx])

                        # second clause is a workaround for nameless interfaces
                        if match_begin and not match_end and not possibly_contains_string_literal(src_lines[scope_end_idx][:match_begin.end(0)]):
                            unmatched_begin_count += 1
                            statement_end_idx = min(statement_end_idx, scope_end_idx - 1)
                        elif src_lines[scope_end_idx].lower().strip().startswith("interface") and not match_end:
                            unmatched_begin_count += 1
                        elif match_end and not possibly_contains_string_literal(src_lines[scope_end_idx][:match_end.end(0)]) and (match_end.group(3) == None or match_end.group(3).lower() not in FORTRAN_NON_SCOPE_END_TYPES):
                            unmatched_begin_count -= 1

                # if this doesn't define a required procedure, remove it in its entirety
                if proc_name not in PROCESSED_SRC_PATHS[this_src_path]['required_proc_symbols']:
                    while j <= scope_end_idx:
                        src_lines[j] = ""
                        j += 1
                else:
                    j = min(scope_end_idx, statement_end_idx)
                continue

            match = derived_type_begin_re.search(src_lines[j])
            if match and not derived_type_var_decl_re.search(src_lines[j]):
                derived_type_name = match.group(2).lower()
                if derived_type_name not in PROCESSED_SRC_PATHS[this_src_path]['required_derived_type_symbols']:
                    
                    # get bounds of the block that defines this derived type to remove it
                    type_start_idx = j
                    type_end_idx = j
                    unmatched_begin_count = 1
                    while unmatched_begin_count > 0:
                        type_end_idx += 1

                        if comment_line(src_lines[type_end_idx]):
                            continue
                        
                        match_begin = derived_type_begin_re.search(src_lines[type_end_idx])
                        match_false_positive = derived_type_var_decl_re.search(src_lines[type_end_idx])
                        match_end = contains_end_statement(src_lines[type_end_idx])

                        if match_begin and not match_false_positive and not possibly_contains_string_literal(src_lines[type_end_idx][:match_begin.end(0)]):
                            unmatched_begin_count += 1
                        elif match_end and not possibly_contains_string_literal(src_lines[type_end_idx][:match_end.end(0)]) and (match_end.group(3) == None or match_end.group(3).lower() not in FORTRAN_NON_SCOPE_END_TYPES):
                            unmatched_begin_count -= 1

                    for k in range(type_start_idx, type_end_idx + 1):
                        src_lines[k] = ""
                    j = type_end_idx
                    
                    continue

            # keep (program|module|submodule|function|subroutine|implicit|contains|public|private|interface)
            match = keep_line_re.search(src_lines[j])
            if match:
                if not comment_line(src_lines[j]):    
                    continue

            # include cpp directives
            match = cpp_directive_re.search(src_lines[j])
            if match:
                continue

            # include any 'end' that does not also match one of the non-proc end types
            match = contains_end_statement(src_lines[j])
            if match and not possibly_contains_string_literal(src_lines[j][:match.end(0)]) and (match.group(3) == None or match.group(3).lower() not in FORTRAN_NON_SCOPE_END_TYPES):
                continue

            # include variable declarations
            match = variable_declaration_re.search(src_lines[j])
            if match and not use_import_re.search(src_lines[j]): # second conditional to handle use statements of the form "USE, INTRINSIC :: IEEE_ARITHMETIC"
                continue

            # any line that has made it to this point will now be removed
            src_lines[j] = ""

        # if there is any src code remaining after the reduction, write it out!
        if len([line for line in src_lines if line.strip()]) == 0:
            continue
        else:
            # save a copy of the original
            shutil.copy(this_src_path, os.path.join(os.environ['PROSE_EXPERIMENT_DIR'], "prose_workspace/original_files", os.path.basename(this_src_path) + ".orig"))

            # lowercase slice except for all non-macro-defining cpp directives and any string literals
            for j in range(len(src_lines)):
                if cpp_define_directive_re.search(src_lines[j]):
                    src_lines[j] = src_lines[j].lower()
                elif cpp_directive_re.search(src_lines[j]):
                    continue
                else:
                    double_quote_count = 0
                    single_quote_count = 0
                    for k in range(len(src_lines[j])):
                        
                        if src_lines[j][k] == '"':
                            
                            # handle escaped quotes
                            if k > 0 and src_lines[j][k-1] == "\\":
                                pass
                            
                            else:
                                double_quote_count += 1

                        elif src_lines[j][k] == "'":

                            # handle apostrophes within a string literal
                            if double_quote_count % 2 == 1:
                                pass

                            # handle escaped quotes
                            elif k > 0 and src_lines[j][k-1] == "\\":
                                pass

                            else:
                                single_quote_count += 1

                        # if we are not in the middle of a string literal, go ahead and lower the char
                        if double_quote_count % 2 == 0 and single_quote_count % 2 == 0:
                            src_lines[j] = src_lines[j][:k] + src_lines[j][k].lower() + src_lines[j][k+1:]

            # write out the slice to the actual src directory
            with open(this_src_path, "w") as f:
                f.writelines(src_lines)

            # also save a copy of the slice
            shutil.copy(this_src_path, os.path.join(os.environ['PROSE_EXPERIMENT_DIR'], "prose_workspace/original_files", os.path.basename(this_src_path) + ".slice"))

    # find which minimal src files are required to process the target files
    i = -1
    target_src_paths += upstream_src_paths
    while i + 1 < len(target_src_paths):
        i += 1
        for depend_src_path in PROCESSED_SRC_PATHS[target_src_paths[i]]['minimal_depends']:
            if depend_src_path not in target_src_paths:
                target_src_paths.append(depend_src_path)
            
    # find the order in which target files should be processed
    process_order = []
    while len(target_src_paths) > 0:
        i = 0
        before = len(process_order)
        while i < len(target_src_paths):
            if PROCESSED_SRC_PATHS[target_src_paths[i]]['minimal_depends'] <= set(process_order).union(set([target_src_paths[i]])):
                process_order.append(target_src_paths[i])
                target_src_paths.pop(i)
            else:
                i += 1
        after = len(process_order)
        if before == after:
            print("\n\nCircular or missing dependencies detected")
            assert(False)

    print(f"\n\n** PASS 6/6 building graphs")
    for i, this_src_path in enumerate(process_order):
        print(f"\t processing file {i+1: >4}/{len(process_order): <3} {this_src_path: <120}", end='\r', flush=True)

        include_dirs = list(set([os.path.relpath(os.path.dirname(os.path.abspath(s)), start=os.path.dirname(this_src_path)) for s in PROCESSED_SRC_PATHS[this_src_path]['minimal_depends']]))
        if len(include_dirs) == 1:
            include_dirs = f"-I{include_dirs[0]}"
        elif len(include_dirs) > 1:
            include_dirs = "-I" + " -I".join(include_dirs)
        else:
            include_dirs = ""

        command = [
            "rose-compiler",
            "-rose:skip_syntax_check",
            "-rose:skipfinalCompileStep",
            f"-rose:plugin_lib {os.environ['PROSE_PLUGIN_PATH']}/ProsePlugin.so",
            "-rose:plugin_action prose-generate-graph",
            f"-rose:plugin_arg_prose-generate-graph {os.environ['PROSE_EXPERIMENT_DIR']}",
            "-DROSE_COMP",
            include_dirs,
            f"-I{os.path.join(os.environ['PROSE_EXPERIMENT_DIR'], 'prose_workspace/rmod_files')}",
            additional_plugin_flags,
            f"{os.path.basename(this_src_path)}",
            "&& rm *postprocessed*",
        ]

        subprocess.run(
            " ".join(command),
            shell=True,
            env=os.environ.copy(),
            check=True,
            cwd = os.path.abspath(os.path.dirname(this_src_path))
        )

        # postprocess and move any generated rmod files
        for rmod_file_name_with_path in glob(os.path.join(os.path.abspath(os.path.dirname(this_src_path)), "*.rmod")):
            rmod_file_name = os.path.basename(rmod_file_name_with_path)
            subprocess.run(
                f"operator_fixup.py {rmod_file_name} && mv {rmod_file_name} {os.path.join(os.environ['PROSE_EXPERIMENT_DIR'], 'prose_workspace/rmod_files')}",
                shell=True,
                env=os.environ.copy(),
                check=True,
                cwd = os.path.abspath(os.path.dirname(this_src_path))
            )


def unslice(transformed_src_path, original_src_path, SETUP):

    global FIXED_FORM_FORTRAN
    FIXED_FORM_FORTRAN = transformed_src_path.lower().endswith(".f") 

    with open(transformed_src_path, "r") as f:
        lines = preprocess(src_lines=f.readlines(), SETUP=SETUP, rose_preprocessing=False, excluded_names=[], fixed_form_fortran=FIXED_FORM_FORTRAN)

    with open(original_src_path, "r") as f:
        llines = preprocess(src_lines=f.readlines(), SETUP=SETUP, rose_preprocessing=False, excluded_names=[], fixed_form_fortran=FIXED_FORM_FORTRAN)

    i = ii = -1
    while i + 1 < len(lines):
        i += 1
    
        if module_begin_re.search(lines[i]):
            break
        else:
            match = procedure_begin(lines[i])
            if match and not contains_end_statement(lines[i]) and not possibly_contains_string_literal(lines[i][:match.end(0)]):
                break

    # the only time this shouldn't be executed is the exceptional case in which the source code file is empty after running the C preprocessor
    if i + 1 < len(lines):
        merge_changes_in_this_scope(lines, i, llines, ii)

    # save transformed slice
    shutil.copy(transformed_src_path, transformed_src_path + ".transformed_slice")
    
    # write out unsliced program with transformations
    with open(transformed_src_path, "w") as f:
        f.writelines(llines)


def merge_changes_in_this_scope(lines, i, llines, ii):

    # get name of this scope
    match_module = module_begin_re.search(lines[i])
    match_proc = procedure_begin(lines[i])
    match_derived_type = derived_type_begin_re.search(lines[i])
    if match_module:
        scope_name = match_module.group(1).lower()
    elif match_proc and not contains_end_statement(lines[i]) and not possibly_contains_string_literal(lines[i][:match_proc.end(0)]):
        scope_name = match_proc.group(3).lower()
    elif match_derived_type and not derived_type_var_decl_re.search(lines[i]):
        scope_name = match_derived_type.group(2).lower()
    else:
        assert("Unexpected scope passed to merge_changes_in_this_scope()" == False)

    # update index into original file lines to match scope
    sscope_name = ""
    ii -= 1     # in case it was already passed a "matching" index
    while ii + 1 < len(llines) and sscope_name != scope_name:
        ii += 1

        if comment_line(llines[ii]):
            continue

        mmatch_module = module_begin_re.search(llines[ii])
        mmatch_proc = procedure_begin(llines[ii])
        mmatch_derived_type = derived_type_begin_re.search(llines[ii])
        if mmatch_module:
            sscope_name = mmatch_module.group(1).lower()
        elif mmatch_proc and not contains_end_statement(llines[ii]) and not possibly_contains_string_literal(llines[ii][:mmatch_proc.end(0)]):
            sscope_name = mmatch_proc.group(3).lower()
        elif mmatch_derived_type and not derived_type_var_decl_re.search(llines[ii]):
            sscope_name = mmatch_derived_type.group(2).lower()

    # if these scopes are an interface, skip it as they are left unchanged
    if match_proc and match_proc.group(2).lower() == "interface":
        while i + 1 < len(lines):
            i += 1

            if comment_line(lines[i]):
                continue

            match_end = contains_end_statement(lines[i])
            if match_end and not possibly_contains_string_literal(lines[i][:match_end.end(0)]) and " interface" in match_end.group(0).lower():
                break
        while ii + 1 < len(llines):
            ii += 1

            if comment_line(llines[ii]):
                continue

            mmatch_end = contains_end_statement(llines[ii])
            if mmatch_end and not possibly_contains_string_literal(llines[ii][:mmatch_end.end(0)]) and " interface" in mmatch_end.group(0).lower():
                break

        return i, ii
    
    # otherwise, find the ends of these scopes
    else:

        scope_end_idx = i
        unmatched_begin_count = 1
        while scope_end_idx + 1 < len(lines) and unmatched_begin_count > 0:
            scope_end_idx += 1

            if comment_line(lines[scope_end_idx]):
                continue

            match_proc = procedure_begin(lines[scope_end_idx])
            match_derived_type = derived_type_begin_re.search(lines[scope_end_idx])
            match_derived_type_false_positive = derived_type_var_decl_re.search(lines[scope_end_idx])
            match_end = contains_end_statement(lines[scope_end_idx])

            # last conditional checks for nameless interfaces
            if not match_end and ((match_proc and not possibly_contains_string_literal(lines[scope_end_idx][:match_proc.end(0)])) or (match_derived_type and not match_derived_type_false_positive and not possibly_contains_string_literal(lines[scope_end_idx][:match_derived_type.end(0)])) or lines[scope_end_idx].lower().strip().startswith("interface") or lines[scope_end_idx].lower().strip().startswith("abstract ")):
                unmatched_begin_count += 1
            elif match_end and not possibly_contains_string_literal(lines[scope_end_idx][:match_end.end(0)]) and (match_end.group(3) == None or match_end.group(3).lower() not in FORTRAN_NON_SCOPE_END_TYPES):
                unmatched_begin_count -= 1

        sscope_end_idx = ii
        unmatched_begin_count = 1
        while sscope_end_idx + 1 < len(llines) and unmatched_begin_count > 0:
            sscope_end_idx += 1

            if comment_line(llines[sscope_end_idx]):
                continue

            if not possibly_contains_string_literal(llines[sscope_end_idx]):

                mmatch_proc = procedure_begin(llines[sscope_end_idx])
                mmatch_derived_type = derived_type_begin_re.search(llines[sscope_end_idx])
                mmatch_derived_type_false_positive = derived_type_var_decl_re.search(llines[sscope_end_idx])
                mmatch_end = contains_end_statement(llines[sscope_end_idx])

                # last conditional checks for nameless interfaces
                if not mmatch_end and ((mmatch_proc and not possibly_contains_string_literal(llines[sscope_end_idx][:mmatch_proc.end(0)])) or (mmatch_derived_type and not mmatch_derived_type_false_positive and not possibly_contains_string_literal(llines[sscope_end_idx][:mmatch_derived_type.end(0)])) or llines[sscope_end_idx].lower().strip().startswith("interface") or llines[sscope_end_idx].lower().strip().startswith("abstract ")):
                    unmatched_begin_count += 1
                elif mmatch_end and not possibly_contains_string_literal(llines[sscope_end_idx][:mmatch_end.end(0)]) and (mmatch_end.group(3) == None or mmatch_end.group(3).lower() not in FORTRAN_NON_SCOPE_END_TYPES):
                    unmatched_begin_count -= 1

    new_import_statements = set()

    # iterate through the lines in the transformed scope
    while i + 1 <= scope_end_idx:
        i += 1

        # skip over 'include' content
        if "!PROSE!" in lines[i]:
            i += 1
            while i + 1 <= scope_end_idx and "!PROSE!" not in lines[i]:
                i += 1
            continue

        if comment_line(lines[i]):
            continue

        # check for nested procedure scopes
        match = procedure_begin(lines[i])
        if match and not contains_end_statement(lines[i]) and not possibly_contains_string_literal(lines[i][:match.end(0)]) and not match.group(2).lower() == "interface":

            # if it is a generated wrapper procedure, we must insert it
            # otherwise, recursively call this procedure
            if prose_wrapper_name_re.search(lines[i]):

                i -= 1
                while i + 1 <= scope_end_idx:
                    i += 1
                    
                    llines[sscope_end_idx-1] = llines[sscope_end_idx-1] + lines[i]

                    if comment_line(lines[i]):
                        continue

                    match_end = contains_end_statement(lines[i])
                    if match_end and not possibly_contains_string_literal(lines[i][:match_end.end(0)]) and (match_end.group(3) == None or match_end.group(3).lower() not in FORTRAN_NON_SCOPE_END_TYPES):
                        break

            else:
                i, ii = merge_changes_in_this_scope(lines, i, llines, ii)
            
            continue

        match_derived_type = derived_type_begin_re.search(lines[i])
        match_derived_type_false_positive = derived_type_var_decl_re.search(lines[i])
        match_end = contains_end_statement(lines[i])
        if not match_end and match_derived_type and not match_derived_type_false_positive and not possibly_contains_string_literal(lines[scope_end_idx][:match_derived_type.end(0)]):
            i, ii = merge_changes_in_this_scope(lines, i, llines, ii)
            continue

        # check for any new use statements
        match = use_import_re.search(lines[i])
        if match and prose_wrapper_name_re.search(lines[i]) and lines[i] not in new_import_statements:
            new_import_statements.add(lines[i])
            llines[ii] = llines[ii] + lines[i]
            continue

        # add all public statements for generated prose wrappers
        match = public_statement_re.search(lines[i])
        if match and prose_wrapper_name_re.search(lines[i]):
            insertion_point = ii
            while insertion_point + 1 <= sscope_end_idx:
                insertion_point += 1
                if public_statement_re.search(llines[insertion_point]):
                    llines[insertion_point] = llines[insertion_point] + lines[i]
                    break
                elif implicit_statement_re.search(llines[insertion_point]):
                    llines[insertion_point] = llines[insertion_point] + lines[i]
                    break
                elif variable_declaration_re.search(llines[insertion_point]) and not use_import_re.search(llines[insertion_point]): # second conditional to handle use statements of the form "USE, INTRINSIC :: IEEE_ARITHMETIC"
                    llines[insertion_point-1] = llines[insertion_point-1] + lines[i]
                    break
                elif contains_statement_re.search(llines[insertion_point]):
                    insertion_point -= 1
                    llines[insertion_point] = llines[insertion_point] + lines[i]
                    break

        # once we encounter the first variable declaration
        # replace all real variable declarations in the scope of the original file
        # with those from the transformed file
        match = variable_declaration_re.search(lines[i])
        if match and not use_import_re.search(lines[i]): # second conditional to handle use statements of the form "USE, INTRINSIC :: IEEE_ARITHMETIC"

            # insert all real variable declarations from the transformed file
            var_declarations_to_insert = {}
            i -= 1
            while i + 1 <= scope_end_idx:
                i += 1

                # skip over 'include' content
                if "!PROSE!" in lines[i]:
                    i += 1
                    while i + 1 <= scope_end_idx and "!PROSE!" not in lines[i]:
                        i += 1
                    continue

                # skip comment lines, blank lines, cpp directives, and end type statements
                if comment_line(lines[i]) or lines[i].strip() == "" or cpp_directive_re.search(lines[i]) or (contains_end_statement(lines[i]) and "type" in lines[i].lower()):
                    continue
                elif variable_declaration_re.search(lines[i]):
                    if real_variable_declaration_re.search(lines[i]):
                        var_name = [mmatch.group(0) for mmatch in valid_fortran_names_re.finditer(lines[i][lines[i].find("::") + len("::"):]) if mmatch][0].lower()
                        var_declarations_to_insert[var_name] = lines[i] 
                else:
                    i -= 1
                    break
            
            # find first variable declaration of this scope in the original file
            # and then replace all real variable declarations
            insertion_idx = -1
            while ii + 1 <= sscope_end_idx:
                ii += 1

                mmatch = variable_declaration_re.search(llines[ii])
                if mmatch and not use_import_re.search(llines[ii]): # second conditional to handle use statements of the form "USE, INTRINSIC :: IEEE_ARITHMETIC"
                    ii -= 1
                    while ii + 1 <= sscope_end_idx:
                        ii += 1

                        # skip comment lines, blank lines, cpp directives, and include directives
                        if comment_line(llines[ii]) or llines[ii].strip() == "" or cpp_directive_re.search(llines[ii]) or include_directive_re.search(llines[ii]) or (contains_end_statement(llines[ii]) and "type" in llines[ii].lower()):
                            continue
                        elif variable_declaration_re.search(llines[ii]):
                            if insertion_idx < 0:
                                insertion_idx = ii
                            if real_variable_declaration_re.search(llines[ii]):
                                old_line = stmt_comment_removal(llines[ii])
                                remaining_var_names = list(var_declarations_to_insert.keys())

                                # only remove lines that contain a variable that we are overwriting
                                for var_name in remaining_var_names:
                                    regex = r"(?<![a-z0-9_])" + var_name + r"(?![a-z0-9_])"
                                    if re.search(regex, old_line, re.IGNORECASE):
                                        llines[ii] = ""
                                        break
                                if not llines[ii]:
                                    for var_name in remaining_var_names:
                                        regex = r"(?<![a-z0-9_])" + var_name + r"(?![a-z0-9_])"
                                        if re.search(regex, old_line, re.IGNORECASE):
                                            llines[ii] = llines[ii] + var_declarations_to_insert[var_name]
                                            var_declarations_to_insert.pop(var_name)
                        else:
                            ii -= 1
                            break
                    break

            continue

        # replace any transformed procedure calls: either those that now have a wrapper or those that involve a "sign" intrinsic
        # or those that involve a literal tagged with a kind 4
        for match in possible_proc_call_re.finditer(lines[i]):
            proc_name = match.group(2).lower()
            
            fixed_sign_call_arg_string_match = None
            if proc_name == "sign":
                fixed_sign_call_arg_string_match = fixed_sign_call_arg_string_re.search(lines[i])
            prose_wrapper_match = prose_wrapper_name_re.search(proc_name)

            arg_string = ""
            j = match.end() - 1
            unmatched_paren_count = 1
            while j + 1 < len(lines[i]):
                j += 1
                arg_string += lines[i][j]
                
                if lines[i][j] == "(":
                    unmatched_paren_count += 1
                elif lines[i][j] == ")":
                    unmatched_paren_count -= 1
                
                if unmatched_paren_count == 0:
                    break
    
            if prose_wrapper_match or fixed_sign_call_arg_string_match or re.search(r"[0-9]_4(?![0-9])", arg_string):                
                arg_string_temp = arg_string
                if fixed_sign_call_arg_string_match:
                    arg_string_temp = arg_string.replace(fixed_sign_call_arg_string_match.group(0), f"{fixed_sign_call_arg_string_match.group(1)},{fixed_sign_call_arg_string_match.group(2)})")

                arg_string_temp = arg_string_temp.lower()
                arg_tokens = [re.sub(r"\s", "", name.lower()) for name in find_valid_fortran_names(arg_string_temp)]
                arg_tokens += [m.group(0) for m in fp_literal_nokind_re1.finditer(arg_string_temp)]
                arg_tokens += [m.group(0) for m in fp_literal_nokind_re2.finditer(arg_string_temp)]
                arg_tokens += [m.group(0) for m in integer_literal_nokind_re.finditer(arg_string_temp)]

                # find the corresponding procedure call in the original src file
                replaced = False
                while ii < sscope_end_idx:
                    for mmatch in possible_proc_call_re.finditer(llines[ii]):
                        pproc_name = mmatch.group(2).lower()

                        if prose_wrapper_match and pproc_name != proc_name[:proc_name.rfind("_wrap")]:
                            continue
                        elif proc_name == "sign" and pproc_name != "sign":
                            continue

                        aarg_string = ""
                        jj = mmatch.end() - 1
                        unmatched_paren_count = 1
                        while jj + 1 < len(llines[ii]):
                            jj += 1
                            aarg_string += llines[ii][jj]
                            
                            if llines[ii][jj] == "(":
                                unmatched_paren_count += 1
                            elif llines[ii][jj] == ")":
                                unmatched_paren_count -= 1
                            
                            if unmatched_paren_count == 0:
                                break

                        aarg_string = aarg_string.lower()
                        aarg_tokens = [re.sub(r"\s", "", name.lower()) for name in find_valid_fortran_names(aarg_string)]
                        aarg_tokens += [m.group(0) for m in fp_literal_nokind_re1.finditer(aarg_string)]
                        aarg_tokens += [m.group(0) for m in fp_literal_nokind_re2.finditer(aarg_string)]
                        aarg_tokens += [m.group(0) for m in integer_literal_nokind_re.finditer(aarg_string)]

                        if sorted(arg_tokens) == sorted(aarg_tokens):
                            try: # subroutine with a "call" statement
                                llines[ii] = llines[ii][:mmatch.start()+len(mmatch.group(1))] + proc_name + "(" + arg_string + llines[ii][mmatch.end() + len(aarg_string):]
                            except TypeError: # function calls without "call" statement
                                llines[ii] = llines[ii][:mmatch.start()] + proc_name + "(" + arg_string + llines[ii][mmatch.end() + len(aarg_string):]
                            replaced = True
                            break
                    
                    if not replaced:
                        ii += 1
                    else:
                        break

    return scope_end_idx, sscope_end_idx

if __name__ == "__main__":
    import sys
    unslice(sys.argv[1],sys.argv[2])
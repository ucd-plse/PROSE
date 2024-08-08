### TODO: need to add instrumentation after fortran "statement functions" (how to identify these??) (MOM_mixed_layer_restrat.F90)
### TODO: fix bug where variable declarations inside of #ifdef blocks are included in prose instrumentation without the same cpp directives; in general blocks like this are not handled
### TODO: add support for sampling across procedure invocations in case arrays change size
### TODO: include "SAVE" keyword

import subprocess
import os
import re

MODULE_NAME_LIST = set()
FORTRAN_DECLARATION_MODIFIERS = ["real", "intent", "optional", "pointer", "target", "allocatable", "dimension"]

module_begin_re = re.compile(r"^\s*module\s", re.IGNORECASE)

use_import_re = re.compile(r"^\s*use\s+", re.IGNORECASE)
use_only_import_re = re.compile(r"^\s*use\s.*only\s*:", re.IGNORECASE)

derived_type_begin_re = re.compile(r"^\s*type[,\s](.*::)?", re.IGNORECASE)
derived_type_end_re = re.compile(r"^\s*end\s*type(\s+|$)", re.IGNORECASE)

procedure_begin_re = re.compile(r"((^|\s)subroutine\s|(^|\s)function\s)", re.IGNORECASE)

interface_begin_re = re.compile(r"^\s*(abstract\s*)?interface", re.IGNORECASE)
interface_end_re = re.compile(r"^\s*end\s+interface(\s+|$)", re.IGNORECASE)

line_continuation_free_form_re = re.compile(r"^.+&\s*(!.*)?$")
line_continuation_fixed_form_re = re.compile(r"\s\s\s\s\s[$&*]")

variable_declaration_re = re.compile(r"((^\s*((real[\s(,\*])|(integer[\s(,\*])|(character[\*\s(,])|(common[\s(,])|(save[\s(,])|(parameter[\s(,])|(logical[\s,])))|(::))", re.IGNORECASE)

# the following regex are to be used on variable declarations that, in the case of multiline variable declarations, any
# & characters, comments, or newlines have been removed and all of the lines are consolidated into a single string to be searched
fp_array_var_declaration_re = re.compile(r"\s*real\s*(\(.*?\))?(([,\s]intent\s*\(.*?\))|([,\s]dimension\s*\(.*?\)))*.*([a-z0-9_]+)\s*\(.*?\)", re.IGNORECASE)
optional_variable_declaration_re = re.compile(r"optional", re.IGNORECASE)

# the following regex is to be used with only the string of dimension info extracted from eithor of the above fp_array_var_declaration regex
deferred_or_assumed_array_shape_re = re.compile(r".*((:\s*,)|(:\s*$))")

# using some ugly parsing here to avoid splitting matching string literal delimiters which have shown up in many if statements
stmt_comment_removal = lambda stmt : stmt[:stmt.find("!")] + "\n" if "!" in stmt and stmt[:stmt.find("!")].count("'")%2 == 0 and stmt[:stmt.find("!")].count('"')%2 == 0 else stmt


def gather_statement_text(lines, i):
 
    # find start of statement
    if line_continuation_free_form_re.match(lines[i - 1]):
        while line_continuation_free_form_re.match(lines[i - 1]):
            i -= 1
    elif line_continuation_fixed_form_re.match(lines[i]):
        while line_continuation_fixed_form_re.match(lines[i]):
            i -= 1

    first_idx = i

    # gather the statement
    gathered_line = stmt_comment_removal(lines[first_idx].rstrip())
    if line_continuation_free_form_re.match(gathered_line):
        while line_continuation_free_form_re.match(gathered_line):
            i += 1
            if lines[i].strip().startswith("!"):
                continue
            else:
                gathered_line += lines[i].strip()
                gathered_line = stmt_comment_removal(gathered_line)

    elif line_continuation_fixed_form_re.match(lines[i+1]):
        while line_continuation_fixed_form_re.match(lines[i+1]) or lines[i+1].strip() == "" or lines[i+1].strip().startswith("!"):
            i += 1
            if not line_continuation_fixed_form_re.match(lines[i]):
                continue
            else:
                gathered_line += lines[i].strip()[1:]
                gathered_line = stmt_comment_removal(gathered_line)
    
    last_idx = i

    # remove any & or newline characters
    gathered_line = re.sub(r"(&|\n)", "", gathered_line)
    return first_idx, last_idx, gathered_line


def semicolons_to_newlines(src_lines):
    
    temp_lines = []
    for i in range(len(src_lines)):
        if ";" in src_lines[i]:
            if "!" in src_lines[i] and src_lines[i].find("!") < src_lines[i].rfind(";"):
                pass
            elif src_lines[i][:src_lines[i].find(";")].count('"')%2 == 1 or src_lines[i][:src_lines[i].find(";")].count("'")%2 == 1:
                pass
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


def line_length_fix(src_lines):
    
    temp_lines = []
    for i in range(len(src_lines)):
        if len(src_lines[i]) > 132 and src_lines[i].strip().startswith("!"):
            whitespace, line_text = tuple(src_lines[i].split("!", 1))
            j = 0
            while j < len(line_text):
                temp_lines.append(whitespace + "!" + line_text[j:min(j+100, len(line_text) + 1)].rstrip() + "\n")
                j += 100
        else:
            temp_lines.append(src_lines[i])

    return temp_lines


def split_single_line_if_stmts(src_lines):

    i = -1
    temp_lines = []
    while i+1 < len(src_lines):
        i += 1

        # process any if statements found
        if_stmt_text = ""
        if not src_lines[i].strip().startswith("!") and (src_lines[i].lower().strip().startswith("if ") or src_lines[i].lower().strip().startswith("if(")):

            _, _, if_stmt_text = gather_statement_text(src_lines, i) 

            # iterate through the original if statement text and build up a multiline version
            j = 0
            split_if_stmt = ""
            while j < len(if_stmt_text):
                split_if_stmt += if_stmt_text[j]
                j += 1

                # we've found the condition
                if split_if_stmt[-1] == "(":
                    
                    # find the end of the condition
                    while j < len(if_stmt_text):
                        split_if_stmt += if_stmt_text[j]
                        j+=1
                        if split_if_stmt.count("(") == split_if_stmt.count(")"):
                            break

                    # if what follows the condition doesn't start with "then"
                    # this is a single line if statement and we must split it up!
                    if not re.match(r"[&\s\)]*then", if_stmt_text[j:].strip().lower(), re.IGNORECASE):
                        if if_stmt_text[j:].lstrip() and if_stmt_text[j:].lstrip()[0] in ["*", "&", "$"]:
                            split_if_stmt += " then\n      " + if_stmt_text[j:].lstrip()[1:] + "\n      endif\n"                                    
                        else:
                            split_if_stmt += " then\n      " + if_stmt_text[j:] + "\n      endif\n"
                        if_stmt_text = split_if_stmt
                    else:
                        if_stmt_text = ""
                    break
        
        if if_stmt_text:
            temp_lines.append(if_stmt_text)
        else:
            temp_lines.append(src_lines[i])

    return temp_lines


def preprocess_project(src_regex="*.f*"):

    print("** Preprocessing DISABLED")

    # result = subprocess.run(
    #     f'find -iname "{src_regex}" | grep -v "prose" | grep -v "/__" | xargs grep -Rie "^\s*module\s"',
    #     stdout=subprocess.PIPE,
    #     text=True,
    #     shell=True
    # )

    # for line in result.stdout.split("\n"):
    #     if line and ":" in line:
    #         matched_text = ":".join(line.lower().split(":")[1:]).split()
    #         MODULE_NAME_LIST.add(matched_text[matched_text.index("module") + 1])

    # result = subprocess.run(
    #     f'find -iname "{src_regex}" | grep -v "prose" | grep -v "/__"',
    #     stdout=subprocess.PIPE,
    #     text=True,
    #     shell=True
    # )
    # for line in result.stdout.split("\n"):
    #     if line:
    #         print(line.strip())
    #         preprocess_src_file(line.strip())


def preprocess_src_file(src_file):
    with open(src_file, "r") as f:
        src_lines = f.readlines()

    src_lines = semicolons_to_newlines(src_lines)
    src_lines = split_single_line_if_stmts(src_lines)
    src_lines = line_length_fix(src_lines)

    statements_to_replace = {}
    statements_to_insert = {}
    statement_idx_to_ignore = []

    i = -1
    while i+1 < len(src_lines):
        i += 1

        if module_begin_re.search(src_lines[i]):
            i = preprocess_module(src_lines, i, statements_to_replace, statements_to_insert, statement_idx_to_ignore)

    with open(src_file, "w") as f:
        i = -1
        while i+1 < len(src_lines):
            i += 1

            if i in statements_to_replace:
                f.write(statements_to_replace[i])
            
            elif i in statements_to_insert:
                # f.write(f"      !BEGIN_PROSE_INSTRUMENTATION\n")

                # for statement in statements_to_insert[i]:
                #     f.write(f"      {statement}\n")
                
                # f.write(f"      !END_PROSE_INSTRUMENTATION\n")
                f.write(src_lines[i])

            elif i not in statement_idx_to_ignore:
                f.write(src_lines[i])


def preprocess_module(src_lines, start_idx, statements_to_replace, statements_to_insert, statement_idx_to_ignore,):

    this_module_name = src_lines[start_idx].split()[1].lower()

    module_variable_declarations_to_insert = set()
    imported_modules_profiling_code = set()

    # Keep track of the index where the modules variable declarations end; we will insert all of our instrumentation after this point.
    # Process all variable declarations found along the way
    temp_statements_to_insert = []
    instrumentation_insert_idx = start_idx + 1
    procedure_insert_idx = -1
    i = start_idx

    while i+1 < len(src_lines):
        i += 1

        # skip interface blocks
        if interface_begin_re.search(src_lines[i]):
            while i+1 < len(src_lines) and not interface_end_re.search(src_lines[i]):
                i+=1
    
        # if we encountered another module, we are done with the current module
        elif module_begin_re.search(src_lines[i]):
            break

        # for every use statement involving a module in MODULE_NAME_LIST, add a statement which
        # will call the prose_profile_module_{module_name}_vars procedure in that module
        elif use_import_re.search(src_lines[i]):
            i = preprocess_use_statement(src_lines, i, imported_modules_profiling_code)
            instrumentation_insert_idx = i + 1

        # handle derived type blocks
        elif derived_type_begin_re.search(src_lines[i]):
            
            # extract derived type name 
            derived_type_name = src_lines[i][derived_type_begin_re.search(src_lines[i]).span()[-1]:].split()[0]
            while i+1 < len(src_lines) and not derived_type_end_re.search(src_lines[i]):
                i += 1

                # if we found a variable declaration...
                if variable_declaration_re.search(src_lines[i]):
                    i = preprocess_module_variable_declaration(src_lines, i, statements_to_replace, temp_statements_to_insert, statement_idx_to_ignore, derived_type_name)

            instrumentation_insert_idx = i + 1

        # make note of the location of the contains statement; we want to insert procedures after this
        elif src_lines[i].strip().lower().startswith("contains"):
            if procedure_insert_idx < 0:
                procedure_insert_idx = i + 1

        # make note of where implicit none may be; we want to insert variables after this
        elif src_lines[i].lower().strip().startswith("implicit "):
            instrumentation_insert_idx = i + 1

        # parse procedures, add instrumentation, reason about what module variables need to be declared
        elif procedure_begin_re.search(src_lines[i]) and not src_lines[i].lower().strip().startswith(("!", "end ")):
            if procedure_insert_idx < 0:
                procedure_insert_idx = i
            i = preprocess_procedure(src_lines, i, statements_to_replace, statements_to_insert, statement_idx_to_ignore, module_variable_declarations_to_insert, this_module_name, imported_modules_profiling_code)

        # manual fix for MOM6, a sneaky variable declaration from an include statement
        elif "version_variable.h" in src_lines[i]:
            instrumentation_insert_idx = i + 1

        # parse module variable declarations
        elif variable_declaration_re.search(src_lines[i]):
            i = preprocess_module_variable_declaration(src_lines, i, statements_to_replace, temp_statements_to_insert, statement_idx_to_ignore)
            instrumentation_insert_idx = i + 1

    derived_type_names = set()
    for statement in temp_statements_to_insert:
        if "%" in statement:
            derived_type_names.add(statement[statement.find("(prose_")+len("(prose_"):statement.find("%")])
    temp_statements_to_insert = [f"type({derived_type_name}) :: prose_{derived_type_name}" for derived_type_name in derived_type_names] + temp_statements_to_insert

    temp_statements_to_insert = [
        f"subroutine prose_profile_module_{this_module_name}_vars()",
        f"\tinteger :: prose_dim_counter",
    ] + [f"\t{statement}" for statement in temp_statements_to_insert if "::" in statement] + [
        f"\tif (prose_module_{this_module_name}_vars_profiled .eq. 0) then",
        f"\t\tprose_module_{this_module_name}_vars_profiled = 1",
        f"\t\tcall prose_profile_{this_module_name}_imported_vars()",
    ] + [f"\t\t{statement}" for statement in temp_statements_to_insert if "::" not in statement] + [
        f"\tendif",
        f"end subroutine",
    ]   

    # construct the procedure to call which contains calls to
    # all the prose_profile_module_{module_name}_vars procedure in the imported modules
    temp_statements_to_insert += [
        f"subroutine prose_profile_{this_module_name}_imported_vars()",
    ] + [f"\t{x}" for x in reversed(sorted(imported_modules_profiling_code))] + [
        "end subroutine"
    ]
    
    statements_to_insert[procedure_insert_idx] = temp_statements_to_insert

    # insert profiling counter module variable declarations
    profiling_counter_variable_declarations = []
    if len(temp_statements_to_insert) > 0:
        profiling_counter_variable_declarations.append(f"integer :: prose_module_{this_module_name}_vars_profiled = 0")
    profiling_counter_variable_declarations += module_variable_declarations_to_insert
    statements_to_insert[instrumentation_insert_idx] = profiling_counter_variable_declarations

    return i - 1


def preprocess_procedure(src_lines, start_idx, statements_to_replace, statements_to_insert, statement_idx_to_ignore, module_variable_declarations_to_insert, this_module_name, imported_modules_profiling_code):

    this_procedure_name = src_lines[start_idx][procedure_begin_re.search(src_lines[start_idx]).span()[-1]:]
    this_procedure_name = this_procedure_name[:this_procedure_name.find("(")].strip()

    # Keep track of the index where the procedure's variable declarations end; we will insert all of our instrumentation after this point.
    # Process all variable declarations found along the way
    temp_statements_to_insert = []
    instrumentation_insert_idx = start_idx + 1
    i = start_idx
    while i+1 < len(src_lines):
        i += 1

        # if we found another procedure declaration, we are done
        if procedure_begin_re.search(src_lines[i]) or module_begin_re.search(src_lines[i]):
            break

        # for every use statement involving a module in MODULE_NAME_LIST, add a statement which
        # will call the prose_profile_module_{module_name}_vars procedure in that module
        elif use_import_re.search(src_lines[i]):
            i = preprocess_use_statement(src_lines, i, imported_modules_profiling_code)
            instrumentation_insert_idx = i + 1
            
        # make note of where implicit none may be; we want to insert variables after this
        elif src_lines[i].lower().strip().startswith("implicit "):
            instrumentation_insert_idx = i + 1

        # manual fix for MOM6, a sneaky variable declaration from an include statement
        elif "version_variable.h" in src_lines[i]:
            instrumentation_insert_idx = i + 1

        # process variable declarations
        elif variable_declaration_re.search(src_lines[i]):
            i = preprocess_procedure_variable_declaration(src_lines, i, statements_to_replace, temp_statements_to_insert, statement_idx_to_ignore, this_procedure_name, this_module_name, module_variable_declarations_to_insert)
            instrumentation_insert_idx = i + 1
        
    if len(temp_statements_to_insert) > 0:
        temp_statements_to_insert += [
            f"if ( prose_module_{this_module_name}_vars_profiled .eq. 0 ) then",
            f"\tcall prose_profile_module_{this_module_name}_vars()",
            f"endif"
        ]
        statements_to_insert[instrumentation_insert_idx] = temp_statements_to_insert

    return i - 1


def preprocess_procedure_variable_declaration(src_lines, start_idx, statements_to_replace, temp_statements_to_insert, statement_idx_to_ignore, this_procedure_name, this_module_name, module_variable_declarations_to_insert):

    start_idx, end_idx, var_decl_text = gather_statement_text(src_lines, start_idx)

    var_instrumentation = []
    optional_var_instrumentation = []

    # process declarations of real arrays
    fp_array_var_decl_match = fp_array_var_declaration_re.match(var_decl_text)
    if fp_array_var_decl_match:

        # extract the variable names
        assignments = []
        var_names = []
        dim_infos = []

        if "::" in var_decl_text:
            j = var_decl_text.find("::") + len("::")
        else:
            j = 0

        while j < len(var_decl_text):
            var_name = ""
            assignment = ""
            dim_info = ""

            # gather variable name
            while (j < len(var_decl_text)) and not (re.match(r"[0-9a-z_]", var_decl_text[j], re.IGNORECASE)):
                j += 1
            while (j < len(var_decl_text)) and (re.match(r"[0-9a-z_]", var_decl_text[j], re.IGNORECASE)):
                var_name += var_decl_text[j]
                j += 1

            # skip white space
            while j < len(var_decl_text) and var_decl_text[j] == " ":
                j += 1

            # gather dimension info
            # case 1: variable declared with the "dimension" keyword
            if "dimension" in var_decl_text.lower():

                # extract dimension info from the dimension attribute
                dim_tokens = re.split(r"([,\)\(])", var_decl_text[var_decl_text.lower().find("dimension") + len("dimension"):var_decl_text.find("::")])
                dim_info = dim_tokens[0]
                k = 0
                while dim_info.count("(") == 0 or (dim_info.count("(") != dim_info.count(")")):
                    k += 1
                    dim_info += dim_tokens[k]

            # case 2: variable declared with the dimension given in parentheses
            else:
                # special case of end of var_decl_text
                if j == len(var_decl_text):
                    if var_name.strip():
                        var_names.append(var_name)
                        dim_infos.append("non-array")
                        assignments.append("")
                    break

                # gather dimension if it has one
                elif var_decl_text[j] == "(":
                    dim_info += var_decl_text[j]
                    while dim_info.count("(") != dim_info.count(")"):
                        j += 1
                        dim_info += var_decl_text[j]
                    j += 1

            # skip white space
            while (j < len(var_decl_text)) and var_decl_text[j] == " ":
                j += 1
            
            # gather initialization if it has one
            if (j < len(var_decl_text)) and var_decl_text[j] == "=":
                assignment += "="
                j += 1

                # save pointer assignment initializations
                if var_decl_text[j] == ">":

                    assignment += ">"
                    j += 1

                    # gather variable name
                    while (j < len(var_decl_text)) and not (re.match(r"[\.%0-9a-z_\(\)]", var_decl_text[j], re.IGNORECASE)):
                        j += 1
                    while (j < len(var_decl_text)) and (re.match(r"[\.%0-9a-z_\(\)]", var_decl_text[j], re.IGNORECASE)):
                        assignment += var_decl_text[j]
                        j += 1

                    if assignment.endswith("NULL") or assignment.endswith("null"):
                        assignment += "()"

                # save array literal assignment initializations
                elif var_decl_text[j:].strip().startswith("(/"):
                    while (j < len(var_decl_text)) and (not assignment.endswith("/)")):
                        assignment += var_decl_text[j]
                        j += 1
                elif var_decl_text[j:].strip().startswith("["):
                    while (j < len(var_decl_text)) and (not assignment.endswith("]")):
                        assignment += var_decl_text[j]
                        j += 1

                # save scalar assignment initializations
                else:
                    while (j < len(var_decl_text)) and not (re.match(r"[\.%0-9a-z_\(\)]", var_decl_text[j], re.IGNORECASE)):
                        j += 1
                    while (j < len(var_decl_text)) and (re.match(r"[\.%0-9a-z_\(\)]", var_decl_text[j], re.IGNORECASE)):
                        assignment += var_decl_text[j]
                        j += 1

            # check to be sure the variable name is not one of the fortran declaration modifiers
            if var_name.lower() in FORTRAN_DECLARATION_MODIFIERS:
                continue
            else:
                var_names.append(var_name)
                assignments.append(assignment)
            if dim_info:
                dim_infos.append(dim_info[1:-1])
            else:
                dim_infos.append("non-array")
            j += 1

        # split the declaration into multiple lines if necessary
        if len(var_names) > 0:
            replacement_statements = ""

            if "::" in var_decl_text:
                prefix = var_decl_text[:var_decl_text.find("::")]                    
            else:
                prefix = var_decl_text[:re.search(r"\s"+var_names[0], var_decl_text).start()]

            for j in range(len(var_names)):
                if "dimension" in prefix.lower() or dim_infos[j] == "non-array":
                    replacement_statements += f"{prefix}:: {var_names[j]} {assignments[j]}\n"
                else:
                    replacement_statements += f"{prefix}, dimension({dim_infos[j]}):: {var_names[j]} {assignments[j]}\n"
            
            statements_to_replace[start_idx] = replacement_statements
            
            # when splitting a multi-line, multi-variable declaration, we need to ignore the extra lines of 
            # the original declaration when printing out the instrumented program
            for idx in range(start_idx + 1, end_idx + 1):
                statement_idx_to_ignore.append(idx)

        # construct instrumentation
        for j in range(len(var_names)):

            var_name = var_names[j]
            dim_info = dim_infos[j]

            # if the variable is a scalar variable in the same declaration as an array, skip
            if (dim_info == "non-array"):
                continue

            if optional_variable_declaration_re.search(var_decl_text):
                optional_var_instrumentation.append(f"if ( present({var_name}) ) then") # execution count for this line is # of proc calls
                optional_var_instrumentation.append(f"\tif ( prose_optional_{var_name}_var_profiled .eq. 0 ) then") # execution count for this line is # of times optional var is present
                optional_var_instrumentation.append(f"\t\tprose_optional_{var_name}_var_profiled = 1")
                optional_var_instrumentation.append(f"\t\tdo prose_dim_counter=1,size({var_name});end do")
                optional_var_instrumentation.append(f"\tendif")
                optional_var_instrumentation.append(f"endif")
                module_variable_declarations_to_insert.add(f"integer :: prose_optional_{var_name}_var_profiled = 0")
            else:
                var_instrumentation.append(f"\tdo prose_dim_counter=1,size({var_name});end do")
                module_variable_declarations_to_insert.add(f"integer :: prose_procedure_{this_procedure_name}_vars_profiled = 0")

    if len(var_instrumentation) + len(optional_var_instrumentation) > 0:
        
        # remove pure keyword for profiling but save the original header in a comment
        if re.search(r"\s?pure\s", src_lines[start_idx], re.IGNORECASE):
            statements_to_replace[start_idx] = src_lines[start_idx].lower().replace("pure", "")
            temp_statements_to_insert += ["!" + src_lines[start_idx]]
        
        if "integer :: prose_dim_counter" not in temp_statements_to_insert:
            temp_statements_to_insert += ["integer :: prose_dim_counter"]
        if len(var_instrumentation) > 0:
            temp_statements_to_insert += [
                f"if ( prose_procedure_{this_procedure_name}_vars_profiled .eq. 0 ) then",
                f"\tprose_procedure_{this_procedure_name}_vars_profiled = 1"
            ] + var_instrumentation + [
                "endif"
            ]
        if len(optional_var_instrumentation) > 0:
            temp_statements_to_insert += optional_var_instrumentation

    return end_idx


def preprocess_module_variable_declaration(src_lines, start_idx, statements_to_replace, temp_statements_to_insert, statement_idx_to_ignore, derived_type_name=""):
    start_idx, end_idx, var_decl_text = gather_statement_text(src_lines, start_idx)

    # process declarations of real arrays
    fp_array_var_decl_match = fp_array_var_declaration_re.match(var_decl_text)
    if fp_array_var_decl_match:

        # extract the variable names
        assignments = []
        var_names = []
        dim_infos = []

        if "::" in var_decl_text:
            j = var_decl_text.find("::") + len("::")
        else:
            j = 0

        while j < len(var_decl_text):
            var_name = ""
            assignment = ""
            dim_info = ""

            # gather variable name
            while (j < len(var_decl_text)) and not (re.match(r"[0-9a-z_]", var_decl_text[j], re.IGNORECASE)):
                j += 1
            while (j < len(var_decl_text)) and (re.match(r"[0-9a-z_]", var_decl_text[j], re.IGNORECASE)):
                var_name += var_decl_text[j]
                j += 1

            # skip white space
            while j < len(var_decl_text) and var_decl_text[j] == " ":
                j += 1

            # gather dimension info
            # case 1: variable declared with the "dimension" keyword
            if "dimension" in var_decl_text.lower():

                # extract dimension info from the dimension attribute
                dim_tokens = re.split(r"([,\)\(])", var_decl_text[var_decl_text.lower().find("dimension") + len("dimension"):var_decl_text.find("::")])
                dim_info = dim_tokens[0]
                k = 0
                while dim_info.count("(") == 0 or (dim_info.count("(") != dim_info.count(")")):
                    k += 1
                    dim_info += dim_tokens[k]

            # case 2: variable declared with the dimension given in parentheses
            else:
                # special case of end of var_decl_text
                if j == len(var_decl_text):
                    if var_name.strip():
                        var_names.append(var_name)
                        dim_infos.append("non-array")
                        assignments.append("")
                    break

                # gather dimension if it has one
                elif var_decl_text[j] == "(":
                    dim_info += var_decl_text[j]
                    while dim_info.count("(") != dim_info.count(")"):
                        j += 1
                        dim_info += var_decl_text[j]
                    j += 1

            # skip white space
            while (j < len(var_decl_text)) and var_decl_text[j] == " ":
                j += 1
            
            # gather initialization if it has one
            if (j < len(var_decl_text)) and var_decl_text[j] == "=":
                assignment += "="
                j += 1

                # save pointer assignment initializations
                if var_decl_text[j] == ">":

                    assignment += ">"
                    j += 1

                    # gather variable name
                    while (j < len(var_decl_text)) and not (re.match(r"[\.%0-9a-z_\(\)]", var_decl_text[j], re.IGNORECASE)):
                        j += 1
                    while (j < len(var_decl_text)) and (re.match(r"[\.%0-9a-z_\(\)]", var_decl_text[j], re.IGNORECASE)):
                        assignment += var_decl_text[j]
                        j += 1

                    if assignment.endswith("NULL") or assignment.endswith("null"):
                        assignment += "()"

                # save array literal assignment initializations
                elif var_decl_text[j:].strip().startswith("(/"):
                    while (j < len(var_decl_text)) and (not assignment.endswith("/)")):
                        assignment += var_decl_text[j]
                        j += 1

                # save scalar assignment initializations
                else:
                    while (j < len(var_decl_text)) and not (re.match(r"[\.%0-9a-z_\(\)]", var_decl_text[j], re.IGNORECASE)):
                        j += 1
                    while (j < len(var_decl_text)) and (re.match(r"[\.%0-9a-z_\(\)]", var_decl_text[j], re.IGNORECASE)):
                        assignment += var_decl_text[j]
                        j += 1

            # check to be sure the variable name is not one of the fortran declaration modifiers
            if var_name.lower() in FORTRAN_DECLARATION_MODIFIERS:
                continue
            else:
                var_names.append(var_name)
                assignments.append(assignment)
            if dim_info:
                dim_infos.append(dim_info[1:-1])
            else:
                dim_infos.append("non-array")
            j += 1

        # split the declaration into multiple lines if necessary
        if len(var_names) > 0:
            replacement_statements = ""

            if "::" in var_decl_text:
                prefix = var_decl_text[:var_decl_text.find("::")]                    
            else:
                prefix = var_decl_text[:re.search(r"\s"+var_names[0], var_decl_text).start()]

            for j in range(len(var_names)):
                if "dimension" in prefix.lower() or dim_infos[j] == "non-array":
                    replacement_statements += f"{prefix} :: {var_names[j]} {assignments[j]}\n"
                else:
                    replacement_statements += f"{prefix}, dimension({dim_infos[j]}) :: {var_names[j]} {assignments[j]}\n"
            
            statements_to_replace[start_idx] = replacement_statements
            
            # when splitting a multi-line, multi-variable declaration, we need to ignore the extra lines of 
            # the original declaration when printing out the instrumented program
            for idx in range(start_idx + 1, end_idx + 1):
                statement_idx_to_ignore.append(idx)

        # construct instrumentation
        for j in range(len(var_names)):

            var_name = var_names[j]
            dim_info = dim_infos[j]

            if dim_info != "non-array":
                if derived_type_name:       
                    temp_statements_to_insert.append(f"do prose_dim_counter=1,size(prose_{derived_type_name}%{var_name});end do")
                else:
                    temp_statements_to_insert.append(f"do prose_dim_counter=1,size({var_name});end do")

    return end_idx


def preprocess_use_statement(src_lines, start_idx, imported_modules_profiling_code):
        
    # gather the text to parse
    _, end_idx, import_text = gather_statement_text(src_lines, start_idx)
    
    # "use-only" handling; need to explicitly import the profiling function
    if use_only_import_re.search(import_text):
        use_only_prefix = import_text[:import_text.find(":")+1].strip()
        module_name = use_only_prefix.split()[1].replace(",", "")
        if module_name.lower() in MODULE_NAME_LIST:
            imported_modules_profiling_code.add(f"{use_only_prefix} prose_profile_module_{module_name}_vars")
            imported_modules_profiling_code.add(f"call prose_profile_module_{module_name}_vars()")

    # "use" handling; possibly many modules in a single "use" statement
    else:
        for module_name in [x.strip() for x in import_text[use_import_re.search(import_text).span()[-1]:].split(",")]:
            if module_name.lower() in MODULE_NAME_LIST:
                imported_modules_profiling_code.add(f"call prose_profile_module_{module_name}_vars()")

    return end_idx


def postprocess_project(project_root):

    # assume we used gcov...
    find_gcno_dir_command = 'find -name "*.gcno" | head -n 1 | xargs dirname | xargs realpath'
    result = subprocess.run(
        find_gcno_dir_command,
        cwd=project_root,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True
    )
    gcno_dir_path = result.stdout.decode("utf-8").strip()

    # ...but if we can't find the gcno files, assume we used intel's codecov...
    if "dirname: missing operand" in gcno_dir_path:
        subprocess.run(["profmerge", "-prof_dir", "prose_workspace/__profiling"])
        subprocess.run(["codecov", "-txtlcov", "-spi", "prose_workspace/__profiling/pgopti.spi", "-dpi", "prose_workspace/__profiling/pgopti.dpi"])
        subprocess.run(["mv", "CodeCoverage", "prose_workspace/__profiling/code_coverage"])
        subprocess.run(["mv", "CODE_COVERAGE.TXT", "prose_workspace/__profiling"])

    # ...otherwise, proceed
    else:
        subprocess.run("mv {}/*.gc* prose_workspace/__profiling".format(gcno_dir_path), shell=True)
        os.mkdir("prose_workspace/__profiling/code_coverage")

        for gcno_file in [os.path.join("prose_workspace/__profiling", filename) for filename in os.listdir("prose_workspace/__profiling") if filename.endswith(".gcno")]:
            result = subprocess.run(
                ["gcov", gcno_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

            stdout = result.stdout.decode("utf-8").split()

            i = 0
            while i < len(stdout) and stdout[i] != 'File':
                i += 1

            if i < len(stdout) and stdout[i] == 'File':
                source_file_path = stdout[i+1].replace("'", "")
                gcov_file = stdout[-1].replace("'", "")

                source_file_path = source_file_path.replace("/", "_")
                source_file_path = source_file_path.replace(".", "_")
                source_file_path = source_file_path.upper()

                os.rename(gcov_file, os.path.join("prose_workspace/__profiling/code_coverage", source_file_path + ".LCOV"))
    
    # process each source file's coverage info
    LCOV_file_paths = [os.path.join("prose_workspace/__profiling/code_coverage",x) for x in os.listdir("prose_workspace/__profiling/code_coverage") if x.endswith(".LCOV")]
    for LCOV_file_path in LCOV_file_paths:

        # parse human-readable LCOV file for original src path and line execution counts
        with open(LCOV_file_path, "r") as f:
            LCOV_lines = f.readlines()
        original_src_path = LCOV_lines[0][LCOV_lines[0].find("/"):].strip()
        execution_counts = []
        for line in LCOV_lines:
            if ":" in line:

                weight = -1
                line_no = -1
                line_collection = line.split(':')

                if line_collection[0].strip() == "#####":
                    weight = 0

                try:
                    line_no = int(line_collection[1].replace(',',''))
                    weight = float(line_collection[0].replace(',','').replace("*",''))
                except ValueError:
                    pass

                if line_no > 0:
                    execution_counts.append(weight)
        
        # read in src code
        with open(original_src_path, "r") as f:
            src_lines = f.readlines()

        src_idx_to_dims = {}
        i = -1
        while i+1 < len(src_lines):
            i += 1
            if module_begin_re.search(src_lines[i]):
                i = postprocess_module(src_lines, i, src_idx_to_dims, execution_counts)

        # save copy of original version of file without added instrumentation
        with open("prose_workspace/__original_files/{}".format(os.path.basename(original_src_path)), "w") as f:
            i = -1
            while i+1 < len(src_lines):
                i += 1
                if "!BEGIN_PROSE_INSTRUMENTATION" in src_lines[i]:
                    while "!END_PROSE_INSTRUMENTATION" not in src_lines[i]:
                        i += 1
                else:
                    f.write(src_lines[i])

        # modify original source file to remove instrumentation and add comments with line numbers
        with open(original_src_path, "w") as f:
            i = -1
            while i+1 < len(src_lines):
                i += 1
                if "!BEGIN_PROSE_INSTRUMENTATION" in src_lines[i]:
                    while "!END_PROSE_INSTRUMENTATION" not in src_lines[i]:
                        i += 1
                else:
                    if src_lines[i].strip().startswith("#"):
                        out_line = src_lines[i]
                    else:
                        out_line = src_lines[i].lower()
                        
                    f.write(f"!PROSE_{i}\n" + out_line)
                    while src_lines[i].strip().endswith("&"):
                        i += 1
                        f.write(src_lines[i].lower())

        # write out Prose-parseable .bcov file with the addition of the calculated array dimensions
        bcov_file_path = LCOV_file_path[:LCOV_file_path.rfind(".")] + ".bcov"
        with open(bcov_file_path, "w") as f:      
            for i in range(len(src_lines)):
                if i in src_idx_to_dims:
                    f.write(f"{i}:{src_idx_to_dims[i]}:{src_lines[i].strip().lower()}\n")
                else:
                    f.write(f"{i}:{execution_counts[i]}:{src_lines[i].strip().lower()}\n")


def postprocess_module(src_lines, start_idx, src_idx_to_dims, execution_counts):
    
    module_array_names_to_dims = {}
    end_idx = len(src_lines) - 1

    # first pass to parse procedures whose execution counts contain array dimension information
    i = start_idx
    while i+1 < len(src_lines):
        i += 1

        # skip interface blocks
        if interface_begin_re.search(src_lines[i]):
            while i+1 < len(src_lines) and not interface_end_re.search(src_lines[i]):
                i += 1

        elif module_begin_re.search(src_lines[i]):
            end_idx = i - 1 
            break

        # find each procedure's array type variables and associate them with the appropriate
        # dimension taken from the execution counts
        elif procedure_begin_re.search(src_lines[i]) and not src_lines[i].lower().strip().startswith(("!", "end ")):
            i = postprocess_procedure(src_lines, i, src_idx_to_dims, module_array_names_to_dims, execution_counts)

    # second pass to associate discovered array dimensions with module variable declarations
    i = start_idx
    while i+1 < len(src_lines):
        i += 1

        # skip interface blocks
        if interface_begin_re.search(src_lines[i]):
            while i+1 < len(src_lines) and not interface_end_re.search(src_lines[i]):
                i += 1

        elif module_begin_re.search(src_lines[i]) or procedure_begin_re.search(src_lines[i]):
            break

        # handle derived type blocks
        elif derived_type_begin_re.search(src_lines[i]):
            
            # extract derived type name 
            derived_type_name = src_lines[i][derived_type_begin_re.search(src_lines[i]).span()[-1]:].split()[0]
            while i+1 < len(src_lines) and not derived_type_end_re.search(src_lines[i]):
                i += 1

                # if we found a variable declaration...
                if variable_declaration_re.search(src_lines[i]):
                    i = postprocess_variable_declaration(src_lines, i, src_idx_to_dims, module_array_names_to_dims, derived_type_name)

        # handle module variable declarations
        elif variable_declaration_re.search(src_lines[i]):
            i = postprocess_variable_declaration(src_lines, i, src_idx_to_dims, module_array_names_to_dims)

    return end_idx


def postprocess_procedure(src_lines, start_idx, src_idx_to_dims, module_array_names_to_dims, execution_counts):
    
    # find start of prose instrumentation for this procedure
    # if the prose instrumentation is the line before prodecure header, this means the
    # whole procedure is added instrumentation and this is a procedure that may contain 
    # profiling information for module variables
    if "!BEGIN_PROSE_INSTRUMENTATION" in src_lines[start_idx - 1]:
        prose_instrumentation_start_idx = start_idx
    else:
        prose_instrumentation_start_idx = -1
        i = start_idx + 1
        while i < len(src_lines):
            if "!BEGIN_PROSE_INSTRUMENTATION" in src_lines[i]:
                prose_instrumentation_start_idx = i
                break
            elif procedure_begin_re.search(src_lines[i]):
                break
            else:
                i += 1
    
        # if there is no instrumentation, move on
        if prose_instrumentation_start_idx < 0:
            return i - 1

    # parse added instrumentation
    procedure_array_names_to_dims = {}
    i = -1
    while "!END_PROSE_INSTRUMENTATION" not in src_lines[prose_instrumentation_start_idx + i + 1]:
        i += 1

        # special case: if there is a commented line here containing the keyword "pure"
        # then we removed it from this procedure's declaration because that hampered our profiling
        # strategy. We add it back now
        if re.search(r"^\s*!(.*\s)?\s?pure\s", src_lines[prose_instrumentation_start_idx + i], re.IGNORECASE):
            src_lines[start_idx] = src_lines[prose_instrumentation_start_idx + i][src_lines[prose_instrumentation_start_idx + i].find("!"):]
        
        # execution counts for these lines will represent the array dimensions we want;
        elif "do prose_dim_counter=1,size(" in src_lines[prose_instrumentation_start_idx + i]:
            array_var_name = src_lines[prose_instrumentation_start_idx + i].split("(")[-1].split(")")[0].strip()
            array_dim = execution_counts[prose_instrumentation_start_idx + i]
            
            # if it is an optional array, we scale the dimension by the percentage of procedure calls
            # in which that array is actually present
            if "if ( present(" in src_lines[prose_instrumentation_start_idx + i - 3]:
                array_dim = array_dim * (execution_counts[prose_instrumentation_start_idx + i - 2] / execution_counts[prose_instrumentation_start_idx + i - 3])
            
            if prose_instrumentation_start_idx == start_idx:
                module_array_names_to_dims[array_var_name] = array_dim
            else:
                procedure_array_names_to_dims[array_var_name] = array_dim

    # find all the original array declarations in this procedure; we will make the array dimension
    # visible in the coverage info as an "execution count"
    i = start_idx
    while i+1 < prose_instrumentation_start_idx:
        i += 1

        if variable_declaration_re.search(src_lines[i]):
            i = postprocess_variable_declaration(src_lines, i, src_idx_to_dims, procedure_array_names_to_dims)

    # find procedure end
    i = start_idx
    while i+1 < len(src_lines):
        i += 1
        if procedure_begin_re.search(src_lines[i]) or module_begin_re.search(src_lines[i]):
            break
    
    return i - 1


def postprocess_variable_declaration(src_lines, start_idx, src_idx_to_dims, array_names_to_dims, derived_type_name=""):
    start_idx, end_idx, var_decl_text = gather_statement_text(src_lines, start_idx)

    # process declarations of real arrays
    fp_array_var_decl_match = fp_array_var_declaration_re.match(var_decl_text)
    if fp_array_var_decl_match:

        # gather variable name
        i = var_decl_text.find("::") + len("::")
        while i < len(var_decl_text):

            var_name = ""
            while (i < len(var_decl_text)) and not (re.match(r"[0-9a-z_]", var_decl_text[i], re.IGNORECASE)):
                i += 1
            while (i < len(var_decl_text)) and (re.match(r"[0-9a-z_]", var_decl_text[i], re.IGNORECASE)):
                var_name += var_decl_text[i]
                i += 1
            while i < len(var_decl_text) and var_decl_text[i] == " ":
                i += 1

            # check to be sure the variable name is not one of the fortran declaration modifiers
            if var_name.lower() in FORTRAN_DECLARATION_MODIFIERS:
                i += 1
                continue

            else:

                if derived_type_name:
                    var_name = f"{derived_type_name}%{var_name}"

                # if we have a dimension from the profiling, make note of what to insert into bcov file 
                if var_name in array_names_to_dims:
                    for idx in range(start_idx, end_idx + 1): 
                        src_idx_to_dims[idx] = array_names_to_dims[var_name]
                    break

    return end_idx

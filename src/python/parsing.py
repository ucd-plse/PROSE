import re
import subprocess
import os

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
    "complex",
    "procedure",
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

fp_literal_re1 = re.compile(r"(?<![a-z0-9_])-?([0-9]+\.[0-9]*|[0-9]+\.[0-9]+)([ed]-?[0-9]*)?(_[a-z0-9_]+)?", re.IGNORECASE)
fp_literal_re2 = re.compile(r"(?<![a-z0-9_])-?([0-9]+)([ed]-?[0-9]+)(_[a-z0-9_]+)?", re.IGNORECASE)
fp_literal_nokind_re1 = re.compile(r"(?<![a-z0-9_])-?(([0-9]+\.[0-9]*)|([0-9]*\.[0-9]+))([ed]-?[0-9]*)?", re.IGNORECASE)
fp_literal_nokind_re2 = re.compile(r"(?<![a-z0-9_])-?([0-9]+)([ed]-?[0-9]+)", re.IGNORECASE)
integer_literal_nokind_re = re.compile(r"(?<![a-z0-9_])-?[0-9]+(?![\.de])")
derived_type_begin_re = re.compile(r"^\s*type(?![a-z0-9_])(?!\s*\()(.*::)?\s*([a-z_][a-z0-9_]*)", re.IGNORECASE)
derived_type_var_decl_re = re.compile(r"^\s*type\s*\(\s*([a-z_][a-z0-9_]*)", re.IGNORECASE)
use_import_re = re.compile(r"^\s*use(\s+|.*::\s*)([a-z_][a-z0-9_]*)", re.IGNORECASE)
procedure_begin_re = re.compile(r"(^|\s)(subroutine|function)\s+([a-z_][a-z0-9_]*)", re.IGNORECASE)
procedure_inside_interface_re = re.compile(r"(^|\s)(subroutine|function|procedure)\s+([a-z_][a-z0-9_]*(?:\s*,\s*[a-z_][a-z0-9_]*)*)", re.IGNORECASE)
program_begin_re = re.compile(r"^\s*program\s+([a-z_][a-z0-9_]*)", re.IGNORECASE)
module_begin_re = re.compile(r"^\s*module\s+([a-z_][a-z0-9_]*)", re.IGNORECASE)
line_continuation_free_form_re = re.compile(r"^.+&\s*(!.*)?$")
line_continuation_fixed_form_re = re.compile(r"^\s\s\s\s\s\S(?!0)")
comment_line_fixed_form_re = re.compile(r"^[cd!*]", re.IGNORECASE) 
comment_line_free_form_re = re.compile(r"^\s*!")
variable_declaration_re = re.compile(r"((^\s*((real[\s(,\*])|(integer[\s(,\*])|(character[\*\s(,])|(common[\s(,])|(save[\s(,](?!\s*$))|(parameter[\s(,])|(type\s*\()|(external[\s(,])|(double\s+precision[\s(,])|(logical[\s,])))|(::))", re.IGNORECASE)
old_pointer_declaration_re = re.compile(r"^\s*pointer\s*\(\s*([a-z_][a-z0-9_]*)\s*,\s*([a-z_][a-z0-9_]*)\s*\)", re.IGNORECASE)
real_variable_declaration_re = re.compile(r"^\s*real[\s(,\*]", re.IGNORECASE)
possible_proc_call_re = re.compile(r"(call\s+)?([a-z_][a-z0-9_]*)\s*\(", re.IGNORECASE)
possible_string_literal_re = re.compile(r"[\"']")
open_paren_re = re.compile(r"\(")
valid_fortran_names_re = re.compile(r"(?<![a-z0-9_\.])([a-z_][a-z0-9_]*(\s*%\s*[a-z_][a-z0-9_]*)*)", re.IGNORECASE)
keep_line_re = re.compile(r"^\s*(?<!!)\s*(implicit|contains|public|private)", re.IGNORECASE)
cpp_directive_re = re.compile(r"^\s*#")
include_directive_re = re.compile(r"^\s*(#\s*)?include\s*[\'\"<]([a-z0-9_\./-]+)", re.IGNORECASE)
cpp_define_directive_re = re.compile(r"^\s*#\s*(if|ifn|un)?def(ine)?\s*", re.IGNORECASE)
cpp_if_directive_re = re.compile(r"^\s*#\s*if\s*", re.IGNORECASE)
cpp_endif_directive_re = re.compile(r"^\s*#\s*endif\s*", re.IGNORECASE)
end_re = re.compile(r"^\s*end((\s*$)|\s+([a-z_][a-z0-9_]*))", re.IGNORECASE)
prose_wrapper_name_re = re.compile(r"([a-zA-Z][a-zA-Z0-9_]*)_(wrapper_(id[0-9]+_)?[048x]+_to_[048x]+|wrap_[0-9]+)")
sign_intrinsic_re = re.compile(r"[,\s]sign[\s/(]", re.IGNORECASE)
public_statement_re = re.compile(r"^\s*public(\s+|\s*::\s*)[a-z_][a-z0-9_]*", re.IGNORECASE)
implicit_statement_re = re.compile(r"^\s*implicit\s", re.IGNORECASE)
contains_statement_re = re.compile(r"^\s*contains\s", re.IGNORECASE)
fixed_sign_call_arg_string_re = re.compile(r"\(REAL\((.*),8\)\),\(REAL\((.*),8\)\)\)")
nameless_interface_re = re.compile(r"^\s*(?:abstract\s+)?interface(?![a-z0-9_])", re.IGNORECASE)
interface_begin_re = re.compile(r"(^|\s)interface\s+([a-z_][a-z0-9_]*)", re.IGNORECASE)
if_statement_re = re.compile(r"^\s*if[\s(].*then$", re.IGNORECASE)
elseif_statement_re = re.compile(r"^\s*else\s*if[\s(].*then$", re.IGNORECASE)
function_return_value_re = re.compile(r"(?<![a-z0-9_])result\s*\(\s*([a-z0-9_]+)\s*\)", re.IGNORECASE)


def is_function(line):
    line = omit_string_literals(line)
    match = procedure_begin_re.search(line)
    if match:
        return match.group(2).lower() == "function"
    else:
        return False

def function_return_value(line):
    line = omit_string_literals(line)
    match = function_return_value_re.search(line)
    if match:
        return match.group(1)
    else:
        return ""

def procedures_in_interface(line):
    line = omit_string_literals(line)
    match = procedure_inside_interface_re.search(line)
    if match and not end_statement(line) and not cpp_directive(line):
        return match.group(3).split(",")
    else:
        return []

def elseif_statement(line):
    if elseif_statement_re.search(line.rstrip()):
        return True
    else:
        return False    

def multiline_if_statement(line):
    if if_statement_re.search(line.rstrip()):
        return True
    else:
        return False    

def is_fixed_form_continuation_line(line):
    line = omit_string_literals(line)
    if line_continuation_fixed_form_re.search(line):
        return True
    else:
        return False 

def cpp_if_directive(line):
    line = omit_string_literals(line)
    if cpp_if_directive_re.search(line):
        return True
    else:
        return False

def cpp_endif_directive(line):
    line = omit_string_literals(line)
    if cpp_endif_directive_re.search(line):
        return True
    else:
        return False

def cpp_directive(line):
    line = omit_string_literals(line)
    if cpp_directive_re.search(line):
        return True
    else:
        return False

def keep_line(line):
    line = omit_string_literals(line)
    if keep_line_re.search(line):
        return True
    else:
        return False

def cpp_define_directive(line):
    line = omit_string_literals(line)
    if cpp_define_directive_re.search(line):
        return True
    else:
        return False

def end_statement(line):
    line = omit_string_literals(line)
    match = end_re.search(line)
    if match and (match.group(3) == None or match.group(3) not in FORTRAN_NON_SCOPE_END_TYPES):
        return match
    else:
        return None

def program_begin(line):
    line = omit_string_literals(line)
    match = program_begin_re.search(line)
    if match:
        return match.group(1)
    else:
        return ""

def module_begin(line):
    line = omit_string_literals(line)
    match = module_begin_re.search(line)
    if match:
        return match.group(1)
    else:
        return ""

def procedure_begin(line):
    line = omit_string_literals(line)
    if not end_statement(line) and not cpp_directive(line):
        match = procedure_begin_re.search(line)
        if match:
            return match.group(3)
    return ""

def derived_type_begin(line):
    line = omit_string_literals(line)
    if not derived_type_var_decl_re.search(line):
        match = derived_type_begin_re.search(line)
        if match:
            return match.group(2)
    return ""

def interface_begin(lines, i):
    line = omit_string_literals(lines[i])
    match = interface_begin_re.search(line)
    if match:
        return match.group(2)

    match = nameless_interface_re.search(line)
    if match:
        return f"nameless_interface@line{i}"
    return ""

def contains_use_statement(line):
    line = omit_string_literals(line)
    match = use_import_re.search(line)
    if match:
        return match.group(2)
    else:
        return ""

def contains_variable_declaration(line):
    line = omit_string_literals(line)
    match = variable_declaration_re.search(line)
    if match:
        if use_import_re.search(line) or public_statement_re.search(line) or procedure_begin_re.search(line) or derived_type_begin_re.search(line):
            return None
    else:
        match = old_pointer_declaration_re.search(line)
    return match

def old_pointer_declaration(line):
    line = omit_string_literals(line)
    return old_pointer_declaration_re.search(line)

def contains_potential_procedure_calls(line):
    line = omit_string_literals(line)
    proc_calls = []
    for match in possible_proc_call_re.finditer(line):
        proc_call = match.group(0)
        i = match.end() - 1
        unbalanced_paren_count = 1
        # try-except added to handle old, multiline FORMAT statements in ADCIRC
        try:
            while unbalanced_paren_count > 0:
                i += 1
                proc_call += line[i]
                if line[i] == "(":
                    unbalanced_paren_count += 1
                elif line[i] == ")":
                    unbalanced_paren_count -= 1
            proc_calls.append(proc_call)
        except:
            continue
    return proc_calls

def is_real_variable_declaration(line):
    line = omit_string_literals(line)
    if real_variable_declaration_re.search(line):
        return True
    else:
        return False

def remove_fp_literals_from_argument_list(text):
    temp_text = omit_parentheticals(text)
    temp_text = omit_string_literals(text)
    fp_literals = []
    fp_literals += [match.group(0) for match in fp_literal_re1.finditer(temp_text) if match]
    fp_literals += [match.group(0) for match in fp_literal_re2.finditer(temp_text) if match]
    for fp_literal in fp_literals:
        text = text.replace(fp_literal, "")
    return text

def find_valid_fortran_names(text, ignore_parentheticals=False):

    valid_fortran_names = set()

    # special handling for derived-type instantiations which take the form "type(foo)"
    # we remove that text from the remaining text which will be searched
    # if we are not ignoring parentheticals, we append a field access symbol and save
    match = re.search(r"(?<![a-z0-9_])type\s*\(\s*([a-z_][a-z0-9_]*)\s*\)", text, flags=re.IGNORECASE)
    if match:
        text = text[:match.start(0)] + text[match.end(0):]
        if not ignore_parentheticals:
            valid_fortran_names.add(match.group(1) + "%")

    # optionally omit parentheticals. This is done when extracting the declared variables in
    # a variable declaration rather than the variables parameterizing that declaration
    if ignore_parentheticals:
        text = omit_parentheticals(text)

    text = omit_string_literals(text)

    # remove intent from variable declarations
    text = re.sub(r"intent\s*\(\s*[(in)|(out)]+\s*\)", " ", text, flags=re.IGNORECASE)

    # remove kind parameters from variable declarations
    text = re.sub(r"real\s*\(\s*kind", " ", text, flags=re.IGNORECASE)

    # remove comparison operators
    text = re.sub(r"\.[a-z]+\.", " ", text, flags=re.IGNORECASE)

    # extract valid Fortran names
    candidate_fortran_names = [match[0] for match in valid_fortran_names_re.findall(text)]
    for i in range(len(candidate_fortran_names)):
        if candidate_fortran_names[i] not in FORTRAN_DECLARATION_MODIFIERS:
            valid_fortran_names.add(candidate_fortran_names[i])

    return list(valid_fortran_names)

def preprocess(src_lines, SETUP, rose_preprocessing, excluded_names, fixed_form_fortran):
    
    i = -1
    while i + 1 < len(src_lines):
        i += 1

        # don't touch non-define cpp directives
        if cpp_directive_re.search(src_lines[i]) and not cpp_define_directive_re.search(src_lines[i]):

            # parse any include files found
            if rose_preprocessing:
                match = include_directive_re.search(src_lines[i])
                if match:
                    include_file_name = match.group(2)
                    src_lines = parse_includes(src_lines, i, include_file_name, SETUP, excluded_names, fixed_form_fortran)

            continue

        src_lines[i] = remove_comments(src_lines[i], fixed_form_fortran)
        if rose_preprocessing:
            src_lines[i] = to_lower(src_lines[i])
        lines = semicolons_to_newlines(src_lines[i])
        if len(lines) > 1:
            src_lines = src_lines[:i] + lines + src_lines[i+1:]

        src_lines[i] = re.sub(r"endsubroutine", "end subroutine", src_lines[i], flags=re.IGNORECASE)
        src_lines[i] = re.sub(r"endmodule", "end module", src_lines[i], flags=re.IGNORECASE)
        src_lines[i] = re.sub(r"endfunction", "end function", src_lines[i], flags=re.IGNORECASE)
        src_lines[i] = re.sub(r"endtype", "end type", src_lines[i], flags=re.IGNORECASE)
        src_lines[i] = re.sub(r"endinterface", "end interface", src_lines[i], flags=re.IGNORECASE)
        if rose_preprocessing:
            src_lines[i] = re.sub(r"character\s*\(\s*len\s*=\s*:\s*\)(\s*,\s*allocatable\s*)?", "character(len=*)", src_lines[i], flags=re.IGNORECASE)

        i, end_idx, consolidated_line = consolidate_multiline_statement(src_lines, i, fixed_form_fortran, rose_preprocessing)
        if consolidated_line != src_lines[i]:
            gathered_lines = consolidated_line.split("\n")
            while i <= end_idx:
                if len(gathered_lines) != 0:
                    src_lines[i] = gathered_lines.pop(0) + "\n"
                else:
                    src_lines[i] = ""
                i += 1
            i = end_idx     

    # add dummy last line
    src_lines.append("")

    # add gptl timing info on a second pass
    if not rose_preprocessing:
        i = -1
        while i + 1 < len(src_lines):
            i += 1
            if procedure_begin(src_lines[i]):
                match = prose_wrapper_name_re.search(src_lines[i])
                if match:
                    src_lines, i = add_gptl_timing(src_lines, i, wrapper_name=match.group(0))

    return src_lines

def add_gptl_timing(src_lines, i, wrapper_name):
    temp = src_lines[i][src_lines[i].find(wrapper_name) + len(wrapper_name) + 1:]
    temp = temp[:temp.find(")")]
    arg_names = temp.split(",")

    i += 1
    
    to_add = [
        "#ifdef GPTL\n",
        "      use gptl\n",
        "#endif\n"
    ]
    src_lines = src_lines[:i] + to_add + src_lines[i:]
    i += len(to_add)

    while not contains_variable_declaration(src_lines[i]):
        i += 1
    while contains_variable_declaration(src_lines[i]):
        i += 1
    
    to_add = [
        "#ifdef GPTL\n",
        "      integer :: gptl_ret, gptl_handle=0\n",
        f'      gptl_ret = gptlstart_handle("<REPLACE>", gptl_handle)\n'
        "#endif\n",
    ]
    src_lines = src_lines[:i] + to_add + src_lines[i:]
    i += len(to_add) - 1

    while True:
        i += 1
        try:
            wrapped_call_args = [x.split("=")[-1].replace("__temp", "") for x in src_lines[i][src_lines[i].find("(")+1:src_lines[i].find(")")].split(",")]
        except IndexError:
            continue
        if set(wrapped_call_args) == set(arg_names):
            wrapped_procedure_name = src_lines[i][:src_lines[i].find("(")].split()[-1].strip()
            j = i
            while "<REPLACE>" not in src_lines[j]:
                j -= 1
            src_lines[j] = src_lines[j].replace("<REPLACE>", wrapped_procedure_name)
            break

    to_add = [
        "#ifdef GPTL\n",
        f"      gptl_ret = gptlstop_handle('{wrapped_procedure_name}', gptl_handle)\n",
        "#endif\n",
    ]
    src_lines = src_lines[:i] + to_add + src_lines[i:]
    i += len(to_add) + 1

    to_add = [
        "#ifdef GPTL\n",
        f"      gptl_ret = gptlstart_handle('{wrapped_procedure_name}', gptl_handle)\n"
        "#endif\n",
    ]
    src_lines = src_lines[:i] + to_add + src_lines[i:]
    i += len(to_add)

    while not end_statement(src_lines[i]):
        i += 1

    to_add = [
        "#ifdef GPTL\n",
        f"      gptl_ret = gptlstop_handle('{wrapped_procedure_name}', gptl_handle)\n"
        "#endif\n",
    ]
    src_lines = src_lines[:i] + to_add + src_lines[i:]
    i += len(to_add)

    return src_lines, i

def advance_idx_to_string_literal_end(line, i):
    """given the index of a string literal delimiter, returns index of the matching
    string literal delimiter that ends the string literal"""
    delimiter = line[i]
    while i + 1 < len(line):
        i += 1 
        if line[i] == delimiter:
            if i + 1 < len(line) and line[i] == line[i + 1] == delimiter:
                i += 1
            elif line[i - 1] == "\\":
                continue
            else:
                break
    return i

def omit_string_literals(text):
    new_text = ""
    i = -1
    while i + 1 < len(text):
        i += 1
        if text[i] in ["'", '"']:
            i = advance_idx_to_string_literal_end(text, i)
        else:
            new_text += text[i]
    return new_text

def omit_parentheticals(text):
    new_text = ""
    i = -1
    unbalanced_paren_count = 0
    while i + 1 < len(text):
        i += 1
        if text[i] == "(":
            unbalanced_paren_count += 1
            while unbalanced_paren_count != 0:
                i += 1

                if text[i] == "(":
                    unbalanced_paren_count += 1
                elif text[i] == ")":
                    unbalanced_paren_count -= 1                        
        else:
            new_text += text[i]
    return new_text

def find_comment_begin_idx(line, fixed_form_fortran):

    if fixed_form_fortran:
        if comment_line_fixed_form_re.search(line):
            return 0

    i = -1
    while i + 1 < len(line):
        i += 1
        if line[i] in ['"', "'"]:
            i = advance_idx_to_string_literal_end(line, i)
        elif line[i] == "!":
            return i

    return -1

def remove_comments(line, fixed_form_fortran):
    comment_begin_idx = find_comment_begin_idx(line, fixed_form_fortran)
    if comment_begin_idx >= 0:
        return line[:comment_begin_idx] + "\n"
    else:
        return line

def find_semicolon_idxs(line):
    idxs = []
    i = -1
    while i + 1 < len(line):
        i += 1
        if line[i] in ['"', "'"]:
            i = advance_idx_to_string_literal_end(line, i)
        elif line[i] == ";":
            idxs.append(i)

    return idxs

def to_lower(line):
    if cpp_directive(line) and not cpp_define_directive(line):
        return line

    i = -1
    new_line = ""
    while i + 1 < len(line):
        i += 1
        if line[i] in ['"', "'"]:
            j = advance_idx_to_string_literal_end(line, i)
            new_line += line[i:j+1]
            i = j                
        else:
            new_line += line[i].lower()
    return new_line

def parse_includes(src_lines, i, include_file_name, SETUP, excluded_names, fixed_form_fortran):
    unique_found = False
    multiple_found = []
    for src_search_path in SETUP['machine']['src_search_paths'].split("|"):
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
            with open(depend_src_path, "r") as f:
                include_lines = f.readlines()
            include_lines = preprocess(src_lines=include_lines, SETUP=SETUP, rose_preprocessing=True, excluded_names=excluded_names, fixed_form_fortran=fixed_form_fortran)
            src_lines = src_lines[:i] + ["      !PROSE!\n"] + include_lines + ["      !PROSE!\n"] + src_lines[i+1:]
            break
        elif len(find_results) > 1:
            multiple_found += find_results
    if not unique_found and multiple_found == []:
        excluded_names.add(include_file_name)
        src_lines[i] = ""
    elif not unique_found and multiple_found != []:
        print(f"multiple matches for {include_file_name} found: {multiple_found}")
        assert(False)

    return src_lines

def semicolons_to_newlines(line):
    semicolon_idxs = find_semicolon_idxs(line)
    if len(semicolon_idxs) > 0:
        line = list(line)
        for semicolon_idx in semicolon_idxs:
            line[semicolon_idx] = "\n      "
        line = "".join(line)
    return [l + "\n" for l in line.split("\n") if l.strip()]

def consolidate_multiline_statement(src_lines, init_idx, fixed_form_fortran, rose_preprocessing):
    
    consolidated_line = ""

    # fixed form
    if fixed_form_fortran:
        if comment_line_fixed_form_re.search(src_lines[init_idx]) or comment_line_free_form_re.search(src_lines[init_idx]) or not line_continuation_fixed_form_re.search(src_lines[init_idx]):
            return init_idx, init_idx, src_lines[init_idx]
        
        # find potential statement end
        # conditional tests for any line continuation characters in column 5, blank lines, cpp directives, or comment lines
        last_idx = init_idx
        while last_idx < len(src_lines) and (line_continuation_fixed_form_re.search(src_lines[last_idx]) or src_lines[last_idx].strip() == "" or comment_line_fixed_form_re.search(src_lines[last_idx]) or comment_line_free_form_re.search(src_lines[last_idx]) or src_lines[last_idx].strip().startswith("#")):
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
            consolidated_line = remove_comments(src_lines[first_idx].strip()[1:], fixed_form_fortran).strip() + consolidated_line
            if rose_preprocessing:
                consolidated_line = to_lower(consolidated_line)

            while first_idx - 1 > 0 and (src_lines[first_idx - 1].strip() == "" or comment_line_fixed_form_re.search(src_lines[first_idx - 1]) or comment_line_free_form_re.search(src_lines[first_idx - 1]) or src_lines[first_idx - 1].strip().startswith("#")):
                first_idx -= 1
                if comment_line_fixed_form_re.search(src_lines[first_idx]) or comment_line_free_form_re.search(src_lines[first_idx]) or src_lines[first_idx].strip() == "":
                    continue
                elif src_lines[first_idx].strip().startswith("#"):
                    consolidated_line = "\n" + src_lines[first_idx].strip() + "\n     & " + consolidated_line
        
        first_idx -= 1
        consolidated_line = remove_comments(src_lines[first_idx].rstrip(), fixed_form_fortran).rstrip() + consolidated_line
        if rose_preprocessing:
            consolidated_line = to_lower(consolidated_line)


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
            trimmed_line = remove_comments(src_lines[last_idx].strip(), fixed_form_fortran).strip()
            if rose_preprocessing:
                trimmed_line = to_lower(trimmed_line)
            if trimmed_line.endswith("&"):
                trimmed_line = trimmed_line[:trimmed_line.rfind("&")]
            if trimmed_line.startswith("&"):
                trimmed_line = trimmed_line[trimmed_line.rfind("&") + 1:]
            consolidated_line = consolidated_line + trimmed_line

            while last_idx + 1 < len(src_lines) and (src_lines[last_idx + 1].strip() == "" or comment_line_free_form_re.search(src_lines[last_idx + 1]) or src_lines[last_idx + 1].strip().startswith("#")):
                last_idx += 1
                if comment_line_free_form_re.search(src_lines[last_idx]) or src_lines[last_idx].strip() == "":
                    continue
                elif src_lines[last_idx].strip().startswith("#"):
                    consolidated_line = consolidated_line + "&\n" + src_lines[last_idx].strip() + "\n"
                
        last_idx += 1
        trimmed_line = remove_comments(src_lines[last_idx].strip(), fixed_form_fortran).strip()
        if rose_preprocessing:
            trimmed_line = to_lower(trimmed_line)
        if trimmed_line.endswith("&"):
            trimmed_line = trimmed_line[:trimmed_line.rfind("&")]
        if trimmed_line.startswith("&"):
            trimmed_line = trimmed_line[trimmed_line.rfind("&") + 1:]
        consolidated_line = consolidated_line + trimmed_line

    return first_idx, last_idx, consolidated_line
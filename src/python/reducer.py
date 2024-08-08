import subprocess
import os
import parsing
import re
from shutil import copy
from glob import glob

SETUP = None
SRC_FILES = set()
TARGETED_SRC_FILES = set()
UPSTREAM_SRC_FILES = set()
DOWNSTREAM_SRC_FILES = set()
PATH_TO_SRC_FILE_MAP = {}
MODULE_NAME_TO_SRC_FILE_MAP = {}
EXCLUDED_NAMES = set()
PROCEDURE_CALL_DEPENDENCIES = set()
SEARCH_SPACE = set()

# TODO: handle statement labels in fixed form fortran

def reduce(setup_dict):
    try:
        return _reduce(setup_dict)
    except:
        import pdb, traceback, sys
        _, _, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)


def _reduce(setup_dict):

    global SETUP
    global SRC_FILES
    global TARGETED_SRC_FILES
    global UPSTREAM_SRC_FILES
    global DOWNSTREAM_SRC_FILES
    global PATH_TO_SRC_FILE_MAP
    global MODULE_NAME_TO_SRC_FILE_MAP
    global EXCLUDED_NAMES
    global PROCEDURE_CALL_DEPENDENCIES
    global SEARCH_SPACE

    SETUP = setup_dict

    # first pass
    # 1. register all the scope symbols, variable declarations, variable symbols, and procedure calls in the targeted src code
    # 2. assign "target" attribute to all targeted FP variables
    print("** Pass 1/3")

    # targeted src code
    for src_path in SETUP['target']['src_files'].split("|"):
        print(f"\t{src_path}")
        s = SourceFile(src_path)
        s.is_targeted = True
        SRC_FILES.add(s)
        TARGETED_SRC_FILES.add(s)
        PATH_TO_SRC_FILE_MAP[s.src_path] = s
        for module_name in s.module_names:
            MODULE_NAME_TO_SRC_FILE_MAP[module_name] = s
    
    with open("prose_logs/__config_template.txt", "w") as f:
        f.writelines(SEARCH_SPACE)

    # src code that is "upstream" from targeted src code
    for src_file in TARGETED_SRC_FILES:
        for src_path in src_file.upstream_src_paths:
            if src_path not in PATH_TO_SRC_FILE_MAP.keys():
                print(f"\t{src_path}")
                s = SourceFile(src_path)
                SRC_FILES.add(s)
                UPSTREAM_SRC_FILES.add(s)
                PATH_TO_SRC_FILE_MAP[s.src_path] = s
                for module_name in s.module_names:
                    MODULE_NAME_TO_SRC_FILE_MAP[module_name] = s

    # src code that is "downstream" from both "targeted" and "upstream" source code
    for src_file in set.union(TARGETED_SRC_FILES, UPSTREAM_SRC_FILES):
        for src_path in src_file.module_name_to_downstream_src_path_map.values():
            if src_path not in PATH_TO_SRC_FILE_MAP.keys():
                print(f"\t{src_path}")
                s = SourceFile(src_path)
                SRC_FILES.add(s)
                DOWNSTREAM_SRC_FILES.add(s)
                PATH_TO_SRC_FILE_MAP[s.src_path] = s
                for module_name in s.module_names:
                    MODULE_NAME_TO_SRC_FILE_MAP[module_name] = s

    # src code that is "downstream" from "downstream" source code
    fixed_point = False
    while not fixed_point:

        new_downstream = set()
        for src_file in DOWNSTREAM_SRC_FILES:
            for downstream_src_path in src_file.module_name_to_downstream_src_path_map.values():
                if downstream_src_path not in PATH_TO_SRC_FILE_MAP.keys():
                    new_downstream.add(downstream_src_path)
        
        for src_path in new_downstream:
            print(f"\t{src_path}")
            s = SourceFile(src_path)
            SRC_FILES.add(s)
            DOWNSTREAM_SRC_FILES.add(s)
            PATH_TO_SRC_FILE_MAP[s.src_path] = s
            for module_name in s.module_names:
                MODULE_NAME_TO_SRC_FILE_MAP[module_name] = s

        fixed_point = len(new_downstream) == 0

    # second pass
    # 1. resolve references
    print("** Pass 2/3")
    for s in SRC_FILES:
        s.resolve_references()
    print("\t source files containing procedures called from targeted source:")
    for src_path in PROCEDURE_CALL_DEPENDENCIES:
        print(f"\t\t{src_path}")

    # third pass
    # 1. propagate taint to all variables referenced in the same procedure call as targeted variables
    # 2. propagate taint to all variables used in the declaration of targeted or tainted variables
    # 3. propagate taint to all variables declared within the same procedure as targeted or tainted variables
    print("** Pass 3/3")
    fixed_point = False
    while not fixed_point:

        new_taint_count = 0
        for s in SRC_FILES:
            new_taint_count += s.propagate_taint()
        
        fixed_point = new_taint_count == 0

    # unparse minimal required source code, which consists of:
    # - any scope containing references to targeted or tainted variables
    # - any declarations of targeted or tainted variables
    # - any procedure call containing references to targeted or tainted variables
    # - appropriate use statements within scopes that contain references to targeted or tainted variables that are not declared within that scope
    # - all interfaces
    to_remove = []
    for s in SRC_FILES:
        if s.is_tainted:
            copy(s.src_path, os.path.join(os.environ['PROSE_EXPERIMENT_DIR'], "prose_workspace/original_files", os.path.basename(s.src_path) + ".orig"))
            s.unparse()
        else:
            to_remove.append(s)

    for s in to_remove:
        del PATH_TO_SRC_FILE_MAP[s.src_path]
        for module_name in s.module_names:
            del MODULE_NAME_TO_SRC_FILE_MAP[module_name]
        SRC_FILES.remove(s)
        del s

    # find the order in which target files should be processed
    SRC_FILES = list(SRC_FILES)
    process_order = []
    while len(SRC_FILES) > 0:
        i = 0
        before = len(process_order)
        while i < len(SRC_FILES):
            if set(SRC_FILES[i].module_name_to_downstream_src_path_map.values()).intersection(PATH_TO_SRC_FILE_MAP.keys()) <= set([s.src_path for s in process_order + [SRC_FILES[i]]]):
                process_order.append(SRC_FILES[i])
                SRC_FILES.pop(i)
            else:
                i += 1
        after = len(process_order)
        if before == after:
            print("\n\nCircular or missing dependencies detected")
            assert(False)

    for i, s in enumerate(process_order):
        print(f"\t processing file {i+1: >4}/{len(process_order): <3} {s.src_path: <120}", end='\r', flush=True)

        include_dirs = list(set([os.path.relpath(os.path.dirname(os.path.abspath(ss)), start=os.path.dirname(ss)) for ss in s.module_name_to_downstream_src_path_map.values()]))
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
            "-Drose_comp",
            include_dirs,
            f"-I{os.path.join(os.environ['PROSE_EXPERIMENT_DIR'], 'prose_workspace/rmod_files')}",
            SETUP['target']['additional_plugin_flags'],
            f"{os.path.basename(s.src_path)}",
            "&& rm *postprocessed*",
        ]

        subprocess.run(
            " ".join(command),
            shell=True,
            env=os.environ.copy(),
            check=True,
            cwd = os.path.abspath(os.path.dirname(s.src_path))
        )

        # postprocess and move any generated rmod files
        for rmod_file_name_with_path in glob(os.path.join(os.path.abspath(os.path.dirname(s.src_path)), "*.rmod")):
            rmod_file_name = os.path.basename(rmod_file_name_with_path)
            subprocess.run(
                f"operator_fixup.py {rmod_file_name} && mv {rmod_file_name} {os.path.join(os.environ['PROSE_EXPERIMENT_DIR'], 'prose_workspace/rmod_files')}",
                shell=True,
                env=os.environ.copy(),
                check=True,
                cwd = os.path.abspath(os.path.dirname(s.src_path))
            )

    return [s.src_path for s in process_order]


class Node:

    def __init__(self, parent):
        self.parent = parent
        self._is_tainted = False


    @property
    def is_tainted(self):
        return self._is_tainted
    

    @is_tainted.setter
    def is_tainted(self, value):
        self._is_tainted = value


    def get_enclosing_src_file(self):
        node = self
        while not isinstance(node, SourceFile):
            node = node.parent
        return node


    def resolve_references(self):
        """To be implemented in the subclass"""
        raise NotImplementedError


    def propagate_taint(self):
        """To be implemented in the subclass"""
        raise NotImplementedError   
    

    def unparse(self):
        """To be implemented in the subclass"""
        raise NotImplementedError   


class Scope(Node):

    def __init__(self, name, start_idx, parent):
        super().__init__(parent)
        self.name = name
        self.text = self.get_enclosing_src_file().src_lines[start_idx]
        self.scoped_name = self.get_scoped_name()
        self.start_idx = start_idx
        self.end_idx = -1
        self.use_statements = []
        self.symbol_table = {}


    def unparse_to_string(self):
        return [line if isinstance(line, str) else line.text for line in self.get_enclosing_src_file().src_lines[self.start_idx:self.end_idx + 1]]


    @Node.is_tainted.getter
    def is_tainted(self):

        # first check if the scope itself has deliberately been marked as tainted
        if self._is_tainted:
            return True
        
        # otherwise, check all nested scopes in the symbol table for taint
        # note that all interfaces are marked as tainted so, if their surrounding scope is unparsed, so are they
        for node in self.symbol_table.values():
            if isinstance(node, Scope) and node.is_tainted:
                return True
        return False


    def get_scoped_name(self):
                
        scoped_name = self.name
        node = self
        while not isinstance(node, SourceFile):
            node = node.parent
            if node.name:
                scoped_name = node.name + "::" + scoped_name

        return "::" + scoped_name


    def unparse(self, src_lines, file):
        
        # for all non-top-level source file scopes, write out the first line of the scope
        if not isinstance(self, SourceFile):
            file.write(src_lines[self.start_idx].text)
        
        i = self.start_idx
        while i + 1 < self.end_idx:
            i += 1

            if isinstance(src_lines[i], str):

                # keep strings with contains, implicit, private, and public modifiers as well as cpp directives
                if parsing.keep_line(src_lines[i]) or parsing.cpp_directive(src_lines[i]):
                    file.write(src_lines[i])

                # if there is a line starting with a fixed form line continuation that is preceded by a cpp directive,
                # we include it if the previous "node"-type is tainted (this happens in ADCIRC for example with multiline
                # variable declarations where some variables are conditional upon CPP directives)
                elif self.get_enclosing_src_file().fixed_form_fortran and parsing.is_fixed_form_continuation_line(src_lines[i]):
                    j = i
                    while not isinstance(src_lines[j], Node):
                        j -= 1
                    if src_lines[j].is_tainted:
                        file.write(src_lines[i])

            elif isinstance(src_lines[i], Scope):
                if src_lines[i].is_tainted:
                    src_lines[i].unparse(src_lines, file)
                i = src_lines[i].end_idx

            elif isinstance(src_lines[i], Statement):
                if src_lines[i].is_tainted:
                    src_lines[i].unparse(file)

        if not isinstance(self, SourceFile):
            file.write(src_lines[self.end_idx])


    def resolve_references(self, src_lines):
        i = self.start_idx
        while i + 1 < self.end_idx:
            i += 1

            if isinstance(src_lines[i], Scope):
                i = src_lines[i].resolve_references(src_lines)
            elif isinstance(src_lines[i], (UseStatement, StatementWithProcedureCalls, VariableDeclaration, ProceduresInInterface)):
                try:
                    src_lines[i].resolve_references()
                except (ExcludedModule, ProcedureWithoutFPArguments, IntrinsicOrOmittedProcedure):
                    to_delete = src_lines[i]
                    src_lines[i] = to_delete.text
                    del to_delete
                
        return self.end_idx


    def propagate_taint(self, src_lines, new_taint_count):
        i = self.start_idx
        while i + 1 < self.end_idx:
            i += 1

            if isinstance(src_lines[i], Scope):
                i, new_taint_count = src_lines[i].propagate_taint(src_lines, new_taint_count)
            elif isinstance(src_lines[i], (UseStatement, StatementWithProcedureCalls, VariableDeclaration)):
                new_taint_count = src_lines[i].propagate_taint(new_taint_count)

        return self.end_idx, new_taint_count


    def _scope_end(self, src_lines, i):
        return parsing.end_statement(src_lines[i])


    def parse(self, src_lines):

        i = self.start_idx
        while not self._scope_end(src_lines, i + 1):
            i += 1

            name = parsing.program_begin(src_lines[i])
            if name:
                program = Program(name, src_lines, i, parent=self)
                src_lines[i] = program
                i = program.end_idx
                continue

            name = parsing.module_begin(src_lines[i])
            if name:
                module = Module(name, src_lines, i, parent=self)
                src_lines[i] = module
                self.symbol_table[module.scoped_name] = module
                i = module.end_idx
                continue

            name = parsing.contains_use_statement(src_lines[i])
            if name:
                use_statement = UseStatement(module_name=name, text=src_lines[i], parent=self)
                src_lines[i] = use_statement
                self.use_statements.append(use_statement)
                continue

            match = parsing.contains_variable_declaration(src_lines[i])
            if match:
                variable_declaration = VariableDeclaration(text=src_lines[i], parent=self)
                src_lines[i] = variable_declaration
                for variable in variable_declaration.resolved_references.values():
                    self.symbol_table[variable.scoped_name] = variable
                continue

            name = parsing.procedure_begin(src_lines[i])
            if name:
                procedure = Procedure(name, src_lines, i, parent=self)
                src_lines[i] = procedure
                self.symbol_table[procedure.scoped_name] = procedure
                i = procedure.end_idx
                continue

            name = parsing.derived_type_begin(src_lines[i])
            if name:
                derived_type = DerivedType(name, src_lines, i, parent=self)
                src_lines[i] = derived_type
                self.symbol_table[derived_type.scoped_name] = derived_type
                i = derived_type.end_idx
                continue

            name = parsing.interface_begin(src_lines, i)
            if name:
                interface = Interface(name, src_lines, i, parent=self)
                src_lines[i] = interface
                self.symbol_table[interface.scoped_name] = interface
                i = interface.end_idx
                continue

            if parsing.contains_potential_procedure_calls(src_lines[i]):
                statement_with_procedure_calls = StatementWithProcedureCalls(text=src_lines[i], parent=self)
                src_lines[i] = statement_with_procedure_calls
                continue

        self.end_idx = i + 1


    def register_downstream_src(self, src_lines, module_name_to_downstream_src_path_map):
    
        global SETUP
        global EXCLUDED_NAMES

        # for any module imported into this source file, find src containing a declaration of that module
        for use_statement in self.use_statements:

            module_name = use_statement.module_name

            if module_name in module_name_to_downstream_src_path_map or module_name in EXCLUDED_NAMES:
                continue

            # if found, it should be unique
            # if not found, mark it as excluded
            unique_found = False
            multiple_found = []
            for src_search_path in SETUP['machine']['src_search_paths'].split("|"):
                result = subprocess.run(
                    f"find {src_search_path} -type f -iname '*.f90' -o -type f -iname '*.f' | xargs -r -n 128 -P 16 grep -lirE '^\s*module\s+{module_name}(\s+|$)'",
                    shell=True,
                    stdout=subprocess.PIPE,
                    text=True,
                    executable="/bin/bash",
                )
                grep_results = [os.path.relpath(x) for x in result.stdout.split("\n") if x]
                if len(grep_results) == 1:
                    unique_found = True
                    module_name_to_downstream_src_path_map[module_name] = grep_results[0]
                    break
                elif len(grep_results) > 1:
                    multiple_found += grep_results
            if not unique_found and multiple_found == []:
                EXCLUDED_NAMES.add(module_name)
                continue
            elif not unique_found and multiple_found != []:
                print(f"multiple matches for {module_name} found: {multiple_found}")
                assert(False)

        i = self.start_idx
        while i + 1 < self.end_idx:
            i += 1
            if isinstance(src_lines[i], Scope):
                src_lines[i].register_downstream_src(src_lines, module_name_to_downstream_src_path_map)


class Statement(Node):

    def __init__(self, text, parent):
        super().__init__(parent)
        self.text = text


    def resolve_references(self):
        pass
    

    def propagate_taint(self, new_taint_count):
        return new_taint_count


    def unparse(self, file):
        if self.is_tainted:
            file.write(self.text)


class Resolver(Statement):

    def __init__(self, text, parent):
        super().__init__(text, parent)
        self.resolved_references = {}
        self.unresolved_references = parsing.find_valid_fortran_names(text)
        

    def resolve_references(self):
        
        i = 0
        while i < len(self.unresolved_references):

            # attempt to resolve the reference
            # TODO: handle nested field accesses
            if "%" in self.unresolved_references[i]:
                tokens = [x for x in self.unresolved_references[i].split("%", 1) if x]
                if len(tokens) == 2:
                    node = self._resolve_derived_type_field_access(derived_type_name=tokens[0], derived_type_field_reference=tokens[1])
                else:
                    node = self._resolve_other_reference(self.unresolved_references[i], reference_types=(DerivedType,))
            else:
                node = self._resolve_other_reference(self.unresolved_references[i], reference_types=(Variable,))
            # save if it was resolved; otherwise, move on
            # we save unresolved references for debugging purposes
            if node:
                self.resolved_references[node.scoped_name] = node
                self.unresolved_references.pop(i)
            else:
                i += 1


    def _resolve_derived_type_field_access(self, derived_type_name, derived_type_field_reference):

        # iterate through parent scopes looking for a match
        # until either a match is found or there are no more scopes to search
        node = self
        while node.parent and not isinstance(node.parent, Scope):
            node = node.parent
        while node.parent:
            node = node.parent

            # check for locally-declared matches
            for child_node in node.symbol_table.values():
                if isinstance(child_node, Variable) and child_node.name == derived_type_name:
                    for reference in child_node.parent.resolved_references.values():
                        if isinstance(reference, DerivedType):
                            for grandchild_node in reference.symbol_table.values():
                                if grandchild_node.name == derived_type_field_reference:
                                    return grandchild_node

            # check for imported matches
            for use_statement in node.use_statements:
                for child_node in use_statement.resolved_references.values():
                    if isinstance(child_node, Variable) and child_node.name == derived_type_name:
                        for reference in child_node.parent.resolved_references.values():
                            if isinstance(reference, DerivedType):
                                for grandchild_node in reference.symbol_table.values():
                                    if grandchild_node.name == derived_type_field_reference:
                                        return grandchild_node

        # no matches found
        # check for a pointer variable which could be pointing to the derived type
        # note that this doesn't proceed to resolve the field reference
        node = self._resolve_other_reference(derived_type_name, reference_types=(Variable,))

        return node


    def _resolve_other_reference(self, target_name, reference_types):
        
        if len(reference_types) == 1 and reference_types[0] == DerivedType:
            target_name = target_name.replace("%", "")

        # iterate through parent scopes looking for a match
        # until either a match is found or there are no more scopes to search
        node = self
        while node.parent and not isinstance(node.parent, Scope):
            node = node.parent
        while node.parent:
            node = node.parent

            # check for locally-declared matches
            for child_node in node.symbol_table.values():
                if isinstance(child_node, reference_types) and child_node.name == target_name:
                    return child_node

            # check for imported matches unless the parent is an interface which should only contain procedure references 
            if not isinstance(self.parent, Interface):
                for use_statement in node.use_statements:
                    for child_node in use_statement.resolved_references.values():
                        if isinstance(child_node, reference_types) and child_node.name == target_name:
                            return child_node

        # no matches found
        return None


class Tainter(Resolver):
    
    def __init__(self, text, parent):
        super().__init__(text, parent)


    def propagate_taint(self, new_taint_count):
        
        # skip this statement if it has already been designated as tainted
        if not self.is_tainted:

            # check if any variable references are tainted or if any procedure references are targeted
            for child_node in self.resolved_references.values():
                if isinstance(child_node, Variable):
                    if self._taint_condition(child_node):
                        self.is_tainted = True
                        break
                elif isinstance(child_node, Procedure):
                    if child_node.is_targeted:
                        self.is_tainted = True
                        break

            # if so, the taint must be propagated
            if self.is_tainted:

                # propagate taint upward from the tainted statement to it's enclosing scope;
                # if the scope is a procedure, also propagate taint downward from the procedure to the dummy variables
                node = self.parent
                while not isinstance(node, Scope):
                    node = node.parent
                node.is_tainted = True
                if isinstance(node, Procedure):
                    for variable in node.dummy_variables:
                        if not variable.is_tainted:
                            variable.is_tainted = True
                            new_taint_count += 1

                # if this reference is a derived type symbol, then propagate taint to the
                # derived type and to all of the fields within that type
                elif isinstance(node, DerivedType):
                    for reference in node.symbol_table.values():
                        if not reference.is_tainted:
                            reference.is_tainted = True
                            new_taint_count += 1

                # propagate taint downward from the statement to all other references in that statement
                for child_node in self.resolved_references.values():    
                    
                    if child_node.is_tainted:
                        continue
                    else:
                        child_node.is_tainted = True

                    # if this reference is a variable, taint it
                    if isinstance(child_node, Variable):
                        new_taint_count += 1

                        # if this newly-tainted variable is declared within a procedure,
                        # propagate taint to the procedure and to all of its dummy variables
                        parent_node = child_node.parent
                        while not isinstance(parent_node, Scope):
                            parent_node = parent_node.parent
                        parent_node.is_tainted = True
                        if isinstance(parent_node, Procedure):
                            for variable in parent_node.dummy_variables:
                                if not variable.is_tainted:
                                    variable.is_tainted = True
                                    new_taint_count += 1
                    
                    # if this reference is a procedure symbol, then propagate taint to the
                    # procedure and to all of the dummy variables in that procedure
                    elif isinstance(child_node, Procedure):
                        for variable in child_node.dummy_variables:
                            if not variable.is_tainted:
                                variable.is_tainted = True
                                new_taint_count += 1

                    # if this reference is a derived type symbol, then propagate taint to the
                    # derived type and to all of the fields within that type
                    elif isinstance(child_node, DerivedType):
                        for reference in child_node.symbol_table.values():
                            if not reference.is_tainted:
                                reference.is_tainted = True
                                new_taint_count += 1

                    # if this reference is an interface symbol, then propagate taint to any
                    # procedures referenced within as well as their dummy variables
                    elif isinstance(child_node, Interface):
                        for reference in child_node.symbol_table.values():
                            reference.is_tainted = True
                            for procedure_reference in reference.resolved_references.values():
                                if not procedure_reference.is_tainted:
                                    procedure_reference.is_tainted = True
                                    for variable in procedure_reference.dummy_variables:
                                        if not variable.is_tainted:
                                            variable.is_tainted = True
                                            new_taint_count += 1

        return new_taint_count


    def _taint_condition(self):
        """To be implemented in the specific subclass"""
        raise NotImplementedError


class SourceFile(Scope):

    def __init__(self, src_path):
        
        self.src_path = os.path.relpath(src_path)
        with open(self.src_path, "r") as f:
            self.src_lines = f.readlines()

        super().__init__(name="", parent=None, start_idx=-1)

        self.upstream_src_paths = set()
        self.module_name_to_downstream_src_path_map = {}
        self.fixed_form_fortran = self.src_path.lower().endswith(".f")

        self.parse()
        self.end_idx = len(self.src_lines) - 1
        self.is_targeted = False


    @property
    def module_names(self):
        return [node.name for node in self.symbol_table.values() if isinstance(node, Module)]


    def _scope_end(self, src_lines, i):
        return i + 1 >= len(src_lines)


    def unparse(self):

        file = open(self.src_path, "w+")
        super().unparse(
            src_lines = self.src_lines,
            file = file
        )
        file.close()
        copy(self.src_path, os.path.join(os.environ['PROSE_EXPERIMENT_DIR'], "prose_workspace/original_files", os.path.basename(self.src_path) + ".slice"))


    def parse(self):

        global SETUP
        global EXCLUDED_NAMES

        # perform preprocessing of source code
        self.src_lines = parsing.preprocess(
            src_lines = self.src_lines,
            SETUP = SETUP,
            rose_preprocessing = True,
            excluded_names=EXCLUDED_NAMES,
            fixed_form_fortran = self.fixed_form_fortran
        )

        # initialize parsing
        super().parse(src_lines=self.src_lines)

        # register paths to upstream and downstream source code
        self.register_upstream_src()
        self.register_downstream_src()


    def register_downstream_src(self):
        super().register_downstream_src(
            src_lines = self.src_lines,
            module_name_to_downstream_src_path_map = self.module_name_to_downstream_src_path_map,
        )


    def resolve_references(self):
        super().resolve_references(self.src_lines)

    
    def propagate_taint(self):
        _, new_taint_count = super().propagate_taint(src_lines=self.src_lines, new_taint_count=0)
        return new_taint_count

        
    def register_upstream_src(self):

        global SETUP

        # for any module declared in this source file, find src containing "use" statements that import that module
        for module_name in self.module_names:    
            for src_search_path in SETUP['machine']['src_search_paths'].split("|"):
                result = subprocess.run(
                    f"find {src_search_path} -type f -iname '*.f90' -o -type f -iname '*.f' | xargs -r -n 128 -P 16 grep -lirE '^\s*use\s+{module_name}(,|\s+|$)'",
                    shell=True,
                    stdout=subprocess.PIPE,
                    text=True,
                    executable="/bin/bash",
                )

                # if we found upstream src files on this src_search_path, parse them and break
                paths = [os.path.relpath(x) for x in result.stdout.split("\n") if x]
                if paths:
                    self.upstream_src_paths = set(paths)
                    break


class Program(Scope):
    def __init__(self, name, src_lines, start_idx, parent):
        super().__init__(name, start_idx, parent)
        self.parse(src_lines)


class Module(Scope):
    def __init__(self, name, src_lines, start_idx, parent):
        super().__init__(name, start_idx, parent)
        self.parse(src_lines)


class DerivedType(Scope):
    def __init__(self, name, src_lines, start_idx, parent):
        super().__init__(name, start_idx, parent)
        self.parse(src_lines)


class Interface(Scope):
    def __init__(self, name, src_lines, start_idx, parent):
        super().__init__(name, start_idx, parent)
        self.parse(src_lines)


    def parse(self, src_lines):

        if "nameless" in self.name:
            super().parse(src_lines)
        else:
            i = self.start_idx
            while not self._scope_end(src_lines, i + 1):
                i += 1
                names = parsing.procedures_in_interface(src_lines[i])
                if names:
                    procedure_references = ProceduresInInterface(names=names, text=src_lines[i], parent=self)
                    self.symbol_table[" ".join(names)] = procedure_references
                    src_lines[i] = procedure_references

            self.end_idx = i + 1


class ProceduresInInterface(Resolver):
    def __init__(self, names, text, parent):
        self.names = names
        super().__init__(text=" ".join(names), parent=parent)
        self.text = text


    def resolve_references(self):
        i = 0
        while i < len(self.unresolved_references):

            # attempt to resolve the reference
            node = self._resolve_other_reference(self.unresolved_references[i], reference_types=(Procedure,))            

            # save if it was resolved; otherwise, move on
            # we save unresolved references for debugging purposes
            if node:
                self.resolved_references[node.scoped_name] = node
                self.unresolved_references.pop(i)
            else:
                i += 1


class Procedure(Scope):
    def __init__(self, name, src_lines, start_idx, parent):
        super().__init__(name, start_idx, parent)
        self.parse(src_lines)
        self.dummy_variables = self._get_dummy_variables(src_lines[start_idx])


    @property
    def is_targeted(self):
        for node in self.dummy_variables:
            if node.is_targeted:
                return True
        return False


    def _get_dummy_variables(self, line):

        arg_text = ""
        text = line[line.find(self.name) + len(self.name):]
        i = text.find("(")
        if i == -1:
            return []
        unbalanced_paren_count = 1
        while unbalanced_paren_count > 0:
            i += 1
            arg_text += text[i]
            if text[i] == "(":
                unbalanced_paren_count += 1
            elif text[i] == ")":
                unbalanced_paren_count -= 1

        dummy_variable_names = parsing.find_valid_fortran_names(arg_text)
        if parsing.is_function(line):
            return_var_name = parsing.function_return_value(text)
            if not return_var_name:
                return_var_name = self.name
            dummy_variable_names.append(return_var_name)

        dummy_variables = []
        for name in dummy_variable_names:
            try:
                dummy_variables.append(self.symbol_table[self.scoped_name + "::" + name])
            except KeyError:
                assert(name == return_var_name)

        return dummy_variables


class StatementWithProcedureCalls(Statement):
    
    def __init__(self, text, parent):
        super().__init__(text, parent)
        self.procedure_calls = []
        for match in parsing.contains_potential_procedure_calls(self.text):
            self.procedure_calls.append(ProcedureCall(text=match, parent=self))


    def resolve_references(self):
        deferred_exceptions = []
        for procedure_call in self.procedure_calls:
            try:
                procedure_call.resolve_references()
            except (ProcedureWithoutFPArguments, IntrinsicOrOmittedProcedure) as e:
                deferred_exceptions.append(e)

        # if all procedure calls raised an exception, raise the latest exception
        if len(deferred_exceptions) == len(self.procedure_calls):
            raise deferred_exceptions[-1]


    def unparse(self, file):
        if parsing.elseif_statement(self.text):
            self.text = self.text.replace("else", "", 1)
        if parsing.multiline_if_statement(self.text):
            self.text = self.text.rstrip() + "; endif\n"
        super().unparse(file)


    def propagate_taint(self, new_taint_count):
        for procedure_call in self.procedure_calls:
            new_taint_count += procedure_call.propagate_taint(new_taint_count)

        return new_taint_count


    @Node.is_tainted.getter
    def is_tainted(self):
        return any([procedure_call.is_tainted for procedure_call in self.procedure_calls])


class ProcedureCall(Tainter):
     
    def __init__(self, text, parent):
        arg_text = parsing.remove_fp_literals_from_argument_list(text[text.find("(") + 1:])
        self.contains_float_arguments = len(arg_text) != len(text[text.find("(") + 1:])
        super().__init__(arg_text, parent)
        self.name = text[:text.find("(")].split()[-1]


    def resolve_references(self):

        global PROCEDURE_CALL_DEPENDENCIES

        # resolve the procedure symbol
        node = self._resolve_other_reference(self.name, reference_types=(Procedure, Interface))
        if not node:
            if self.name != "sign":
                raise IntrinsicOrOmittedProcedure
        else:
            self.resolved_references[node.scoped_name] = node
            if self.get_enclosing_src_file().is_targeted:
                PROCEDURE_CALL_DEPENDENCIES.add(node.get_enclosing_src_file().src_path)

        # resolve the variable symbols
        super().resolve_references()

        # if there are no floating point arguments, signal this via exception
        # this node will be deleted
        if not self.contains_float_arguments:
            for child in self.resolved_references.values():
                if isinstance(child, Variable) and child.is_float_type:
                    self.contains_float_arguments = True
                    break
            if not self.contains_float_arguments:
                raise ProcedureWithoutFPArguments

    def _taint_condition(self, variable):
        return variable.is_targeted


class VariableDeclaration(Tainter):

    def __init__(self, text, parent):
        super().__init__(text, parent)
        self.is_float_type = parsing.is_real_variable_declaration(self.text)
        self.declared_variables = self._process_declared_variables()


    def _taint_condition(self, variable):
        return (variable.is_targeted or variable.is_tainted) and variable.scoped_name in self.declared_variables


    def _process_declared_variables(self):

        declared_variables = {}

        # handle case of the old-school pointer declaration (https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnan/index.html)
        match = parsing.old_pointer_declaration(self.text)
        if match:
            parseable_text = re.sub(r"((\()|([a-z0-9_]+\s*\)))", " ", self.text, flags=re.IGNORECASE)
        else:
            parseable_text = self.text

        # extract and construct new variable objects for each variable in this variable declaration
        # we omit parentheticals because, while they contain the names of variables used in the 
        # declaration, they do not contain the names of the actual variables declared in the declaration
        for name in parsing.find_valid_fortran_names(parseable_text, ignore_parentheticals=True):
            variable = Variable(parent=self, name=name, is_float_type=self.is_float_type)
            self.resolved_references[variable.scoped_name] = variable
            declared_variables[variable.scoped_name] = variable
            self.unresolved_references.remove(name)

        return declared_variables


class UseStatement(Resolver):
    
    def __init__(self, module_name, text, parent):
        super().__init__(text, parent)
        self.module_name = module_name


    @Node.is_tainted.getter
    def is_tainted(self):
        external_module = MODULE_NAME_TO_SRC_FILE_MAP[self.module_name].symbol_table["::" + self.module_name]
        return external_module.is_tainted
    

    def resolve_references(self):

        global EXCLUDED_NAMES

        if self.module_name in EXCLUDED_NAMES:
            raise ExcludedModule

        global MODULE_NAME_TO_SRC_FILE_MAP

        external_module = MODULE_NAME_TO_SRC_FILE_MAP[self.module_name].symbol_table["::" + self.module_name]
        for node in external_module.symbol_table.values():
            self.resolved_references[node.scoped_name] = node


class Variable(Node):
    
    def __init__(self, parent, name, is_float_type):
        super().__init__(parent)
        self.name = name
        self.scoped_name = self._get_scoped_name()
        self.is_float_type = is_float_type
        self.is_targeted = self.check_if_targeted()
        self.is_tainted = self.is_targeted


    def _get_scoped_name(self):
        node = self.parent
        while not isinstance(node, Scope):
            node = node.parent
        return node.get_scoped_name() + "::" + self.name


    def check_if_targeted(self):

        global SETUP
        global SEARCH_SPACE

        targeted = 0
        if self.is_float_type:

            for search_pattern in SETUP['target']['search_patterns'].split("|"):
                if search_pattern and re.search(r"{}".format(search_pattern), self.scoped_name, re.IGNORECASE):
                    targeted = 1
                    break

            for ignore_pattern in SETUP['target']['ignore_patterns'].split("|"):
                if ignore_pattern and re.match(r"{}".format(ignore_pattern), self.scoped_name, re.IGNORECASE):
                    targeted = -1
                    break

        if targeted == 1:
            SEARCH_SPACE.add(f"{self.scoped_name},$kind$\n")

        return targeted == 1


class ProcedureWithoutFPArguments(Exception):
    pass


class IntrinsicOrOmittedProcedure(Exception):
    pass


class ExcludedModule(Exception):
    pass
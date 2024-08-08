#include "rose.h"
#include "plugin.h"
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <cmath>
#include <cctype>
#include <boost/algorithm/string.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/filesystem.hpp>

#define DEFAULT_KIND 8
#define OPTIONAL_ARG_NOT_PROVIDED -1
#define NON_FLOAT_TYPE 0
#define NOT_THE_SAME_PROCEDURE 1
#define UNRESOLVED_FUNCTION -2
#define GENERATE_PROJECT_SOURCE_LIST 0

#define VERBOSE 0

using namespace std;
using namespace Rose;


/************************************
 * BEGIN requirements for boost graph
 ************************************/

// graph internal properties
struct edge_vars_t {
  typedef boost::edge_property_tag kind;
  static size_t const num;
};
size_t const edge_vars_t::num = (size_t)&edge_vars_t::num;

struct vertex_kind_t {
  typedef boost::vertex_property_tag kind;
  static size_t const num;
};
size_t const vertex_kind_t::num = (size_t)&vertex_kind_t::num;

// G_proc property typedefs
typedef boost::property<boost::vertex_name_t, string> GP_VertexProperty;
typedef boost::property<boost::edge_name_t, string> GP_EdgeProperty;

// G_var property typedefs
typedef boost::property<boost::vertex_name_t, string,
                boost::property<vertex_kind_t, int>> GV_VertexProperty;
typedef boost::property<boost::edge_weight_t, double> GV_EdgeProperty;

// typedef for G_proc and G_var graphs
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, GP_VertexProperty, GP_EdgeProperty> G_proc_t;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, GV_VertexProperty, GV_EdgeProperty> G_var_t;

// G_proc property maps
typedef boost::property_map<G_proc_t, boost::edge_name_t>::type GP_EdgeNameMap; //TODO is this necessary?
typedef boost::property_map<G_proc_t, boost::vertex_name_t>::type GP_VertexNameMap;

// G_var property maps
typedef boost::property_map<G_var_t, boost::vertex_name_t>::type GV_VertexNameMap;
typedef boost::property_map<G_var_t, vertex_kind_t>::type GV_VertexKindMap;
typedef boost::property_map<G_var_t, boost::edge_weight_t>::type GV_EdgeWeightMap;

/**********************************
 * END requirements for boost graph
 **********************************/

struct variable_binding {
  double weight;
  vector<pair<string,int>> binding;
};

string WORKING_DIR;
vector<double> EXECUTION_COUNTS;
vector<string> ORIGINAL_CODE_TEXT;

string get_scoped_name( SgNode* node );
SgType* get_expression_type( SgExpression* expression );
int get_array_rank( SgExpression* expression, SgArrayType* arrayExpressionType );
void set_real_kind( SgInitializedName* targetNode, SgType* replacementType );
SgScopeStatement* get_specified_scope( SgProject* n, string targetScopeString );
SgTypeFloat* is_underlying_real_type( SgNode* node );
void check_args( SgFunctionCallExp* functionCall, SgProcedureHeaderStatement* calleeProcHeader, bool* mismatchPtr = NULL, map<string,pair<int, int>>* argKindConfigsPtr = NULL, 
  vector<variable_binding>* accumulatedBoundVariablesPtr = NULL, vector<string>* accumulatedConstantListPtr = NULL );
vector<SgInitializedName*> get_dummy_var_decl_order( SgProcedureHeaderStatement* procHeader );
void preserve_aliasing( vector<SgNode*> functionCalls );
SgInitializedName* is_constant( SgNode* node );
unordered_set<SgSourceFile*> parse_source_files_into_AST( SgProject* n, vector<string> filePaths, string unparseLocation );
void intraprocedural_variable_reasoning( vector<SgNode*> binaryOps, vector<variable_binding> & intraboundVariables );
double get_edgeWeight_from_profiling_info( SgNode* targetNode );
void init_edgeWeight_info( SgSourceFile* sourceFile );
double calc_similarity( string s1, string s2 );
SgInitializedName* find_imported_variable_declaration( SgInitializedName* var, SgSourceFile* originalFile );

/*
 * GenerateGraph is called once for each source code file in the entire project
 *
 * It takes one command line argument:
 *    (1) the path to the working directory
 *
 * It does the following:
 *    (1) incrementally constructs a graph G_proc, in which each vertex represents
 *          a procedure and each edge represents a procedure call
 *    (2) incrementally constructs a graph G_var, in which each vertex
 *          represents an fp variable and each edge represents a fp dataflow
 *          between variables
 *    (3) writes out a collection of all interprocedural fp dataflow
 *          and all constants for constant propagation in the python framework
 *    (4) writes out information on each source file used in the
 *          python framework for the construction of ProseSourceTransformers:
 *          path to the original source and all the scopes in that file
 */
class GenerateGraph : public Rose::PluginAction {
  public:
    GenerateGraph() {}
    ~GenerateGraph() {}

    SgScopeStatement* topLevelScopePtr;
    map<SgFunctionCallExp*, vector<SgFunctionDefinition*>> genericProcedureCallToSpecificProcedure;

    bool ParseArgs(const std::vector<std::string> &args){
      assert(args.size() == 1);
      WORKING_DIR = args[0];

      if ( WORKING_DIR[WORKING_DIR.size() - 1] != '/' ){
        WORKING_DIR = WORKING_DIR + '/';
      }

      return true;
    }

    void process (SgProject* n) {

      #if VERBOSE >= 1
      cout << "Start Action: GenerateGraph" << endl;
      #endif

      // declare necessary graph objects
      G_proc_t G_proc(0);
      GP_EdgeNameMap GP_EdgeName = boost::get(boost::edge_name, G_proc);
      GP_VertexNameMap GP_VertexName = boost::get(boost::vertex_name, G_proc);
      map<string, int> GP_vertexMap;
      map<string, int> GP_callMap;

      G_var_t G_var(0);
      GV_EdgeWeightMap GV_EdgeWeight = boost::get(boost::edge_weight, G_var);
      GV_VertexNameMap GV_VertexName = boost::get(boost::vertex_name, G_var);
      GV_VertexKindMap GV_VertexKind = boost::get(vertex_kind_t(), G_var);
      map<string, int> GV_vertexMap;

      // get global scope of target file; it should be the single
      // filename in the file list for the project that doesn't end
      // with the ".rmod" extension.
      for ( const auto& x : n->get_files() ){
        SgSourceFile* sourceFile = isSgSourceFile(x);
        assert(sourceFile);

        if ( ! boost::algorithm::ends_with(sourceFile->getFileName(), ".rmod") ){

          #if VERBOSE >= 1
          cout << "Processing file: " << sourceFile->getFileName() << endl;
          #endif

          topLevelScopePtr = sourceFile->get_globalScope();
          init_edgeWeight_info(sourceFile);
          break;
        }
      }
      assert(topLevelScopePtr != NULL);

      // gather all scopes by querying for SgFunctionDefinition nodes
      // and SgClassDefinition nodes
      vector<SgNode*> scopes;
      vector<SgInitializedName*> fp_variables;
      auto procedureDefinitions = NodeQuery::querySubTree(topLevelScopePtr, V_SgFunctionDefinition);
      auto classDefinitions = NodeQuery::querySubTree(topLevelScopePtr, V_SgClassDefinition);
      scopes.insert(scopes.end(), procedureDefinitions.begin(), procedureDefinitions.end());
      scopes.insert(scopes.end(), classDefinitions.begin(), classDefinitions.end());

      #if VERBOSE >= 1
      cout << "Creating vertices in G_proc for all discovered scopes" << endl;
      cout << "AND Creating vertices in G_var for all floating-point variables declared in those scopes" << endl;
      #endif

      // for each scope, make sure there is a corresponding vertex in
      // G_proc; for each floating point variable in each scope, make
      // sure there is a corresponding vertex in G_var
      for ( const auto& x : scopes ){

        SgScopeStatement* scope = isSgScopeStatement(x);
        assert( scope != NULL );
        string scopeName = get_scoped_name(scope);

        #if VERBOSE >= 2
        cout << "\t Processing scope: " << scopeName << endl;
        #endif

        // get symbols from this scope's symbol table and its body's
        // symbol table as well if it is a procedure definition
        auto symbolsInScope = scope->get_symbol_table()->get_symbols();
        if ( SgFunctionDefinition* procDef = isSgFunctionDefinition(scope) ){
          auto temp = procDef->get_body()->get_symbol_table()->get_symbols();
          symbolsInScope.insert(temp.begin(), temp.end());
        }

        // search GP_vertexMap to see if a vertex for this scope has
        // already been created. If not, create one.
        G_proc_t::vertex_descriptor fromScope_vertexD;
        if ( GP_vertexMap.count(scopeName) > 0 ){
          fromScope_vertexD = boost::vertex(GP_vertexMap[scopeName], G_proc);

          #if VERBOSE >= 3
          cout << "\t\t vertexD " << fromScope_vertexD << " already exists in G_proc" << endl;
          #endif
        }
        else{

          #if VERBOSE >= 3
          cout << "\t\t Adding new vertexD " << fromScope_vertexD << " to G_proc" << endl;
          #endif

          fromScope_vertexD = boost::add_vertex(G_proc);
          GP_vertexMap[scopeName] = fromScope_vertexD;
          boost::put(GP_VertexName, fromScope_vertexD, scopeName);
        }

        // process each floating point variable in this scope
        for ( auto& symbol : symbolsInScope ){
          if ( SgVariableSymbol* varSymbol = isSgVariableSymbol(symbol) ){
            if ( SgInitializedName* var = isSgInitializedName(varSymbol->get_declaration()) ){
              if ( SgTypeFloat* varFloatType = is_underlying_real_type(var) ){

                // make sure that the variable was declared in this file
                if ( SageInterface::getEnclosingSourceFile(topLevelScopePtr) == SageInterface::getEnclosingSourceFile(var->get_declaration()) ){

                  // save it for writing out information at the end of this plugin action
                  fp_variables.push_back(var);

                  string scopedVarName = get_scoped_name(var);

                  // manual workaround for bug revealed by ADCIRC in which variables
                  // declared with EXTERNAL are parsed into the AST as REAL but, when
                  // getting their scoped name, they don't have a prefix scope
                  if ( !boost::algorithm::starts_with(scopedVarName, "::") ){
                    continue;
                  }

                  #if VERBOSE >= 4
                  cout << "\t\t\t Processing variable: " << scopedVarName << endl;
                  #endif

                  // search GV_vertexMap to see if a vertex for this var has
                  // already been created. If not, create one.
                  G_var_t::vertex_descriptor var_vertexD;
                  if ( GV_vertexMap.count(scopedVarName) > 0 ){
                    var_vertexD = boost::vertex(GV_vertexMap[scopedVarName], G_var);

                    #if VERBOSE >= 5
                    cout << "\t\t\t\t vertexD " << var_vertexD << " already exists in G_var" << endl;
                    #endif
                  }
                  else{

                    #if VERBOSE >= 5
                    cout << "\t\t\t\t Adding new vertexD " << var_vertexD << " to G_var" << endl;
                    #endif

                    var_vertexD = boost::add_vertex(G_var);
                    GV_vertexMap[scopedVarName] = var_vertexD;

                    int kind = DEFAULT_KIND;
                    if ( SgIntVal* intVal = isSgIntVal(varFloatType->get_type_kind()) ){
                      kind = intVal->get_value();
                    }

                    boost::put(GV_VertexName, var_vertexD, scopedVarName);
                    boost::put(GV_VertexKind, var_vertexD, kind);
                  }
                }
              }
            }
          }
        }
      } // end construction of vertices in G_var and G_proc

      // structures for keeping track of information required for
      // constant propagation (to be written out and eventually read
      // into the python framework which performs the constant
      // propagation)
      vector<variable_binding> interboundVariables;
      vector<string> constantList;

      // query for all procedure calls
      auto functionCalls = NodeQuery::querySubTree(topLevelScopePtr, V_SgFunctionCallExp);

      #if VERBOSE >= 1
      cout << "Creating edges in G_proc for all discovered procedure calls" << endl;
      cout << "AND Creating edges in G_var for all discovered interprocedural floating-point dataflow" << endl;
      #endif

      // process all procedure calls
      for ( const auto& x : functionCalls ){
        if ( SgFunctionCallExp* functionCall = isSgFunctionCallExp(x) ){

          #if VERBOSE >= 2
          cout << "\t Processing procedure call: " << functionCall->unparseToString() << endl;
          #endif

          // make sure that this procedure call has at least one floating point argument before continuing
          bool floatFlag = false;
          for ( const auto& exp : functionCall->get_args()->get_expressions() ){
            if ( isSgTypeFloat(get_expression_type(exp)->findBaseType()) ){
              floatFlag = true;
              break;
            }
          }
          if ( !floatFlag ){

            #if VERBOSE >= 3
            cout << "\t\t SKIPPING (no floats)" << endl;
            #endif

            continue;
          }

          // get vertex in G_proc corresponding to the fromScope
          // this vertex will have been added above
          SgScopeStatement* fromScope = SageInterface::getEnclosingScope(functionCall);
          while ( !(isSgFunctionDefinition(fromScope) || isSgClassDefinition(fromScope)) ){
            fromScope = SageInterface::getEnclosingScope(fromScope);
          }
          string fromScopeName = get_scoped_name(fromScope);
          assert( GP_vertexMap.count(fromScopeName) > 0 );
          G_proc_t::vertex_descriptor fromScope_vertexD = boost::vertex(GP_vertexMap[fromScopeName], G_proc);

          // get vertex in G_proc correponding to the toScope
          // start by trying to resolve the call
          vector<SgFunctionDefinition*> toScopes = resolve_function_call_to_definition(functionCall, genericProcedureCallToSpecificProcedure);

          // if there are no candidate resolutions found, the function
          // is considere obscured, i.e., its definition is not
          // visible to Prose. In this case, SKIP
          if ( toScopes.size() == 0 ){

            #if VERBOSE >= 3
            cout << "\t\t SKIPPING (call to obscured function " <<  isSgFunctionRefExp(functionCall->get_function())->get_symbol()->get_name() << " )" << endl;
            #endif

            continue;
          }
          // if there are multiple candidate resolutions, we have a
          // sound overapproximation. This is likely due to some
          // shortcomings in our ability to uniquely resolve the call
          // statically, probably due to a lack of good info in the
          // AST represantation. Log such cases!
          else if ( toScopes.size() > 1 ){
            ofstream errFile;
            errFile.open(WORKING_DIR + "prose_logs/failed_name_resolutions.txt", ios_base::app);
            if ( errFile.is_open() ){
              errFile << "\tfile: " <<  SageInterface::getEnclosingSourceFile(topLevelScopePtr)->get_sourceFileNameWithPath() << endl;
              errFile << "\tline:" << functionCall->unparseToString() << endl;
              errFile.close();
            }
            else{
              assert(false);
            }
          }

          // (at this program point, there will only be a single toScope if resolve_function_call_to_definition was able to
          // resolve the functionCall to a unique SgProcedureHeaderStatement; otherwise, it's a generic interface where we
          // couldn't resolve for whatever reason so we overapproximate and include all of the specific procedures in the interface)
          for ( const auto& toScope : toScopes ){
            string toScopeName = get_scoped_name(toScope);

            // Finally, get the vertex descriptor of the toScope.
            // Search GP_vertexMap to see if a vertex for this toScope
            // has already been created. If not, create one.
            G_proc_t::vertex_descriptor toScope_vertexD;
            if ( GP_vertexMap.count(toScopeName) ){
              toScope_vertexD = boost::vertex(GP_vertexMap[toScopeName], G_proc);
            }
            else{
              toScope_vertexD = boost::add_vertex(G_proc);
              GP_vertexMap[toScopeName] = toScope_vertexD;
              boost::put(GP_VertexName, toScope_vertexD, toScopeName);
            }

            #if VERBOSE >= 3
            cout << "\t\t updating G_proc with edge from " << fromScopeName << " to vertex " << toScopeName << endl;
            cout << "\t\t\t (" << fromScope_vertexD << ") -> (" << toScope_vertexD << ")" << endl;                
            #endif

            // to G_proc, add a new edge between the vertices
            // representing the toScope and the fromScope for this
            // particular procedure call. Must call preserve_aliasing
            // on this single function call in order to save the
            // correct scopedCallName into the G_proc edge
            G_proc_t::edge_descriptor GP_edgeD;
            bool success;
            boost::tie(GP_edgeD, success) = boost::add_edge(fromScope_vertexD, toScope_vertexD, G_proc);
            assert( success );
            preserve_aliasing({isSgNode(functionCall)});
            string scopedCallName = get_scoped_name(functionCall);
            boost::put(GP_EdgeName, GP_edgeD, scopedCallName);

            // if the resolution of the procedure call to the actual
            // SgProcedureHeaderStatement was unique, record this
            // resolution by mapping the unique procedure call name to
            // the vertex descriptor of the toScope
            if ( toScopes.size() == 1 ){
              GP_callMap[scopedCallName] = toScope_vertexD;
            }

            #if VERBOSE >= 3
            cout << "\t\t checking for interprocedural fp dataflow and constants" << endl;
            #endif

            // invoke check_args function to make note of any
            // interprocedural fp dataflow and any constants in this
            // procedure call
            SgProcedureHeaderStatement* calleeProcHeader = isSgProcedureHeaderStatement(toScope->get_declaration());
            assert( calleeProcHeader != NULL );
            map<string,pair<int, int>> argKindConfigs;
            vector<variable_binding> temp_interboundVariables;
            bool mismatch;
            check_args(functionCall, calleeProcHeader, &mismatch, &argKindConfigs, &temp_interboundVariables, &constantList);
            interboundVariables.insert(interboundVariables.end(), temp_interboundVariables.begin(), temp_interboundVariables.end());

            // update G_var vertices and edges based on discovered
            // interprocedural floating-point dataflow. This means
            // inspecting the bindings discovered by check_args,
            // making sure that all variables have corresponding
            // vertices in G_var, and then making sure that there is
            // an edge between the vertices representing the variables
            // in the binding
            for ( const auto& temp : temp_interboundVariables ){
              double edgeWeight = temp.weight;
              auto binding = temp.binding;

              // note that the toVar is always first in the bindings.
              // In the case of interprocedural dataflow, the toVars
              // in each binding are the dummy arguments of the called
              // procedure 
              string toVarScopedName;
              int toVarKind;
              boost::tie(toVarScopedName, toVarKind)= binding.front();

              // manual workaround for bug revealed by ADCIRC in which variables
              // declared with EXTERNAL are parsed into the AST as REAL but, when
              // getting their scoped name, they don't have a prefix scope
              if ( !boost::algorithm::starts_with(toVarScopedName, "::") ){
                continue;
              }

              // get the vertex descriptor of the vertex in G_var corresponding
              // to this toVar; create it if it does not yet exist
              G_var_t::vertex_descriptor toVar_vertexD;
              if ( GV_vertexMap.count(toVarScopedName) > 0 ){
                toVar_vertexD = boost::vertex(GV_vertexMap[toVarScopedName], G_var);
              }
              else{
                toVar_vertexD = boost::add_vertex(G_var);
                GV_vertexMap[toVarScopedName] = toVar_vertexD;

                boost::put(GV_VertexName, toVar_vertexD, toVarScopedName);
                boost::put(GV_VertexKind, toVar_vertexD, toVarKind);
              }

              // now process each of the fromVars in this binding
              // recall that we skip the first variable in the binding
              // since this is the toVar
              auto it = binding.begin();
              while ( ++it != binding.end() ){

                string fromVarScopedName;
                int fromVarKind;
                boost::tie(fromVarScopedName, fromVarKind) = *it;

                // manual workaround for bug revealed by ADCIRC in which variables
                // declared with EXTERNAL are parsed into the AST as REAL but, when
                // getting their scoped name, they don't have a prefix scope
                if ( !boost::algorithm::starts_with(fromVarScopedName, "::") ){
                  continue;
                }

                // get the vertex descriptor of the vertex in G_var corresponding
                // to this fromVar; create it if it does not yet exist
                G_var_t::vertex_descriptor fromVar_vertexD;
                if ( GV_vertexMap.count(fromVarScopedName) > 0 ){
                  fromVar_vertexD = boost::vertex(GV_vertexMap[fromVarScopedName], G_var);
                }
                else{
                  fromVar_vertexD = boost::add_vertex(G_var);
                  GV_vertexMap[fromVarScopedName] = fromVar_vertexD;

                  boost::put(GV_VertexName, fromVar_vertexD, fromVarScopedName);
                  boost::put(GV_VertexKind, fromVar_vertexD, fromVarKind);
                }

                #if VERBOSE >= 4
                cout << "\t\t\t updating G_var with edge between " << toVarScopedName << " and " << fromVarScopedName << endl;
                cout << "\t\t\t\t (" << toVar_vertexD << ") -> (" << fromVar_vertexD << ")" << endl;                
                #endif

                // to G_var, add a new edge between the vertices
                // representing the toVar and the fromVar in this
                // particular instance of interprocedural fp dataflow
                G_var_t::edge_descriptor GV_edgeD;
                bool success;
                boost::tie(GV_edgeD, success) = boost::add_edge(fromVar_vertexD, toVar_vertexD, G_var);
                if ( success ) {
                  boost::put(GV_EdgeWeight, GV_edgeD, edgeWeight);
                }
                else {
                  boost::put(GV_EdgeWeight, GV_edgeD, boost::get(GV_EdgeWeight, GV_edgeD) + edgeWeight);
                }
              }
            }
          }
        }
      }

      // preserve aliasing in any procedure calls that were not
      // processed above
      preserve_aliasing(functionCalls);

      #if VERBOSE >= 1
      cout << "Creating edges in G_var for all discovered intraprocedural floating-point dataflow" << endl;
      #endif

      vector<variable_binding> intraboundVariables;
      intraprocedural_variable_reasoning(NodeQuery::querySubTree(topLevelScopePtr, V_SgBinaryOp), intraboundVariables);

      // for all possible pairs in each intraprocudural binding,
      // increment the weight of the edges between the vertices
      // representing the variables in that pair (instantiating new
      // vertices and edges as necessary)
      for ( int index = 0; index < intraboundVariables.size(); index++) {
        auto binding = intraboundVariables[index].binding;
        double edgeWeight = intraboundVariables[index].weight;
        auto it1 = binding.begin() - 1;
        while ( ++it1 != binding.end() ){

          string scopedVarName1;
          int varKind1;
          boost::tie(scopedVarName1, varKind1) = *it1;

          // manual workaround for bug revealed by ADCIRC in which variables
          // declared with EXTERNAL are parsed into the AST as REAL but, when
          // getting their scoped name, they don't have a prefix scope
          if ( !boost::algorithm::starts_with(scopedVarName1, "::") ){
            continue;
          }

          auto it2 = it1;
          while ( ++it2 != binding.end() ){

            string scopedVarName2;
            int varKind2;
            boost::tie(scopedVarName2, varKind2) = *it2;

            // manual workaround for bug revealed by ADCIRC in which variables
            // declared with EXTERNAL are parsed into the AST as REAL but, when
            // getting their scoped name, they don't have a prefix scope
            if ( !boost::algorithm::starts_with(scopedVarName2, "::") ){
              continue;
            }

            // get vertex 1 or create it
            G_var_t::vertex_descriptor var1_vertexD;
            if ( GV_vertexMap.count(scopedVarName1) > 0 ){
              var1_vertexD = boost::vertex(GV_vertexMap[scopedVarName1], G_var);
            }
            else{
              var1_vertexD = boost::add_vertex(G_var);
              GV_vertexMap[scopedVarName1] = var1_vertexD;

              boost::put(GV_VertexName, var1_vertexD, scopedVarName1);
              boost::put(GV_VertexKind, var1_vertexD, varKind1);
            }

            // get vertex 2 or create it
            G_var_t::vertex_descriptor var2_vertexD;
            if ( GV_vertexMap.count(scopedVarName2) > 0 ){
              var2_vertexD = boost::vertex(GV_vertexMap[scopedVarName2], G_var);
            }
            else{
              var2_vertexD = boost::add_vertex(G_var);
              GV_vertexMap[scopedVarName2] = var2_vertexD;

              boost::put(GV_VertexName, var2_vertexD, scopedVarName2);
              boost::put(GV_VertexKind, var2_vertexD, varKind2);
            }

            #if VERBOSE >= 2
            cout << "\t Updating G_var with edge between " << scopedVarName1 << " and " << scopedVarName2 << endl;
            cout << "\t\t\t (" << var1_vertexD << ") -> (" << var2_vertexD << ")" << endl;                
            #endif

            // to G_var, add a new edge between the vertices
            // representing the toVar and the fromVar in this
            // particular instance of intraprocedural fp dataflow
            G_var_t::edge_descriptor GV_edgeD;
            bool success;
            boost::tie(GV_edgeD, success) = boost::edge(var1_vertexD, var2_vertexD, G_var);
            if( !success ){
              boost::tie(GV_edgeD, success) = boost::add_edge(var1_vertexD, var2_vertexD, G_var);
              boost::put(GV_EdgeWeight, GV_edgeD, edgeWeight);
            }
            else{
              boost::put(GV_EdgeWeight, GV_edgeD, boost::get(GV_EdgeWeight, GV_edgeD) + edgeWeight);
            }
          }
        }
      }

      // write out binary of G_proc and its vertexMap
      ofstream outFile;
      string fileNameWithoutPathExt = SageInterface::getEnclosingSourceFile(topLevelScopePtr)->get_sourceFileNameWithoutPath();
      fileNameWithoutPathExt = fileNameWithoutPathExt.substr(0, fileNameWithoutPathExt.find_last_of("."));
      outFile.open(WORKING_DIR + "prose_workspace/temp_graph/" + fileNameWithoutPathExt + "_G_proc.graph", ios_base::binary);
      if ( outFile.is_open() ){
        boost::archive::binary_oarchive writeArchive(outFile);
        writeArchive << G_proc;
        outFile.close();
      }
      else{
        assert(false);
      }
      outFile.open(WORKING_DIR + "prose_workspace/temp_graph/" + fileNameWithoutPathExt + "_G_proc_vertexMap.map", ios_base::binary);
      if ( outFile.is_open() ){
        boost::archive::binary_oarchive writeArchive(outFile);
        writeArchive << GP_vertexMap;
        outFile.close();
      }
      else{
        assert(false);
      }
      outFile.open(WORKING_DIR + "prose_workspace/temp_graph/" + fileNameWithoutPathExt + "_G_proc_callMap.map", ios_base::binary);
      if ( outFile.is_open() ){
        boost::archive::binary_oarchive writeArchive(outFile);
        writeArchive << GP_callMap;
        outFile.close();
      }
      else{
        assert(false);
      }

      // write out binary of G_var and its vertexMap
      outFile.open(WORKING_DIR + "prose_workspace/temp_graph/" + fileNameWithoutPathExt + "_G_var.graph", ios_base::binary);
      if ( outFile.is_open() ){
        boost::archive::binary_oarchive writeArchive(outFile);
        writeArchive << G_var;
        outFile.close();
      }
      else{
        assert(false);
      }
      outFile.open(WORKING_DIR + "prose_workspace/temp_graph/" + fileNameWithoutPathExt + "_G_var_vertexMap.map", ios_base::binary);
      if ( outFile.is_open() ){
        boost::archive::binary_oarchive writeArchive(outFile);
        writeArchive << GV_vertexMap;
        outFile.close();
      }
      else{
        assert(false);
      }

      // write out any interprocedural bindings (to be read in by
      // python framework to propagate constants)
      outFile.open(WORKING_DIR + "prose_workspace/temp_graph/" + fileNameWithoutPathExt + "_inter_bound_variables.txt", ios_base::app);
      if ( outFile.is_open() ){
        for ( const variable_binding& x : interboundVariables ){
          for ( const pair<string,int>& y : x.binding ){
            outFile << y.first << ";";
          }
          outFile << endl;
        }
        outFile.close();
      }
      else{
        assert(false);
      }

      // write out any discovered constants (to be read in by python
      // framework to be propagated)
      outFile.open(WORKING_DIR + "prose_workspace/temp_graph/" + fileNameWithoutPathExt + "_constant_list.txt", ios_base::app);
      if (outFile.is_open()) {
        for ( const string& x : constantList) {
          outFile << x << endl;
        }
        outFile.close();
      }
      else {
        assert(false);
      }

      // write out info on this source file for future use
      vector<string> containedScopes;
      for ( const auto& x : scopes ){

        // get info on this scope
        string fromScopeName;
        if ( SgFunctionDefinition* procDef = isSgFunctionDefinition(x) ){
          fromScopeName = get_scoped_name(procDef);
        }
        else if ( SgClassDefinition* classDef = isSgClassDefinition(x) ){
          fromScopeName = get_scoped_name(classDef);
        }

        containedScopes.push_back(fromScopeName);
      }
      outFile.open(WORKING_DIR + "prose_workspace/__source_data/" + fileNameWithoutPathExt, ios_base::out);
      if ( outFile.is_open() ){
        outFile << SageInterface::getEnclosingSourceFile(topLevelScopePtr)->get_sourceFileNameWithPath() << endl;

        for ( auto itr = containedScopes.begin(); itr != containedScopes.end(); ++itr ){
          outFile << *itr << "+";
        }

        outFile << endl;

        for (const auto& x : fp_variables) {
          if ( SgInitializedName* var = isSgInitializedName(x) ){
            if ( (isSgVariableDeclaration(var->get_declaration())) || (isSgProcedureHeaderStatement(var->get_declaration())) ){

              string typeString;
              SgType* varType;

              if ( SgArrayType* arrayType = isSgArrayType(var->get_type()) ){
                typeString = "array";
                varType = arrayType->findBaseType();

                double n_elements = 1.0;//get_edgeWeight_from_profiling_info(isSgVariableDeclaration(var->get_declaration()));

                if ( n_elements > 0 ){
                  typeString = typeString + "(" + to_string(n_elements) + ")";
                }
              }
              else if ( SgPointerType* pointerType = isSgPointerType(var->get_type()) ){
                typeString = "pointer";
                varType = pointerType->findBaseType();
              }
              else{
                typeString = "scalar";
                varType = var->get_type();
              }

              if ( SgTypeFloat* floatVarType = isSgTypeFloat(varType) ){
                if ( SgIntVal* kindType = isSgIntVal(floatVarType->get_type_kind()) ){
                  outFile << get_scoped_name(var) << ",variableType=" << typeString << ",kind=" << kindType->get_value() << endl;
                }else{
                  outFile << get_scoped_name(var) << ",variableType=" << typeString << "," << "kind=" << DEFAULT_KIND << endl;
                }
              }
            }
          }
        }

        outFile.close();
      }
      else{
        assert(false);
      }

      // exit without unnecessarily unparsing
      exit(0);
    } // end process()


private:


  vector<SgFunctionDefinition*> resolve_function_call_to_definition( SgFunctionCallExp* functionCall, map<SgFunctionCallExp*, vector<SgFunctionDefinition*>> &genericProcedureCallToSpecificProcedure ){

    vector<SgFunctionDefinition*> resolvedNames;

    if ( genericProcedureCallToSpecificProcedure.count(functionCall) > 0 ){
      resolvedNames = genericProcedureCallToSpecificProcedure[functionCall];
    }
    else{
      // If the function declaration node does not have a corresponding function definition
      // node readily returned by the get_definition() method, it is either:
      //    1. an "obscured" function (i.e., intrinsic or defined in a module outside of the target piece of software) -or-
      //    2. a procedure call with a generic name to an interface that is defined within the target file
      //
      // In the first case, an empty list will be returned. In the second case, the return will either be
      // a list with a single SgProcedureHeaderStatement pointer (if the name resolution was successful) or a list with
      // all candidate SgProcedureHeaderStatement* pointers contained in the corresponding generic interface
      if ( functionCall->getAssociatedFunctionDeclaration()->get_definition() == NULL ){
        resolvedNames = resolve_function_call_to_definition_helper(functionCall);
        if ( !resolvedNames.empty() ){
          genericProcedureCallToSpecificProcedure[functionCall] = resolvedNames;
        }
      }

      // On the other hand, if the function declaration node has a corresponding function definition node
      // readily returned by the get_definition() method, it is either:
      //    1. a procedure call with a generic name to an interface that is defined outside of the target file -or-
      //    2. a normal procedure call with a specific name
      else{
        if ( isSgRenameSymbol(isSgFunctionRefExp(functionCall->get_function())->get_symbol()) ){
          resolvedNames = resolve_function_call_to_definition_helper(functionCall);
        }
        else{
          resolvedNames = {functionCall->getAssociatedFunctionDeclaration()->get_definition()};
          assert(resolvedNames.front() != NULL);
        }
      }

      if ( !resolvedNames.empty() ){
        genericProcedureCallToSpecificProcedure[functionCall] = resolvedNames;
      }
    }

    return resolvedNames;
  }


  // OVERAPPROXIMATES
  vector<SgFunctionDefinition*> resolve_function_call_to_definition_helper( SgFunctionCallExp* functionCall ){

    #if VERBOSE >= 11
    cout << "\n[BEGIN] attempting to resolve generic name: " << isSgFunctionRefExp(functionCall->get_function())->get_symbol()->get_name() << endl;
    cout << "\t source code: " << functionCall->unparseToString() << endl;
    #endif

    SgFunctionSymbol* functionSymbol = isSgFunctionRefExp(functionCall->get_function())->get_symbol();
    string calledName = functionSymbol->get_name();
    vector<SgFunctionDefinition*> resolvedNames;
    SgInterfaceStatement* matchingInterface = isSgInterfaceStatement(functionCall); // initialize to NULL

    SgModuleStatement* enclosingModule = isSgModuleStatement(SageInterface::getEnclosingClassDeclaration(functionCall->getAssociatedFunctionDeclaration()));
    if ( enclosingModule != NULL ){
      vector<SgInterfaceStatement*> previouslyExistingInterfaces = enclosingModule->get_interfaces();

      // search through previously existing interfaces to see if the name matches the name of the funcDecl
      auto interfaces_it = previouslyExistingInterfaces.begin();
      while ( (matchingInterface == NULL) && (interfaces_it != previouslyExistingInterfaces.end()) ){   

        // if we found a matching generic name...
        if ( (*interfaces_it)->get_name() == calledName ){
          matchingInterface = *interfaces_it;
        }
        ++interfaces_it;
      }
    }

    // if we found a matching generic interface, iterate through the specific names declared within that interface
    if ( matchingInterface ){
      vector<SgInterfaceBody*> interfaceBodyList = matchingInterface->get_interface_body_list();
      for ( const auto& interfaceBody : interfaceBodyList ){
        if ( SgProcedureHeaderStatement* specificNameProcHeader = isSgProcedureHeaderStatement(interfaceBody->get_functionDeclaration()) ){

          bool mismatch;
          check_args(functionCall, specificNameProcHeader, &mismatch);

          if ( !mismatch ){
            resolvedNames.clear();
            resolvedNames.push_back(specificNameProcHeader->get_definition());

            #if VERBOSE >= 11
            cout << "[END] SUCCESSFUL attempt to resolve generic name; Resolved generic name to " << resolvedNames.front()->unparseToString() << endl;
            cout << endl;
            #endif

            return resolvedNames;
          }
          else{
            resolvedNames.push_back(specificNameProcHeader->get_definition());
          }
        }
      }
    }

    #if VERBOSE >= 11
    cout << "[END] FAILED attempt to resolve generic name; Identified " << resolvedNames.size() << " candidate resolutions" << endl;
    cout << endl;
    #endif

    return resolvedNames;
  }
}; // end GenerateGraph


/*
 * LinkGraph is called once on the source code file designated as the top-level entrypoint in setup.ini
 *
 * It takes one command line argument:
 *    (1) the path to prose_workspace which contains the temporary
 *          graph files to be linked (generated by the GenerateGraph action)
 *
 * It does the following:
 *    (1) combines all of the files generated by the GenerateGraph action
 *          which was called on all source code files in the project,
 *          possibly in parallel. These include fragments of G_var, G_proc,
 *          and text files of the constants and interprocedural bindings
 *          used for constant prop in the Python framework
 */
class LinkGraph : public Rose::PluginAction {
  public:
    LinkGraph() {}
    ~LinkGraph() {}

    string WORKING_DIR;

    bool ParseArgs(const std::vector<std::string> &args){
      assert(args.size() == 1);
      WORKING_DIR = args[0];

      if ( WORKING_DIR[WORKING_DIR.size() - 1] != '/' ){
        WORKING_DIR = WORKING_DIR + '/';
      }

      return true;
    }

    void process (SgProject* n) {

      #if VERBOSE >= 1
      cout << "Start Action: LinkGraph" << endl;
      #endif

      // declare structure for keeping track of bound variables
      vector<vector<string>> intraboundVariables;

      // declare structure for named constants
      vector<string> constantList;

      // declare necessary graph objects
      G_proc_t merged_G_proc(0);
      GP_EdgeNameMap merged_GP_EdgeName = boost::get(boost::edge_name, merged_G_proc);
      GP_VertexNameMap merged_GP_VertexName = boost::get(boost::vertex_name, merged_G_proc);
      map<string, int> merged_GP_vertexMap;
      map<string, int> merged_GP_callMap;

      G_var_t merged_G_var(0);
      GV_EdgeWeightMap merged_GV_EdgeWeight = boost::get(boost::edge_weight, merged_G_var);
      GV_VertexNameMap merged_GV_VertexName = boost::get(boost::vertex_name, merged_G_var);
      GV_VertexKindMap merged_GV_VertexKind = boost::get(vertex_kind_t(), merged_G_var);
      map<string, int> merged_GV_vertexMap;

      // look for all graph files in the temp_graph directory; extract
      // their names. Each name corresponds to multiple partial pieces
      // of data that must be merged
      vector<string> fileNames;
      boost::filesystem::path temp_graph_path(WORKING_DIR + "temp_graph/");
      for (auto& file : boost::make_iterator_range(boost::filesystem::directory_iterator(temp_graph_path), {})){
        std::string fileName = file.path().filename().string();
        size_t fileNameBeforeExtPos = fileName.find("_G_proc.graph");
        if (boost::filesystem::is_regular_file(file) && (fileNameBeforeExtPos != std::string::npos)) {
          fileNames.push_back(fileName.substr(0, fileNameBeforeExtPos));
        }
      }

      ifstream inFile;
      ofstream outFile;

      // merge partial G_proc and G_var graphs and their respective
      // vertexMaps into single entities for the whole project
      for (auto& fileName : fileNames) {

        // declare and load partial vertexMaps
        map<string,int> partial_GP_vertexMap;
        map<string,int> partial_GP_callMap;
        map<string,int> partial_GV_vertexMap;

        inFile.open(WORKING_DIR + "temp_graph/" + fileName + "_G_proc_vertexMap.map", ios::binary);
        if ( inFile.is_open() ){
          boost::archive::binary_iarchive readArchive(inFile);
          readArchive >> partial_GP_vertexMap;
          inFile.close();
        }
        else {
          assert(false);
        }
        inFile.open(WORKING_DIR + "temp_graph/" + fileName + "_G_proc_callMap.map", ios::binary);
        if ( inFile.is_open() ){
          boost::archive::binary_iarchive readArchive(inFile);
          readArchive >> partial_GP_callMap;
          inFile.close();
        }
        else {
          assert(false);
        }
        inFile.open(WORKING_DIR + "temp_graph/" + fileName + "_G_var_vertexMap.map", ios::binary);
        if ( inFile.is_open() ){
          boost::archive::binary_iarchive readArchive(inFile);
          readArchive >> partial_GV_vertexMap;
          inFile.close();
        }
        else {
          assert(false);
        }

        // declare and load partial graphs
        G_proc_t partial_G_proc(0);
        GP_EdgeNameMap partial_GP_EdgeName = boost::get(boost::edge_name, partial_G_proc);
        GP_VertexNameMap partial_GP_VertexName = boost::get(boost::vertex_name, partial_G_proc);

        G_var_t partial_G_var(0);
        GV_EdgeWeightMap partial_GV_EdgeWeight = boost::get(boost::edge_weight, partial_G_var);
        GV_VertexNameMap partial_GV_VertexName = boost::get(boost::vertex_name, partial_G_var);
        GV_VertexKindMap partial_GV_VertexKind = boost::get(vertex_kind_t(), partial_G_var);

        inFile.open(WORKING_DIR + "temp_graph/" + fileName + "_G_proc.graph", ios::binary);
        if ( inFile.is_open() ){
          boost::archive::binary_iarchive readArchive(inFile);
          readArchive >> partial_G_proc;
          inFile.close();
        }
        else{
          assert(false);
        }
        inFile.open(WORKING_DIR + "temp_graph/" + fileName + "_G_var.graph", ios::binary);
        if ( inFile.is_open() ){
          boost::archive::binary_iarchive readArchive(inFile);
          readArchive >> partial_G_var;
          inFile.close();
        }
        else{
          assert(false);
        }

        // merge all vertices from partial graphs into merged graphs
        vector<G_proc_t::vertex_descriptor> p2m_GP_vertexDMap(partial_GP_vertexMap.size());
        for ( const auto &x : partial_GP_vertexMap ) {

          string vertexName = x.first;
          G_proc_t::vertex_descriptor partial_GP_vertexD = boost::vertex(x.second, partial_G_proc);

          // if this vertex has already been merged in, get its merged descriptor
          G_proc_t::vertex_descriptor merged_GP_vertexD;
          if ( merged_GP_vertexMap.count(vertexName) > 0 ) {
            merged_GP_vertexD = boost::vertex(merged_GP_vertexMap[vertexName], merged_G_proc);
          }

          // otherwise, create a new vertex in the merged graph
          else{
            merged_GP_vertexD = boost::add_vertex(merged_G_proc);
            merged_GP_vertexMap[vertexName] = merged_GP_vertexD;

            // copy over all internal properties from the vertex in the partial graph
            boost::put(merged_GP_VertexName, merged_GP_vertexD, boost::get(partial_GP_VertexName, partial_GP_vertexD));
          }

          p2m_GP_vertexDMap[partial_GP_vertexD] = merged_GP_vertexD;
        }
        vector<G_var_t::vertex_descriptor> p2m_GV_vertexDMap(partial_GV_vertexMap.size());
        for ( const auto&x : partial_GV_vertexMap ) {

          string vertexName = x.first;
          G_var_t::vertex_descriptor partial_GV_vertexD = x.second;

          // if this vertex has already been merged in, get its merged descriptor
          G_var_t::vertex_descriptor merged_GV_vertexD;
          if ( merged_GV_vertexMap.count(vertexName) > 0 ) {
            merged_GV_vertexD = boost::vertex(merged_GV_vertexMap[vertexName], merged_G_var);
          }

          // otherwise, create a new vertex in the merged graph
          else{
            merged_GV_vertexD = boost::add_vertex(merged_G_var);
            merged_GV_vertexMap[vertexName] = merged_GV_vertexD;

            // copy over all internal properties from the vertex in the partial graph
            boost::put(merged_GV_VertexName, merged_GV_vertexD, boost::get(partial_GV_VertexName, partial_GV_vertexD));
            boost::put(merged_GV_VertexKind, merged_GV_vertexD, boost::get(partial_GV_VertexKind, partial_GV_vertexD));
          }

          p2m_GV_vertexDMap[partial_GV_vertexD] = merged_GV_vertexD;
        }

        for ( auto& keyValPair : partial_GP_callMap ){
          merged_GP_callMap[keyValPair.first] = p2m_GP_vertexDMap[keyValPair.second];
        }

        // merge all edges from partial G_proc graph into merged G_proc graph
        auto partial_GP_edge_range = boost::edges(partial_G_proc);
        for (auto it = partial_GP_edge_range.first; it != partial_GP_edge_range.second; ++it) {

          G_proc_t::edge_descriptor partial_GP_edgeD = *it;
          G_proc_t::edge_descriptor merged_GP_edgeD;
          bool success;
          boost::tie(merged_GP_edgeD, success) = boost::add_edge(p2m_GP_vertexDMap[boost::source(partial_GP_edgeD, partial_G_proc)], p2m_GP_vertexDMap[boost::target(partial_GP_edgeD, partial_G_proc)], merged_G_proc);
          assert(success);

          // copy over all internal properties from the edge in the partial graph
          boost::put(merged_GP_EdgeName, merged_GP_edgeD, boost::get(partial_GP_EdgeName, partial_GP_edgeD));
        }

        // merge all edges from partial G_var graph into merged G_var graph
        auto partial_GV_edge_range = boost::edges(partial_G_var);
        for (auto it = partial_GV_edge_range.first; it != partial_GV_edge_range.second; ++it) {

          G_var_t::edge_descriptor partial_GV_edgeD = *it;
          G_var_t::edge_descriptor merged_GV_edgeD;
          bool success;
          boost::tie(merged_GV_edgeD, success) = boost::add_edge(p2m_GV_vertexDMap[boost::source(partial_GV_edgeD, partial_G_var)], p2m_GV_vertexDMap[boost::target(partial_GV_edgeD, partial_G_var)], merged_G_var);
          if (success) {
            // copy over all internal properties from the edge in the partial graph
            boost::put(merged_GV_EdgeWeight, merged_GV_edgeD, boost::get(partial_GV_EdgeWeight, partial_GV_edgeD));            
          }
          else {
            double weight = boost::get(merged_GV_EdgeWeight, merged_GV_edgeD);
            boost:put(merged_GV_EdgeWeight, merged_GV_edgeD, weight + boost::get(partial_GV_EdgeWeight, partial_GV_edgeD));
          }
        }
      }

      #if VERBOSE >= 3
      cout << "start printing edge weights" << endl;

      auto merged_GV_edge_range = boost::edges(merged_G_var);
      for (auto it = merged_GV_edge_range.first; it != merged_GV_edge_range.second; ++it) {
        G_var_t::edge_descriptor merged_GV_edgeD = *it;
        double weight = boost::get(merged_GV_EdgeWeight, merged_GV_edgeD);
        int sourceIndex = boost::source(merged_GV_edgeD, merged_G_var);
        int targetIndex = boost::target(merged_GV_edgeD, merged_G_var);
        string sourceName = "";
        string targetName = "";
        for (auto it = merged_GV_vertexMap.begin(); it != merged_GV_vertexMap.end(); ++it) {
          if (it->second == sourceIndex)
            sourceName = it->first;
          if (it->second == targetIndex)
            targetName = it->first;
        }
        cout << sourceName << "," << targetName << ":" << boost::get(merged_GV_EdgeWeight, merged_GV_edgeD) << endl;
      }
      #endif

      // combine all discovered interprocedural variable bindings into
      // single text file to be read in by the python framework
      vector<string> outputFileLines;
      for (auto& fileName : fileNames) {
        std::string inputBuffer;
        inFile.open(WORKING_DIR + "temp_graph/" + fileName + "_inter_bound_variables.txt");
        if ( inFile.is_open() ){
          while (getline(inFile, inputBuffer)){
            outputFileLines.push_back(inputBuffer);
          }
          inFile.close();
        }
        else {
          assert(false);
        }
      }
      outFile.open(WORKING_DIR + "__inter_bound_variables.txt", ios_base::app);
      if (outFile.is_open() ) {
        for (auto& line : outputFileLines) {
          outFile << line << endl;
        }
        outFile.close();
      }
      outputFileLines.clear();

      // combine all discovered constant variables into
      // single text file to be read in by the python framework
      for (auto& fileName : fileNames) {
        std::string inputBuffer;
        inFile.open(WORKING_DIR + "temp_graph/" + fileName + "_constant_list.txt");
        if ( inFile.is_open() ){
          while (getline(inFile, inputBuffer)){
            outputFileLines.push_back(inputBuffer);
          }
          inFile.close();
        }
        else {
          assert(false);
        }
      }
      outFile.open(WORKING_DIR + "constant_list.txt", ios_base::app);
      if (outFile.is_open() ) {
        for (auto& line : outputFileLines) {
          outFile << line << endl;
        }
        outFile.close();
      }
      outputFileLines.clear();

      // save the merged G_proc and its vertexMap
      outFile.open(WORKING_DIR + "G_proc.graph", ios_base::binary);
      if ( outFile.is_open() ){
        boost::archive::binary_oarchive writeArchive(outFile);
        writeArchive << merged_G_proc;
        outFile.close();
      }
      else{
        assert(false);
      }
      outFile.open(WORKING_DIR + "GP_vertexMap.map", ios_base::binary);
      if ( outFile.is_open() ){
        boost::archive::binary_oarchive writeArchive(outFile);
        writeArchive << merged_GP_vertexMap;
        outFile.close();
      }
      else{
        assert(false);
      }
      outFile.open(WORKING_DIR + "GP_callMap.map", ios_base::binary);
      if ( outFile.is_open() ){
        boost::archive::binary_oarchive writeArchive(outFile);
        writeArchive << merged_GP_callMap;
        outFile.close();
      }
      else{
        assert(false);
      }

      // write out a dot file of merged G_proc to be read in by the
      // Python framework
      outFile.open("prose_logs/__G_proc.dot", ios_base::trunc);
      if ( outFile.is_open() ){
        boost::write_graphviz(outFile, merged_G_proc, boost::make_label_writer(merged_GP_VertexName), boost::make_label_writer(merged_GP_EdgeName));
        outFile.close();
      }
      else{
        assert(false);
      }

      // write out a dot file of merged G_var to be read in by the
      // python framework
      outFile.open("prose_logs/__G_var.dot", ios_base::trunc);
      if ( outFile.is_open() ){
        boost::write_graphviz(outFile, merged_G_var, boost::make_label_writer(merged_GV_VertexName), boost::make_label_writer(merged_GV_EdgeWeight));
        outFile.close();
      }
      else{
        assert(false);
      }

    }
}; // end LinkGraph


 class ApplyConfiguration : public Rose::PluginAction {
  public:
    ApplyConfiguration() {}
    ~ApplyConfiguration() {}

    string WORKING_DIR;
    string CONFIGURATION_DIR;
    map<string, int> targetConfig;
    vector<string> filePaths;
    map<int, SgType*> replacementKinds;
    vector<string> ignoreScopes;

    // declare necessary graph objects
    G_proc_t G_proc;
    GP_EdgeNameMap GP_EdgeName;
    GP_VertexNameMap GP_VertexName;
    map<string, int> GP_vertexMap;
    map<string, int> GP_callMap;

    unsigned long long typeConversions = 0;
    unsigned long long floatingPointOp = 0;
    unsigned long long opInLowPrecision = 0;
    unsigned long long divOps = 0;
    unsigned long long intrinsicOps = 0;

    // command line args are:
    // 1. path to the working directory
    // 2. path to the log directory for the configuration being applied
    bool ParseArgs(const std::vector<std::string> &args){
      assert( args.size() == 2 );

      fstream inFile;
      string inputBuffer;
      vector<string> splitBuffer;
      vector<string> transformations;

      WORKING_DIR = args[0];
      CONFIGURATION_DIR = args[1];

      if ( WORKING_DIR[WORKING_DIR.size() - 1] != '/' ){
        WORKING_DIR = WORKING_DIR + '/';
      }
      if ( CONFIGURATION_DIR[CONFIGURATION_DIR.size() - 1] != '/' ){
        CONFIGURATION_DIR = CONFIGURATION_DIR + '/';
      }


      // open config file
      // split each line in the config file into the name of the var and the kind it is to be assigned; store them in targetConfig
      inFile.open(CONFIGURATION_DIR + "config", ios::in);
      if ( inFile.is_open() ){
        while (getline(inFile, inputBuffer)){
          if (inputBuffer != ""){
            boost::split(splitBuffer, inputBuffer, boost::is_any_of(","));
            string varName = splitBuffer[0];
            varName.erase(remove_if(varName.begin(), varName.end(), ::isspace), varName.end()); // remove whitespace
            targetConfig[varName] = atoi(splitBuffer[1].c_str());
          }
        }
      }
      else{
        assert(false);
      }
      inFile.close();

      // read in file paths
      inFile.open(WORKING_DIR + "prose_workspace/__target_files.txt", ios::in);
      if (inFile.is_open()){
        while (getline(inFile, inputBuffer)){
          filePaths.push_back(inputBuffer);
        }
      }
      else{
        assert(false);
      }
      inFile.close();

      // read in ignoreScopes
      inFile.open("prose_workspace/ignore_scopes.txt", ios::in);
      if (inFile.is_open()){
        while (getline(inFile, inputBuffer)){
          ignoreScopes.push_back(inputBuffer);
        }
      }
      else{
        assert(false);
      }
      inFile.close();

      return true;
    }

    void process (SgProject* n) {

      #if VERBOSE >= 1
      cout << "Start Action: Apply Configuration" << endl;
      #endif

      // instantiate necessary graph objects
      GP_EdgeName = boost::get(boost::edge_name, G_proc);
      GP_VertexName = boost::get(boost::vertex_name, G_proc);
      load_graph();

      // parse all given files into the project AST
      unordered_set<SgSourceFile*> sourceFilesToBeUnparsed = parse_source_files_into_AST(n, filePaths, CONFIGURATION_DIR);

      // construct replacement kinds
      replacementKinds[4] = SgTypeFloat::createType(SageBuilder::buildIntVal(4));
      replacementKinds[8] = SgTypeFloat::createType(SageBuilder::buildIntVal(8));
      replacementKinds[10] = SgTypeFloat::createType(SageBuilder::buildIntVal(10));
      replacementKinds[16] = SgTypeFloat::createType(SageBuilder::buildIntVal(16));

      for ( const auto& sourceFile : sourceFilesToBeUnparsed ){
        preprocess(sourceFile);
        apply_configuration(sourceFile);
      }

      // data structure to keep track of what wrappers have already
      // been generated to avoid redundancy
      set<string> wrapperNames;

      for ( const auto& sourceFile : sourceFilesToBeUnparsed ){
        wrapper_fix(sourceFile, wrapperNames);
        literal_fix(sourceFile);
        sign_intrinsic_fix(sourceFile);
      }

      for ( const auto& sourceFile : sourceFilesToBeUnparsed ){

        #if VERBOSE >= 2
        cout << "\t unparsing " << sourceFile->get_unparse_output_filename() << endl;
        #endif

        sourceFile->unparse();
      }

      exit(0);
    } // end process()

  private:


    void literal_fix( const SgSourceFile* sourceFile ){

      #if VERBOSE >= 11
      cout << "[BEGIN] literal fix for " << sourceFile->getFileName() << endl;
      #endif

      infer_kinds_of_literals(NodeQuery::querySubTree(sourceFile->get_globalScope(), V_SgFloatVal), false);
    }


    // fixing function calls where there is a kind mismatch between real parameters
    // currently, only concerned with the "sign" intrinsic 
    void sign_intrinsic_fix( const SgSourceFile* sourceFile ){

      #if VERBOSE >= 11
      cout << "[BEGIN] function call fix for " << sourceFile->getFileName() << endl;
      #endif

      auto functionCalls = NodeQuery::querySubTree(sourceFile->get_globalScope(), V_SgFunctionCallExp);

      for (const auto& y : functionCalls ) {
        if ( SgFunctionCallExp* functionCall = isSgFunctionCallExp(y) ){

          // force sign function argument types to match if they are floats
          if ( boost::algorithm::to_lower_copy(functionCall->getAssociatedFunctionSymbol()->get_name().getString()) == "sign" ) {
            vector<SgExpression*> givenArguments = functionCall->get_args()->get_expressions();

            SgType* type1 = get_expression_type(givenArguments[0]);
            SgType* type2 = get_expression_type(givenArguments[1]);

            if (type1->class_name() == "SgTypeFloat" && type2->class_name() == "SgTypeFloat") {
              SgScopeStatement* fromScope = SageInterface::getEnclosingScope(functionCall);
              while ( !(isSgFunctionDefinition(fromScope) || isSgClassDefinition(fromScope)) ){
                fromScope = SageInterface::getEnclosingScope(fromScope);
              }

              vector<SgExpression*> realConversionArgList1;
              realConversionArgList1.push_back(SageInterface::copyExpression(givenArguments[0]));
              realConversionArgList1.push_back(SageBuilder::buildIntVal(DEFAULT_KIND));

              vector<SgExpression*> realConversionArgList2;
              realConversionArgList2.push_back(SageInterface::copyExpression(givenArguments[1]));
              realConversionArgList2.push_back(SageBuilder::buildIntVal(DEFAULT_KIND));

              SgFunctionCallExp* realConversionCall1 = SageBuilder::buildFunctionCallExp("REAL", type1, SageBuilder::buildExprListExp(realConversionArgList1), fromScope);
              SageInterface::replaceExpression(givenArguments[0], realConversionCall1, false);
              SgFunctionCallExp* realConversionCall2 = SageBuilder::buildFunctionCallExp("REAL", type2, SageBuilder::buildExprListExp(realConversionArgList2), fromScope);
              SageInterface::replaceExpression(givenArguments[1], realConversionCall2, false);
            }
          }
        }
      }
    }


    void wrapper_fix( const SgSourceFile* sourceFile, set<string>& wrapperNames ){

      #if VERBOSE >= 11
      cout << "[BEGIN] wrapper fix for " << sourceFile->getFileName() << endl;
      #endif

      auto functionCalls = NodeQuery::querySubTree(sourceFile->get_globalScope(), V_SgFunctionCallExp);

      for ( int i = 0; i < functionCalls.size(); i++ ){
        if ( SgFunctionCallExp* functionCall = isSgFunctionCallExp(functionCalls[i]) ){
          if ( !functionCall->get_containsTransformation() ){

            #if VERBOSE >= 11
            cout << "\t checking " << functionCall->unparseToString();
            #endif

            string scopedCallName = get_scoped_name(functionCall);

            #if VERBOSE >= 11
            cout << "( " << scopedCallName << " )" << endl;
            #endif

            if ( GP_callMap.count(scopedCallName) > 0 ){
              G_proc_t::vertex_descriptor GP_toScope_vertexD = boost::vertex(GP_callMap[scopedCallName], G_proc);
              string toScopeName = boost::get(GP_VertexName, GP_toScope_vertexD);

              #if VERBOSE >= 11
              cout << "\t\t => " << toScopeName << endl;
              #endif

              bool ignore = false;
              for ( string ignoreScope : ignoreScopes ){
                if ( boost::algorithm::to_lower_copy(toScopeName) == ignoreScope ){
                  ignore = true;
                  break;
                }
              }
              if ( ignore ){

                #if VERBOSE >= 12
                cout << "\t\t\t blacklisted scope! skipping" << endl;
                #endif

                continue;
              }

              SgFunctionDefinition* toScope = isSgFunctionDefinition(get_specified_scope(SageInterface::getProject(), toScopeName));

              if ( toScope != NULL ){
                if ( SgProcedureHeaderStatement* calleeProcHeader = isSgProcedureHeaderStatement(toScope->get_declaration()) ){

                  #if VERBOSE >= 11
                  cout << "\t\t calleeProcHeader: " << calleeProcHeader->unparseToString() << endl;
                  #endif

                  // now that we've found the callee procedure, check if
                  // a wrapper is necessary
                  map<string,pair<int, int>> argKindConfigs;
                  bool mismatch;
                  check_args(functionCall, calleeProcHeader, &mismatch, &argKindConfigs);
                  if ( mismatch ){

                    #if VERBOSE >= 11
                    cout << "\t\t\t MISMATCH" << endl;
                    #endif

                    wrap_procedure( functionCall, calleeProcHeader, argKindConfigs, wrapperNames );
                  }
                }
              }
            }
            else {
              #if VERBOSE >= 11
              cout << "\t\t toScope obscured, possibly an intrinsic function. Inferring kinds of literals." << endl;
              #endif

              // here, we infer the kinds of float literals present as
              // arguments to intrinsic functions in order to ensure compilability
              infer_kinds_of_literals(NodeQuery::querySubTree(functionCall->get_args(), V_SgFloatVal), true);
            }
          }
        }
      }

    #if VERBOSE >= 11
    cout << "[END] wrapper fix for " << sourceFile->getFileName() << endl;
    #endif

    }


    void wrap_procedure( SgFunctionCallExp* functionCall,
                        SgProcedureHeaderStatement* calleeProcHeader,
                        map<string,pair<int, int>>& argKindConfigs,
                        set<string>& wrapperNames ){

      // if wrapper name exceeds fortran maximum of 63 characters, tag it with _nameTrunc and a random integer
      // to maintain uniqueness
      // also add an optional tag for any procedures whose specific names
      // are different from their called name due to generic interfaces
      hash<string> hash_string;
      string calleeName = functionCall->getAssociatedFunctionSymbol()->get_name();
      string optionalTag = "";
      if (boost::algorithm::to_lower_copy(calleeProcHeader->get_name().getString()) != boost::algorithm::to_lower_copy(calleeName)){
        optionalTag = to_string(hash_string(calleeProcHeader->get_name().getString())%1000);
      }
      string expectedKindConfigString = generate_expected_kind_config_string(calleeProcHeader);
      string givenKindConfigString = generate_given_kind_config_string(functionCall, calleeProcHeader);
      string wrapperName = calleeName + "_wrapper_";
      if ( optionalTag != "" ){
        wrapperName = wrapperName + "id" + optionalTag + "_";
      }
      wrapperName = wrapperName + givenKindConfigString + "_to_" + expectedKindConfigString;
      if ( wrapperName.size() > 63 ){
        wrapperName = calleeName + "_wrap_" + to_string(hash_string(givenKindConfigString)%1000);
      }

      // get destination scope
      bool NO_CONTAINING_MODULE = false;
      SgScopeStatement* destinationScope = SageInterface::getEnclosingClassDefinition(calleeProcHeader);
      if ( destinationScope == NULL ){
        NO_CONTAINING_MODULE = true;
        destinationScope = SageInterface::getGlobalScope(calleeProcHeader);
      }
      assert(destinationScope != NULL);

      // if it hasn't already been generated, generate a new wrapper
      if ( wrapperNames.count(wrapperName) == 0 ){
        wrapperNames.insert(wrapperName);

        // add a "CONTAINS" statement if there is not one already
        auto containsStatements = NodeQuery::querySubTree(destinationScope, V_SgContainsStatement);
        if ( (containsStatements.size() == 0) && (!NO_CONTAINING_MODULE) ){
          SgContainsStatement* containsStatement = new SgContainsStatement(Sg_File_Info::generateDefaultFileInfoForTransformationNode());
          containsStatement->set_firstNondefiningDeclaration(containsStatement);
          containsStatement->set_definingDeclaration(containsStatement);
          SageInterface::appendStatement(containsStatement, destinationScope);
        }

        #if VERBOSE >= 11
        cout << "\n[BEGIN] building wrapper " << wrapperName << endl;
        #endif

        //////////////////////////////////////////////////// BEGIN LOCAL VAR DEDUCTION
        #if VERBOSE >= 11
        cout << "\t deducing local variables..." << endl;
        #endif

        // first value in the pair is the input parameters to the wrapper procedure
        // second value is the temporary variable the input parameters is cast to before the function call (if any)
        map<string,pair<SgInitializedName*, SgInitializedName*>> localVarInitNames;

        vector<SgInitializedName*> dummyVars_paramListOrder = calleeProcHeader->get_parameterList()->get_args();
        vector<SgInitializedName*> dummyVars_declListOrder = get_dummy_var_decl_order(calleeProcHeader);

        // go through all of the dummyVars_paramListOrder
        for ( const auto& originalArg : dummyVars_paramListOrder ){

          string varName = originalArg->get_name().str();
          boost::algorithm::to_lower(varName);

          if ( argKindConfigs[varName].first == OPTIONAL_ARG_NOT_PROVIDED ){
            continue;
          }

          #if VERBOSE >= 11
          cout << "\t\t processing " << varName << endl;
          #endif

          SgInitializedName* var = SageBuilder::buildInitializedName(varName, originalArg->get_type());

          SgInitializedName* tempVar = isSgInitializedName(calleeProcHeader); // intializes to NULL
          string tempVarName = varName + "__temp";
          SgType* tempVarType;

          // if the kind has changed from what was expected OR the dummy argument is optional a cast must be made.
          // First, a pair of variable declarations must be created: one for the dummy argument and one for the
          // temporary variable it is cast to before being fed as an argument to the function call.
          if ( argKindConfigs[varName].first != argKindConfigs[varName].second){
            // In the case of an optional dummy argument, the temporary variable is a pointer
            if ( originalArg->get_declaration()->get_declarationModifier().get_typeModifier().isOptional() ){

              // if the original arg is already a pointer, we can just
              // reuse the type
              if ( isSgPointerType(originalArg->get_type()) ){
                tempVarType = originalArg->get_type();
              }

              // if the original arg does not contain internal types, the temp variable can just be a
              // pointer to the original type
              else if ( !originalArg->get_type()->containsInternalTypes() ){
                tempVarType = SageBuilder::buildPointerType(originalArg->get_type());
              }

              // however, if the original arg does contain internal types, we make deep copies of each type node
              // in the type chain and, for each array, ensure that it is a deferred shape array type
              else{

                // get original base type, i.e., the last type in the type chain
                tempVarType = originalArg->get_type()->findBaseType();

                // iterate through internal types in reverse, making deep copies of each type node and
                // constructing a cloned type chain
                auto internalTypes = originalArg->get_type()->getInternalTypes();
                for ( auto internalTypes_it = internalTypes.rbegin() + 1; internalTypes_it != internalTypes.rend(); ++internalTypes_it ){

                  SgType* deepCopiedType = isSgType(SageInterface::deepCopyNode(*internalTypes_it));

                  // if it's an array type, make it a deferred shape array before adding it to the chain
                  if ( SgArrayType* arrayType = isSgArrayType(deepCopiedType) ){

                    auto dimExpressionList = arrayType->get_dim_info()->get_expressions();
                    for ( auto dimExpressionList_it = dimExpressionList.begin(); dimExpressionList_it != dimExpressionList.end(); ++dimExpressionList_it ){
                      SgColonShapeExp* colonExpression = new SgColonShapeExp(Sg_File_Info::generateDefaultFileInfoForTransformationNode());
                      colonExpression->set_endOfConstruct(Sg_File_Info::generateDefaultFileInfoForTransformationNode());
                      SageInterface::replaceExpression(*dimExpressionList_it, colonExpression);
                    }
                    arrayType->set_base_type(tempVarType);
                  }

                  // else, if it's a pointer type, we can just add it to the chain as-is
                  else if ( SgPointerType* pointerType = isSgPointerType(deepCopiedType) ){
                    pointerType->set_base_type(tempVarType);
                  }

                  tempVarType = deepCopiedType;
                }

                // seems roundabout but if not done in this way, some behind-the-scenes pointer changes done
                // by SageBuilder messes up the rank of the resulting underlying array type
                SgPointerType* pointerType = SageBuilder::buildPointerType(tempVarType);
                pointerType->set_base_type(tempVarType);
                tempVarType = pointerType;
              }
            }

            // otherwise, if it's not optional, we can just reuse the existing type...
            else{

              tempVarType = originalArg->get_type();

              // ... unless it's an assumed shape array in which case we
              // must make the temp var a deferred shape array
              if ( SgArrayType* arr = isSgArrayType(tempVarType) ){

                vector<SgExpression*> dims = arr->get_dim_info()->get_expressions();
                if ( dims.size() > 0 ){
                  if ( boost::algorithm::ends_with(dims.front()->unparseToString(), ":") ){

                    // get original base type, i.e., the last type in the type chain
                    tempVarType = originalArg->get_type()->findBaseType();

                    // iterate through internal types in reverse, making deep copies of each type node and
                    // constructing a cloned type chain
                    auto internalTypes = originalArg->get_type()->getInternalTypes();
                    for ( auto internalTypes_it = internalTypes.rbegin() + 1; internalTypes_it != internalTypes.rend(); ++internalTypes_it ){

                      SgType* deepCopiedType = isSgType(SageInterface::deepCopyNode(*internalTypes_it));

                      // if it's an array type, make it a deferred shape array before adding it to the chain
                      if ( SgArrayType* arrayType = isSgArrayType(deepCopiedType) ){

                        auto dimExpressionList = arrayType->get_dim_info()->get_expressions();
                        for ( auto dimExpressionList_it = dimExpressionList.begin(); dimExpressionList_it != dimExpressionList.end(); ++dimExpressionList_it ){
                          SgColonShapeExp* colonExpression = new SgColonShapeExp(Sg_File_Info::generateDefaultFileInfoForTransformationNode());
                          colonExpression->set_endOfConstruct(Sg_File_Info::generateDefaultFileInfoForTransformationNode());
                          SageInterface::replaceExpression(*dimExpressionList_it, colonExpression);
                        }
                        arrayType->set_base_type(tempVarType);
                      }

                      // else, if it's a pointer type, we can just add it to the chain as-is
                      else if ( SgPointerType* pointerType = isSgPointerType(deepCopiedType) ){
                        pointerType->set_base_type(tempVarType);
                      }

                      tempVarType = deepCopiedType;
                    }
                  }
                }
              }
            }

            // build the temporary variable
            tempVar = SageBuilder::buildInitializedName(tempVarName, tempVarType);
          }

          if ( tempVar != NULL ){
            set_real_kind(var, replacementKinds[argKindConfigs[varName].first]);
            set_real_kind(tempVar, replacementKinds[argKindConfigs[varName].second]);
          }

          localVarInitNames[varName] = make_pair(var, tempVar);
        }

        SgInitializedName* returnVar = isSgInitializedName(calleeProcHeader); // initializes to NULL
        if ( calleeProcHeader->isFunction() ){
          #if VERBOSE >= 11
          cout << "\t\t because this is a function-type procedure, adding a return value..." << endl;
          #endif

          returnVar = SageBuilder::buildInitializedName(wrapperName, calleeProcHeader->get_result_name()->get_type());
        }

        //////////////////////////////////////////////////// END LOCAL VAR DEDUCTION

        //////////////////////////////////////////////////// BEGIN CREATION OF WRAPPER PROCEDURE HEADER
        #if VERBOSE >= 11
        cout << "\t creating procedure header..." << endl;
        #endif

        // construct function parameter list for procedure wrapper
        SgFunctionParameterList* parameterList = SageBuilder::buildFunctionParameterList();
        for ( const auto &originalArg : dummyVars_paramListOrder ){

          string varName = originalArg->get_name();
          boost::algorithm::to_lower(varName);

          if ( argKindConfigs[varName].first == OPTIONAL_ARG_NOT_PROVIDED ){
            continue;
          }
          else{
            SageInterface::appendArg(parameterList, localVarInitNames[varName].first);
          }
        }

        // create header statement for procedure wrapper
        SgProcedureHeaderStatement* wrapperProcHeader = SageBuilder::buildProcedureHeaderStatement(wrapperName, calleeProcHeader->get_orig_return_type(), parameterList, calleeProcHeader->get_subprogram_kind(), destinationScope);
        if ( calleeProcHeader->get_functionModifier().isElemental() ){
          wrapperProcHeader->get_functionModifier().setElemental();
        }
        if ( calleeProcHeader->get_functionModifier().isPure() ){
          wrapperProcHeader->get_functionModifier().setPure();
        }

        //////////////////////////////////////////////////// END CREATION OF WRAPPER PROCEDURE HEADER


        //////////////////////////////////////////////////// BEGIN POPULATION OF WRAPPER PROCEDURE BODY

        SgBasicBlock* funcBody = wrapperProcHeader->get_definition()->get_body();

        #if VERBOSE >= 11
        cout << "\t populating procedure body with variable declarations..." << endl;
        #endif

        // add any use statements from the original procedure in case they are needed for the variable declarations
        for ( const auto& useStatement : NodeQuery::querySubTree(calleeProcHeader, V_SgUseStatement) ){
          SageInterface::appendStatement(isSgUseStatement(useStatement), funcBody);
        }

        for ( const auto& originalArg : dummyVars_declListOrder ){

          string varName = originalArg->get_name();
          boost::algorithm::to_lower(varName);

          // skip if the arg was optional and not provided in the function call
          if ( argKindConfigs[varName].first == OPTIONAL_ARG_NOT_PROVIDED ){
            continue;
          }

          // append dummy variable declaration
          SgInitializedName* var = localVarInitNames[varName].first;
          SgVariableDeclaration* varDecl = SageBuilder::buildVariableDeclaration(var->get_name(), var->get_type(), var->get_initializer(), isSgScopeStatement(funcBody));
          varDecl->get_declarationModifier().get_typeModifier().set_modifierVector(originalArg->get_declaration()->get_declarationModifier().get_typeModifier().get_modifierVector());
          varDecl->get_declarationModifier().get_accessModifier().setUndefined();
          if ( !is_underlying_real_type(originalArg) ){
            varDecl->get_declarationModifier().get_typeModifier().unsetOptional();
          }
          SageInterface::appendStatement(varDecl, funcBody);

          // append declaration for temporary variable if it exists
          var = localVarInitNames[varName].second;
          if ( var != NULL ){

            // if it is a pointer, we need to initialize to NULL()
            if ( isSgPointerType(var->get_type()) ){
              var->set_initializer(SageBuilder::buildAssignInitializer(SageBuilder::buildIntVal()));
            }

            varDecl = SageBuilder::buildVariableDeclaration(var->get_name(), var->get_type(), var->get_initializer(), isSgScopeStatement(funcBody));
            varDecl->get_declarationModifier().get_typeModifier().set_modifierVector(originalArg->get_declaration()->get_declarationModifier().get_typeModifier().get_modifierVector());
            varDecl->get_declarationModifier().get_accessModifier().setUndefined();
            varDecl->get_declarationModifier().get_typeModifier().unsetIntent_in();
            varDecl->get_declarationModifier().get_typeModifier().unsetIntent_out();
            varDecl->get_declarationModifier().get_typeModifier().unsetIntent_inout();
            varDecl->get_declarationModifier().get_typeModifier().unsetOptional();
            SageInterface::appendStatement(varDecl, funcBody);

            // if it is a deferred-shape array, we need to mark temp as allocatable...
            if ( SgArrayType* arrayVar = isSgArrayType(var->get_type()) ){
              vector<SgExpression*> dimList = arrayVar->get_dim_info()->get_expressions();
              if ( dimList.size() > 0 ){
                if ( dimList.front()->unparseToString() == ":" ){

                  // set the declaration for the temp array to be allocatable
                  varDecl->get_declarationModifier().get_typeModifier().setAllocatable();
                }
              }
            }
          }
        }

        // append variable declaration for return value if there is one
        if ( returnVar != NULL ){
          wrapperProcHeader->set_result_name(returnVar);
          SgVariableDeclaration* returnVarDecl = SageBuilder::buildVariableDeclaration(returnVar->get_name(), returnVar->get_type(), returnVar->get_initializer(), isSgScopeStatement(funcBody));
          
          // needed to remove default "PUBLIC" modifier
          returnVarDecl->get_declarationModifier().get_accessModifier().setUndefined();
          
          // needed so that we don't try to unparse the return type into the procedure header statement which fails for array types
          returnVar->set_definition(returnVarDecl);

          SageInterface::appendStatement(returnVarDecl, funcBody);
        }

        #if VERBOSE >= 11
        cout << "\t constructing all the pre-function-call statements..."<< endl;
        #endif

        // construct all the pre-function-call statements
        for ( const auto &x : localVarInitNames ){

          SgInitializedName* varFrom = x.second.first;
          SgInitializedName* varTo = x.second.second;

          // continue if there is no cast required
          if ( varTo == NULL ){
            continue;
          }

          #if VERBOSE >= 12
          cout << "\t\t processing... " << varFrom->get_name() << ", " << varTo->get_name() << endl;
          #endif

          SgIfStmt* ifStatement = NULL;
          SgExprStatement* castViaAssignmentStatement = NULL;

          // if varFrom is not intent(out), we must include a pre-function-call casting statement
          if ( (!varFrom->get_declaration()->get_declarationModifier().get_typeModifier().isIntent_out()) ){
            SgVarRefExp* lhs = SageBuilder::buildVarRefExp(varTo, funcBody);
            SgVarRefExp* rhs = SageBuilder::buildVarRefExp(varFrom, funcBody);
            SgAssignOp* castOp = SageBuilder::buildAssignOp(lhs, rhs);
            castViaAssignmentStatement = SageBuilder::buildExprStatement(castOp);
          }

          // if varTo is a pointer, we must allocate space for it
          if ( SgPointerType* varToPtr = isSgPointerType(varTo->get_type()) ){

            // BEGIN BUILD ALLOCATE STATEMENT
            // if the pointer is to an array, the argument to the allocate statement must be of
            // of the form tempVar(SIZE(originalVar, rank_1),SIZE(originalVar, rank_2),...) for each rank_i in rank
            vector<SgExpression*> allocateArgList;
            SgArrayType* underlyingArrayType = isSgArrayType(varToPtr->get_base_type());
            if ( underlyingArrayType != NULL ){

              vector<SgExpression*> arrayRefArgList;
              int rank = underlyingArrayType->get_rank();

              for ( int i = 1; i <= rank; ++i ){
                vector<SgExpression*> shapeArgList;
                shapeArgList.push_back(SageBuilder::buildVarRefExp(varFrom, funcBody));
                shapeArgList.push_back(SageBuilder::buildIntVal(i));
                SgFunctionCallExp* shapeCall = SageBuilder::buildFunctionCallExp("SIZE", SageBuilder::buildIntType(), SageBuilder::buildExprListExp(shapeArgList), funcBody);
                arrayRefArgList.push_back(shapeCall);
              }

              allocateArgList.push_back(SageBuilder::buildPntrArrRefExp(SageBuilder::buildVarRefExp(varTo, funcBody), SageBuilder::buildExprListExp(arrayRefArgList)));
            }
            // otherwise, we can just allocate space for the single var
            else{
              allocateArgList.push_back(SageBuilder::buildVarRefExp(varTo, funcBody));
            }

            // build allocate statement
            SgAllocateStatement* allocateStatement = new SgAllocateStatement(Sg_File_Info::generateDefaultFileInfoForTransformationNode());
            allocateStatement->set_parent(funcBody);
            allocateStatement->set_expr_list(SageBuilder::buildExprListExp(allocateArgList));
            // END BUILD ALLOCATE STATEMENT

            // BEGIN BUILD conditional statement
            // a call to the ASSOCIATED intrinsic to check for the association status of varFrom if varFrom is also a pointer
            if ( isSgPointerType(varFrom->get_type()) ){
              vector<SgExpression*> associatedArgList;
              associatedArgList.push_back(SageBuilder::buildVarRefExp(varFrom, funcBody));
              SgFunctionCallExp* associatedCall = SageBuilder::buildFunctionCallExp("ASSOCIATED", SageBuilder::buildBoolType(), SageBuilder::buildExprListExp(associatedArgList), funcBody);

              // build if statement
              if ( castViaAssignmentStatement != NULL ){
                ifStatement = SageBuilder::buildIfStmt(associatedCall, SageBuilder::buildBasicBlock(allocateStatement, castViaAssignmentStatement), NULL);
              }
              else {
                ifStatement = SageBuilder::buildIfStmt(associatedCall, SageBuilder::buildBasicBlock(allocateStatement), NULL);
              }
            }

            // if varFrom is also optional, we must include the if statement constructed above inside of the
            // true branch of an "if present" construct
            if ( varFrom->get_declaration()->get_declarationModifier().get_typeModifier().isOptional() ){

              // build conditional, a call to PRESENT intrinsic to check for presence of optional arg
              vector<SgExpression*> presentArgList;
              presentArgList.push_back(SageBuilder::buildVarRefExp(varFrom, funcBody));
              SgFunctionCallExp* presentCall = SageBuilder::buildFunctionCallExp("PRESENT", SageBuilder::buildBoolType(), SageBuilder::buildExprListExp(presentArgList), funcBody);

              if ( ifStatement != NULL ){
                ifStatement = SageBuilder::buildIfStmt(presentCall, SageBuilder::buildBasicBlock(ifStatement), NULL);
              }
              else if ( castViaAssignmentStatement != NULL ){
                ifStatement = SageBuilder::buildIfStmt(presentCall, SageBuilder::buildBasicBlock(allocateStatement, castViaAssignmentStatement), NULL);
              }
              else{
                ifStatement = SageBuilder::buildIfStmt(presentCall, SageBuilder::buildBasicBlock(allocateStatement), NULL);
              }
            }
          }

          if ( ifStatement != NULL ){
            SageInterface::appendStatement(ifStatement, funcBody);
          }
          else if ( castViaAssignmentStatement != NULL ){
            SageInterface::appendStatement(castViaAssignmentStatement, funcBody);
          }
          else{
            continue;
          }
        }

        #if VERBOSE >= 11
        cout << "\t constructing parameter list for the function call..." <<endl;
        #endif

        // construct parameter list for function call
        vector<SgExpression*> expressions;
        for ( const auto &originalArg : dummyVars_paramListOrder ){

          string varName = originalArg->get_name();
          boost::algorithm::to_lower(varName);

          if ( localVarInitNames.count(varName) == 0 ){
            continue;
          }

          SgInitializedName* argumentVar;
          pair<SgInitializedName*,SgInitializedName*> varPair = localVarInitNames[varName];

          // if there isn't a cast, use the original var as the argument.
          // otherwise, use the temp variable that holds the result of the cast
          if ( varPair.second == NULL ) {
            argumentVar = varPair.first;
          }
          else{
            argumentVar = varPair.second;
          }

          // build an expression referencing that variable and associate it with its keyword argument if it is optional and this isn't MAIN
          SgExpression* expr = SageBuilder::buildVarRefExp(argumentVar, funcBody);
          if ( !NO_CONTAINING_MODULE ){
            expr = SageBuilder::buildActualArgumentExpression(originalArg->get_name(), expr);
          }
          expressions.push_back(expr);
        }

        #if VERBOSE >= 11
        cout << "\t constructing function call..." << endl;
        #endif

        // construct function call to originally-named expression
        SgFunctionSymbol* calleeSymbol = SageInterface::lookupFunctionSymbolInParentScopes(calleeProcHeader->get_name(), SageInterface::getEnclosingScope(calleeProcHeader));
        assert( calleeSymbol != NULL );
        SgFunctionRefExp* funcRef = SageBuilder::buildFunctionRefExp(calleeSymbol);
        SgExprListExp* argumentList = SageBuilder::buildExprListExp(expressions);
        SgFunctionCallExp* funcCall = SageBuilder::buildFunctionCallExp(funcRef, argumentList);

        // if this procedure is a function with a return value, we need to assign to it
        SgExprStatement* newStatement;
        if ( calleeProcHeader->isFunction() ){
          SgVarRefExp* result = SageBuilder::buildVarRefExp(wrapperName, funcBody);
          SgAssignOp* assignToReturn = SageBuilder::buildAssignOp(result, funcCall);
          newStatement = SageBuilder::buildExprStatement(assignToReturn);
        }
        else{
          newStatement = SageBuilder::buildExprStatement(funcCall);
        }
        SageInterface::appendStatement(newStatement, funcBody);

        #if VERBOSE >= 11
        cout << "\t constructing all post-function-call statements..." << endl;
        #endif

        // construct all the post-function-call statements
        for ( const auto &x : localVarInitNames ){

          SgInitializedName* varTo = x.second.first;
          SgInitializedName* varFrom = x.second.second;

          // continue if there is no cast required
          if ( (varTo == NULL) || (varFrom == NULL) ){
            continue;
          }

          // possible statements...
          SgExprStatement* castViaAssignmentStatement = NULL;
          SgExprStatement* cleanupTemporaryVariableStatement = NULL;
          SgIfStmt* ifStatement = NULL;
          
          // if the original variable is not intent(in), perform post-function-call casting
          if ( !varTo->get_declaration()->get_declarationModifier().get_typeModifier().isIntent_in() ){
            SgVarRefExp* lhs = SageBuilder::buildVarRefExp(varTo, funcBody);
            SgVarRefExp* rhs = SageBuilder::buildVarRefExp(varFrom, funcBody);
            SgAssignOp* castOp = SageBuilder::buildAssignOp(lhs, rhs);
            castViaAssignmentStatement = SageBuilder::buildExprStatement(castOp);
          }

          // if the temporary variable (aka, varFrom) is a pointer, we must clean up the memory
          // (either because the original variable is also a pointer, is optional, or both)
          // then we must eventually nullify that pointer
          if ( isSgPointerType(varFrom->get_type()) ){
            if ( isSgPointerType(varTo->get_type()) ){
              SgFunctionCallExp* rhs = SageBuilder::buildFunctionCallExp("NULL", varTo->get_type(), NULL, funcBody);
              SgVarRefExp* lhs = SageBuilder::buildVarRefExp(varFrom, funcBody);
              SgPointerAssignOp* pointerAssignOp = new SgPointerAssignOp(Sg_File_Info::generateDefaultFileInfoForTransformationNode(), lhs, rhs, NULL);
              cleanupTemporaryVariableStatement = SageBuilder::buildExprStatement(pointerAssignOp);
            }
            else{
              vector<SgExpression*> argList;
              argList.push_back(SageBuilder::buildVarRefExp(varFrom, funcBody));
              SgFunctionCallExp* deallocCall = SageBuilder::buildFunctionCallExp("DEALLOCATE", SageBuilder::buildIntType(), SageBuilder::buildExprListExp(argList), funcBody);
              cleanupTemporaryVariableStatement = SageBuilder::buildExprStatement(deallocCall);
            }
          }


          // if the original variable (aka, varTo) is also a pointer, we need a block of the following form
          /*
            *  ! case of deallocation or no change in wrapped procedure
            *  if ( .not. associated(var_temp) ) then
            *        var => NULL()
            *  else
            * 
            *    ! case of reallocation in wrapped procedure
            *    if ( associated(var) ) then
            *      if (size(var) /= size(var_temp) ) then
            *        var => NULL()
            *      endif
            *    endif
            * 
            *    ! case of reallocation OR initialization in wrapped procedure
            *    if ( .not. associated(var) ) then
            *          allocate( var(SIZE(var_temp,1),SIZE(var_temp,2),SIZE(var_temp,3)) )
            *    endif
            * 
            *    ! case of reallocation OR initialization OR simple change in wrapped procedure
            *    var = var_temp
            *    var_temp => NULL() // OR IS IT DEALLOCATE as per MOM6 wrappers?
            * endif
             */
          if ( isSgPointerType(varTo->get_type()) ){

            // if the original pointer is to an array, the argument to the allocate statement must be of
            // of the form toVar(SIZE(fromVar, rank_1),SIZE(fromVar, rank_2),...) for each rank_i in rank
            // BEGIN BUILD ALLOCATE BRANCH FOR ORIGINAL VARIABLE IN THE CASE OF REALLOC OR INIT
            vector<SgExpression*> allocateArgList;
            SgArrayType* underlyingArrayType = isSgArrayType(isSgPointerType(varFrom->get_type())->get_base_type());
            if ( underlyingArrayType != NULL ){

              vector<SgExpression*> arrayRefArgList;
              int rank = underlyingArrayType->get_rank();

              for ( int i = 1; i <= rank; ++i ){
                vector<SgExpression*> shapeArgList;
                shapeArgList.push_back(SageBuilder::buildVarRefExp(varFrom, funcBody));
                shapeArgList.push_back(SageBuilder::buildIntVal(i));
                SgFunctionCallExp* shapeCall = SageBuilder::buildFunctionCallExp("SIZE", SageBuilder::buildIntType(), SageBuilder::buildExprListExp(shapeArgList), funcBody);
                arrayRefArgList.push_back(shapeCall);
              }
              allocateArgList.push_back(SageBuilder::buildPntrArrRefExp(SageBuilder::buildVarRefExp(varTo, funcBody), SageBuilder::buildExprListExp(arrayRefArgList)));
            }
            // otherwise, we can just allocate space for the single var
            else{
              allocateArgList.push_back(SageBuilder::buildVarRefExp(varTo, funcBody));
            }

            SgAllocateStatement* allocateStatement = new SgAllocateStatement(Sg_File_Info::generateDefaultFileInfoForTransformationNode());
            allocateStatement->set_parent(funcBody);
            allocateStatement->set_expr_list(SageBuilder::buildExprListExp(allocateArgList));

            vector<SgExpression*> associatedArgList;
            associatedArgList.push_back(SageBuilder::buildVarRefExp(varTo, funcBody));
            SgFunctionCallExp* associatedCall = SageBuilder::buildFunctionCallExp("ASSOCIATED", SageBuilder::buildBoolType(), SageBuilder::buildExprListExp(associatedArgList), funcBody);
            SgNotOp* notAssociatedCall = SageBuilder::buildNotOp(associatedCall);
            SgIfStmt* ifStmtAllocOriginalVar = SageBuilder::buildIfStmt(notAssociatedCall, SageBuilder::buildBasicBlock(allocateStatement), NULL);
            // END BUILD ALLOCATE BRANCH FOR ORIGINAL VARIABLE IN THE CASE OF REALLOC OR INIT

            // nullying original pointer occurs in the case of dealloc or realloc in the wrapped procedure
            SgFunctionCallExp* rhs = SageBuilder::buildFunctionCallExp("NULL", varTo->get_type(), NULL, funcBody);
            SgVarRefExp* lhs = SageBuilder::buildVarRefExp(varTo, funcBody);
            SgPointerAssignOp* pointerAssignOp = new SgPointerAssignOp(Sg_File_Info::generateDefaultFileInfoForTransformationNode(), lhs, rhs, NULL);
            SgExprStatement* nullifyOriginalPointerStatement = SageBuilder::buildExprStatement(pointerAssignOp);

            // if the original pointer is to an array, we also must add a conditional to check for a size mismatch in the case of a 
            // reallocation in the wrapped procedure so we can nullify the original variable pointer
            // BEGIN BUILD NULLIFY BRANCH FOR ORIGINAL VARIABLE IN THE CASE OF REALLOC
            SgIfStmt* ifStmtNullifyOriginalVarDueToRealloc = NULL;
            if (underlyingArrayType != NULL){

              vector<SgExpression*> shapeArgList;
              shapeArgList.push_back(SageBuilder::buildVarRefExp(varFrom, funcBody));
              SgFunctionCallExp* shapeCall1 = SageBuilder::buildFunctionCallExp("SIZE", SageBuilder::buildIntType(), SageBuilder::buildExprListExp(shapeArgList), funcBody);
              shapeArgList[0] = SageBuilder::buildVarRefExp(varTo, funcBody);
              SgFunctionCallExp* shapeCall2 = SageBuilder::buildFunctionCallExp("SIZE", SageBuilder::buildIntType(), SageBuilder::buildExprListExp(shapeArgList), funcBody);
              SgNotEqualOp* notEqualOp = SageBuilder::buildNotEqualOp(shapeCall1, shapeCall2);
              ifStmtNullifyOriginalVarDueToRealloc = SageBuilder::buildIfStmt(notEqualOp, SageBuilder::buildBasicBlock(nullifyOriginalPointerStatement), NULL);

              associatedArgList[0] = SageBuilder::buildVarRefExp(varTo, funcBody);
              associatedCall = SageBuilder::buildFunctionCallExp("ASSOCIATED", SageBuilder::buildBoolType(), SageBuilder::buildExprListExp(associatedArgList), funcBody);
              ifStmtNullifyOriginalVarDueToRealloc = SageBuilder::buildIfStmt(associatedCall, SageBuilder::buildBasicBlock(ifStmtNullifyOriginalVarDueToRealloc), NULL);
            }
            // END BUILD NULLIFY BRANCH FOR ORIGINAL VARIABLE IN THE CASE OF REALLOC

            // BEGIN CHECK FOR ASSOCIATION OF TEMP VARIABLE TO HANDLE CASE OF DEALLOCATION IN WRAPPED PROCEDURE
            associatedArgList[0] = SageBuilder::buildVarRefExp(varFrom, funcBody);
            associatedCall = SageBuilder::buildFunctionCallExp("ASSOCIATED", SageBuilder::buildBoolType(), SageBuilder::buildExprListExp(associatedArgList), funcBody);
            notAssociatedCall = SageBuilder::buildNotOp(associatedCall);

            SgBasicBlock* falseBody = NULL;
            if (ifStmtNullifyOriginalVarDueToRealloc != NULL){
              if ( castViaAssignmentStatement != NULL ){
                falseBody = SageBuilder::buildBasicBlock(ifStmtNullifyOriginalVarDueToRealloc, ifStmtAllocOriginalVar, castViaAssignmentStatement, cleanupTemporaryVariableStatement);
              }
              else{
                falseBody = SageBuilder::buildBasicBlock(ifStmtNullifyOriginalVarDueToRealloc, ifStmtAllocOriginalVar, cleanupTemporaryVariableStatement);
              }
            }
            else{
              if ( castViaAssignmentStatement != NULL ){
                falseBody = SageBuilder::buildBasicBlock(ifStmtAllocOriginalVar, castViaAssignmentStatement, cleanupTemporaryVariableStatement);
              }
              else{
                falseBody = SageBuilder::buildBasicBlock(ifStmtAllocOriginalVar, cleanupTemporaryVariableStatement);
              }
            }
            ifStatement = SageBuilder::buildIfStmt(notAssociatedCall, SageBuilder::buildBasicBlock(nullifyOriginalPointerStatement), falseBody);
          }

          // if varTo is optional, we must include all post-function-call processing inside of the
          // true branch of an "if present" construct
          if ( varTo->get_declaration()->get_declarationModifier().get_typeModifier().isOptional() ){

            // build call to present intrinsic to check for presence of optional arg
            vector<SgExpression*> presentArgList;
            presentArgList.push_back(SageBuilder::buildVarRefExp(varTo, funcBody));
            SgFunctionCallExp* presentCall = SageBuilder::buildFunctionCallExp("PRESENT", SageBuilder::buildBoolType(), SageBuilder::buildExprListExp(presentArgList), funcBody);

            if ( ifStatement != NULL ){
              ifStatement = SageBuilder::buildIfStmt(presentCall, SageBuilder::buildBasicBlock(ifStatement), NULL);
            }
            else if ( castViaAssignmentStatement != NULL ){
              ifStatement = SageBuilder::buildIfStmt(presentCall, SageBuilder::buildBasicBlock(castViaAssignmentStatement, cleanupTemporaryVariableStatement), NULL);
            }
            else{
              ifStatement = SageBuilder::buildIfStmt(presentCall, SageBuilder::buildBasicBlock(cleanupTemporaryVariableStatement), NULL);
            }
          }

          if ( ifStatement != NULL ){
            SageInterface::appendStatement(ifStatement, funcBody);
          }
          else{
            if ( castViaAssignmentStatement != NULL ){
              SageInterface::appendStatement(castViaAssignmentStatement, funcBody);
            }
            if ( cleanupTemporaryVariableStatement != NULL ){
              SageInterface::appendStatement(cleanupTemporaryVariableStatement, funcBody);
            }
          }
        }

        //////////////////////////////////////////////////// END POPULATION OF WRAPPER PROCEDURE BODY

        #if VERBOSE >= 11
        cout << "[END] building wrapper " << wrapperName << endl;
        cout << endl;
        #endif

        // insert procedure wrapper
        SageInterface::appendStatement(wrapperProcHeader, destinationScope);

        // mark as transformation
        for ( const auto& x : NodeQuery::querySubTree(wrapperProcHeader, V_SgLocatedNode) ){
          x->set_containsTransformation(true);
        }

        //////////////////////////////////////////////////// BEGIN SYMBOL EXPORT

        if ( !NO_CONTAINING_MODULE ){

          auto statementList = destinationScope->generateStatementList();
          auto it = statementList.begin();

          // create public export statement
          SgAttributeSpecificationStatement* publicStatement = new SgAttributeSpecificationStatement(Sg_File_Info::generateDefaultFileInfoForTransformationNode());
          publicStatement->set_attribute_kind(SgAttributeSpecificationStatement::e_accessStatement_public);
          publicStatement->set_firstNondefiningDeclaration(publicStatement);
          publicStatement->set_definingDeclaration(publicStatement);
          publicStatement->get_name_list().push_back(wrapperName);

          // place public export statement
          it = statementList.begin();
          auto insertion_ptr = it;
          while ( it != statementList.end() ){
            if ( SgImplicitStatement* implicitStatement = isSgImplicitStatement(*it) ){
              insertion_ptr = it;
            }
            else if ( SgAttributeSpecificationStatement* existingPublicStatement = isSgAttributeSpecificationStatement(*it) ) {
              insertion_ptr = it;
            }
            else if ( SgUseStatement* existingUseStatement = isSgUseStatement(*it) ) {
              while ( ++it != statementList.end() && isSgUseStatement(*(it)) ){
                insertion_ptr = it;
              }
              continue;
            }
            ++it;
          }

          if ( insertion_ptr == statementList.begin() ){
            SageInterface::insertStatementBefore(*insertion_ptr, publicStatement);
          }
          else{
            SageInterface::insertStatementAfter(*insertion_ptr, publicStatement);
          }
        }
        //////////////////////////////////////////////////// END SYMBOL EXPORT
      }

      // replace the function symbol for the original call with the function symbol of the wrapper procedure

      #if VERBOSE >= 11
      cout << "\t\t\t\t querying for the the symbol " << wrapperName << " in the parent symbol tables of " << SageInterface::getEnclosingScope(functionCall) << endl;
      #endif

      SgFunctionSymbol* wrapperSymbol = SageInterface::lookupFunctionSymbolInParentScopes(wrapperName, SageInterface::getEnclosingScope(functionCall));
      if ( wrapperSymbol == NULL){

        #if VERBOSE >= 11
        cout << "\t\t\t\t\t not found; importing from " << isSgClassDefinition(destinationScope)->get_declaration()->get_name() << endl;
        #endif

        // import wrapper name
        SgStatement* use_or_implicit_statement = NULL;
        SgNode* current_node = functionCall;
        vector<SgNode*> node_list;
        while ( use_or_implicit_statement == NULL ){
          SgScopeStatement* enclosing_scope = SageInterface::getEnclosingScope(current_node);
          if ( enclosing_scope == NULL ){
            assert(false);
          }
          node_list = NodeQuery::querySubTree(enclosing_scope, V_SgUseStatement);
          if ( node_list.size() > 0 ){
            use_or_implicit_statement = isSgStatement(node_list.front());
            break;
          }
          node_list = NodeQuery::querySubTree(enclosing_scope, V_SgImplicitStatement);
          if ( node_list.size() > 0 ){
            use_or_implicit_statement = isSgStatement(node_list.front());
            break;
          }
          current_node = enclosing_scope;
        }

        SgUseStatement* newUseStatement = new SgUseStatement(Sg_File_Info::generateDefaultFileInfoForTransformationNode(), isSgClassDefinition(destinationScope)->get_declaration()->get_name(), true);
        SgRenamePair* newRenamePair = new SgRenamePair(Sg_File_Info::generateDefaultFileInfoForTransformationNode(), wrapperName, wrapperName);
        newUseStatement->get_rename_list().push_back(newRenamePair);
        newRenamePair->set_parent(newUseStatement);
        SageInterface::insertStatementBefore(use_or_implicit_statement, newUseStatement);

        wrapperSymbol = SageInterface::lookupFunctionSymbolInParentScopes(wrapperName, destinationScope);
        assert( wrapperSymbol != NULL );
      }

      #if VERBOSE >= 11
      cout << "\t\t\t\t\t replacing original call with call to wrapper" << endl;
      #endif

      SgFunctionRefExp* oldProcRef = isSgFunctionRefExp(functionCall->get_function());
      SgFunctionRefExp* newProcRef = SageBuilder::buildFunctionRefExp(wrapperSymbol);
      SageInterface::replaceExpression(oldProcRef, newProcRef, false);
    }


    void apply_configuration( const SgSourceFile* sourceFile ){

      // get global scope of target file
      SgScopeStatement* topLevelScopePtr = sourceFile->get_globalScope();

      // query for all variables in the scope file
      for ( const auto& x : NodeQuery::querySubTree(topLevelScopePtr, V_SgInitializedName) ){
        if ( SgInitializedName* var = isSgInitializedName(x) ){

          // if they are floating-point type and they are in the
          // target configuration, set their kind accordingly
          if ( is_underlying_real_type(var) ){
            string scopedVarName = get_scoped_name(var);
            if ( targetConfig.count(scopedVarName) > 0 ){
              set_real_kind(var, replacementKinds[targetConfig[scopedVarName]]);
            }
          }
        }
      }
    }


    void preprocess( const SgSourceFile* sourceFile ){

      // get global scope of target file
      SgScopeStatement* topLevelScopePtr = sourceFile->get_globalScope();

      // query for all function calls in the file and preserve their aliasing
      preserve_aliasing(NodeQuery::querySubTree(topLevelScopePtr, V_SgFunctionCallExp));

      // load constant info
      ifstream inFile;
      unordered_set<string> constantSet;
      inFile.open(WORKING_DIR + "prose_workspace/constant_list.txt");
      if ( inFile.is_open() ){
        string const_var;
        while (getline(inFile, const_var)) {
          constantSet.insert(const_var);
        }
        inFile.close();
      }
      else{
        assert(false);
      }

      // gather all nested scopes and process each scope's symbol table for floating-point variables
      auto nestedScopes =  NodeQuery::querySubTree(topLevelScopePtr, V_SgScopeStatement);
      for ( const auto& x : nestedScopes ){
        if ( SgScopeStatement* scopePtr = isSgScopeStatement(x) ){

          // gather all InitializedNames from this scope's symbol table and then preprocess accordingly
          auto symbolSet = scopePtr->get_symbol_table()->get_symbols();
          for ( auto& symbol : symbolSet ){
            if ( SgVariableSymbol* varSymbol = isSgVariableSymbol(symbol) ){
              if ( SgInitializedName* var = isSgInitializedName(varSymbol->get_declaration()) ){
                if ( SgTypeFloat* floatType = is_underlying_real_type(var) ){
                  if ( SgVariableDeclaration* varDecl = isSgVariableDeclaration(var->get_parent()) ){
                    preprocess_var_declaration(var, varDecl, constantSet);
                  }
                }
              }
            }
          }
        }
      }
    }


    void preprocess_var_declaration( SgInitializedName* var, SgVariableDeclaration* targetNode, unordered_set<string>& constantSet ){

      auto variableList = targetNode->get_variables();

      if ( variableList.size() > 1) {
        // get parent of that variable declaration
        if ( SgStatement* variableDeclarationParent = isSgStatement(targetNode->get_parent()) ){

          Rose_STL_Container<SgStatement*> replacementStatements;

          // ...give each variable it's own declaration node
          for ( const auto& k : variableList ){

            if ( SgInitializedName* var = isSgInitializedName(k) ){

              // make sure pointers to reals that are not dummy variables are set to NULL()
              if ( (isSgPointerType(var->get_type())) && (!isSgFunctionDefinition(SageInterface::getEnclosingScope(var))) ){
                var->set_initializer(SageBuilder::buildAssignInitializer(SageBuilder::buildIntVal()));
              }

              SgVariableDeclaration* newDeclaration = new SgVariableDeclaration(Sg_File_Info::generateDefaultFileInfoForTransformationNode());
              newDeclaration->initializeData(Sg_File_Info::generateDefaultFileInfoForTransformationNode(), var);
              newDeclaration->set_parent(variableDeclarationParent);

              // self reference is required for some reason
              newDeclaration->set_definingDeclaration(newDeclaration);

              // copies access modifiers from original variable declaration
              newDeclaration->get_declarationModifier().get_typeModifier().set_modifierVector(targetNode->get_declarationModifier().get_typeModifier().get_modifierVector());

              // changes default C++ "PUBLIC" modifier to fortran default
              newDeclaration->get_declarationModifier().get_accessModifier().setUndefined();

              // set intent(in) when a variable could be propagated from a named constant
              add_constant_intent(var, newDeclaration, constantSet);

              var->set_parent(newDeclaration);
              replacementStatements.push_back(newDeclaration);
            }
          }
          variableDeclarationParent->replace_statement(targetNode,replacementStatements);
        }
      }
      else {
        // only one variable, still performs constant intent(in) check
        add_constant_intent(var, targetNode, constantSet);

        // still make sure pointers to reals that are not dummy variables are set to NULL()
        if ( (isSgPointerType(var->get_type())) && (!isSgFunctionDefinition(SageInterface::getEnclosingScope(var))) ){
          var->set_initializer(SageBuilder::buildAssignInitializer(SageBuilder::buildIntVal()));
        }
      }
    }


    void add_constant_intent(SgInitializedName* var, SgVariableDeclaration* varDecl, unordered_set<string>& constantSet) {
      string name = get_scoped_name(var);
      unordered_set<string>::iterator it = constantSet.find(name);

      if (it != constantSet.end()) {
        varDecl->get_declarationModifier().get_typeModifier().setIntent_in();
      }
    }


    void load_graph(){

      // load G_proc
      ifstream inFile;
      inFile.open(WORKING_DIR + "prose_workspace/G_proc.graph", ios::binary);
      if ( inFile.is_open() ){
        boost::archive::binary_iarchive readArchive(inFile);
        readArchive >> G_proc;
        inFile.close();
      }
      else{
        assert(false);
      }
      inFile.open(WORKING_DIR + "prose_workspace/GP_vertexMap.map", ios::binary);
      if ( inFile.is_open() ){
        boost::archive::binary_iarchive readArchive(inFile);
        readArchive >> GP_vertexMap;
        inFile.close();
      }
      else{
        assert(false);
      }
      inFile.open(WORKING_DIR + "prose_workspace/GP_callMap.map", ios::binary);
      if ( inFile.is_open() ){
        boost::archive::binary_iarchive readArchive(inFile);
        readArchive >> GP_callMap;
        inFile.close();
      }
      else{
        assert(false);
      }
    }


    void infer_kinds_of_literals( vector<SgNode*> floatLiterals, bool standalone ){

      #if VERBOSE >= 11
      cout << "\t inferring kinds of " << floatLiterals.size() << " float literals" << endl;
      #endif

      for (const auto& x : floatLiterals ){
        SgFloatVal* floatLiteral = isSgFloatVal(x);
        assert(floatLiteral != NULL);

        #if VERBOSE >= 11
        cout << "\t\t inferring type of " << floatLiteral->get_valueString() << endl;
        #endif

        // traverse through containing expressions until we've
        // discovered fp non-literals with which we can infer the
        // literal's kind
        vector<int> nonLiteralFloatKinds;
        SgNode* temp = floatLiteral->get_parent();
        if ( standalone ){
          while ( !(isSgExpression(temp) || isSgScopeStatement(temp) ) ){
            temp = temp->get_parent();
          }
        }
        else{
          while ( !(isSgBinaryOp(temp) || isSgScopeStatement(temp) ) ){
            temp = temp->get_parent();

            // if the literal is an argument to a procedure, do not
            // use the types of any non-literals outside of that call
            // to infer its type
            if ( isSgFunctionCallExp(temp) ){
              temp = SageInterface::getGlobalScope(temp); // will be cast to NULL
            }
          }
        }

        SgExpression* containingExpression = isSgExpression(temp);

        while ( !((containingExpression == NULL) || (nonLiteralFloatKinds.size() > 0)) ){

          #if VERBOSE >= 11
          cout << "\t\t\t searching containing expression " << containingExpression->unparseToString() << " for variable references" << endl;
          #endif

          if ( SgAssignInitializer* initAssign = isSgAssignInitializer(containingExpression) ){
            if ( SgInitializedName* var = isSgInitializedName(initAssign->get_parent()) ){
              if ( SgTypeFloat* thisFloatType = isSgTypeFloat(var->get_type()->findBaseType()) ){

                int kind = DEFAULT_KIND;
                if ( SgIntVal* thisKindVal = isSgIntVal(thisFloatType->get_type_kind()) ){
                  kind = thisKindVal->get_value();
                }

                #if VERBOSE >= 11
                cout << "\t\t\t\t discovered " << var->get_name() << " of kind " << kind << endl;
                #endif

                nonLiteralFloatKinds.push_back(kind);
              }
            }
            break;
          }
          else{
            vector<SgVariableSymbol*> varSymbols = SageInterface::getSymbolsUsedInExpression(containingExpression);
            for ( const auto& varSymbol : varSymbols ){
              SgType* varType = varSymbol->get_type()->findBaseType();
              
              // handle imported variables whose declarations are inexplicably variable reference expressions
              if ( SgInitializedName* var = isSgInitializedName(varSymbol->get_declaration()) ){
                if ( boost::algorithm::ends_with(SageInterface::getEnclosingSourceFile(var->get_declaration())->getFileName(), ".rmod") ){
                  var = find_imported_variable_declaration(var, SageInterface::getEnclosingSourceFile(containingExpression));
                  assert ( var != NULL );
                  varType = var->get_type();
                }
              }

              if ( SgTypeFloat* thisFloatType = isSgTypeFloat(varType) ){

                int kind = DEFAULT_KIND;
                if ( SgIntVal* thisKindVal = isSgIntVal(thisFloatType->get_type_kind()) ){
                  kind = thisKindVal->get_value();
                }

                #if VERBOSE >= 11
                cout << "\t\t\t\t discovered" << varSymbol->get_name() << " of kind " << kind << endl;
                #endif

                nonLiteralFloatKinds.push_back(kind);
              }
            }
            temp = containingExpression->get_parent();
            if ( standalone ){
              while ( !(isSgExpression(temp) || isSgScopeStatement(temp) ) ){
                temp = temp->get_parent();
              }
            }
            else{
              while ( !(isSgBinaryOp(temp) || isSgScopeStatement(temp) ) ){
                temp = temp->get_parent();

                // if the literal is an argument to a procedure, do not
                // use the types of any non-literals outside of that call
                // to infer its type
                if ( isSgFunctionCallExp(temp) ){
                  temp = SageInterface::getGlobalScope(temp); // will be cast to NULL
                }
              }
            }
            containingExpression = isSgExpression(temp);
          }
        }

        // if possible, infer and set the kind of this literal
        if ( nonLiteralFloatKinds.size() > 0 ){

          vector<int> histogram (17, 0);
          for ( int kind : nonLiteralFloatKinds ){
            ++histogram[kind];
          }

          string literalKindString = to_string(distance(histogram.begin(), max_element(histogram.begin(), histogram.end())));
          string literalValueString = floatLiteral->get_valueString();
          literalValueString = literalValueString.substr(0, literalValueString.find_last_of("_"));

          // if the literal uses the exponent form with 'd', this locks it in as a double precision float and 
          // adding a precision via the underscore notation will cause a compiler error. Change it to an 'e'
          if ( literalValueString.find_last_of("d") != string::npos ){
            literalValueString = literalValueString.substr(0, literalValueString.find_last_of("d")) + "e" + literalValueString.substr(literalValueString.find_last_of("d")+1, string::npos);
          }
          else if ( literalValueString.find_last_of("D") != string::npos ){
            literalValueString = literalValueString.substr(0, literalValueString.find_last_of("D")) + "E" + literalValueString.substr(literalValueString.find_last_of("D")+1, string::npos);
          }

          // account for float literals that already have a kind specifier
          if ( literalValueString.find_last_of("_") != string::npos ){
            literalValueString = literalValueString.substr(0, literalValueString.find_last_of("_"));
          }

          #if VERBOSE >= 11
          cout << "\t\t\t\t\t changing " << floatLiteral->get_valueString() << " to " << literalValueString + "_" + literalKindString << endl;
          #endif

          floatLiteral->set_valueString(literalValueString + "_" + literalKindString);
        }
      }
      return;
    }


    char kind_to_char(int kind){
      if ( kind == OPTIONAL_ARG_NOT_PROVIDED ){
        return 'x';
      }
      else if ( kind == NON_FLOAT_TYPE ){
        return '0';
      }
      else{
        vector<char> kindToCharMap = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f', 'g'};
        assert(kind <= kindToCharMap.size());
        return kindToCharMap[kind];
      }
    }


    string generate_expected_kind_config_string(SgProcedureHeaderStatement* procHeader){
      string expectedKindConfigString;

      for ( const auto& var : get_dummy_var_decl_order(procHeader) ){
        int kind;
        SgType* baseType = var->get_type()->findBaseType();

        if ( SgTypeFloat* floatType = isSgTypeFloat(baseType) ){
          kind = DEFAULT_KIND;
          if ( SgIntVal* temp = isSgIntVal(floatType->get_type_kind()) ){
            kind = temp->get_value();
          }
        }
        else{
          kind = NON_FLOAT_TYPE;
        }

        expectedKindConfigString.push_back(kind_to_char(kind));
      }

      return expectedKindConfigString;
    }


    string generate_given_kind_config_string(SgFunctionCallExp* functionCall, SgProcedureHeaderStatement* calleeProcHeader){

      string givenKindConfigString;
      vector<SgExpression*> givenArguments = functionCall->get_args()->get_expressions();
      vector<SgInitializedName*> expectedArguments = calleeProcHeader->get_parameterList()->get_args();
      map<string,int> givenVarNamesAndKinds;

      for ( int i = 0; i < givenArguments.size(); ++i ){

        // get the variable name
        string varName;
        if ( SgActualArgumentExpression* keywordAssignment = isSgActualArgumentExpression(givenArguments[i]) ){
            varName = keywordAssignment->get_argument_name().str();
        }
        else{
          varName = expectedArguments[i]->get_name();
        }

        boost::algorithm::to_lower(varName);

        // if it's a float, get its kind and save the name kind pair
        // otherwise, mark it as a NON_FLOAT_TYPE
        if ( SgTypeFloat* givenFloat = isSgTypeFloat(get_expression_type(givenArguments[i])->findBaseType()) ){

          int givenKind = DEFAULT_KIND;
          if ( SgIntVal* temp1 = isSgIntVal(givenFloat->get_type_kind()) ){
            givenKind = temp1->get_value();
          }
          else {
            auto floatLiterals = NodeQuery::querySubTree(givenArguments[i], V_SgFloatVal);
            for ( const auto& x : floatLiterals ){
              if ( SgFloatVal* floatLiteral = isSgFloatVal(x) ){
                string valueString = floatLiteral->get_valueString();
                if ( valueString.find_last_of("_") != string::npos ){
                  string kindString = valueString.substr(valueString.find_last_of("_") + 1, 2);
                  if ( kindString.length() > 0 ){
                    int i = 0;
                    do{
                      if ( !isdigit(kindString[i]) ){
                        break;
                      }
                    }
                    while (++i < kindString.length());

                    if ( i == kindString.length() ){
                      givenKind = stoi(kindString);
                    }
                  }
                }
              }
            }
          }
          givenVarNamesAndKinds[varName] = givenKind;
        }
        else{
          givenVarNamesAndKinds[varName] = NON_FLOAT_TYPE;
        }
      }

      // build the givenKindConfigString in dummy var declaration order
      for ( const auto& originalArg : get_dummy_var_decl_order(calleeProcHeader) ){
        string originalName = originalArg->get_name();
        boost::algorithm::to_lower(originalName);
        if ( givenVarNamesAndKinds.count(originalName) > 0 ){
          givenKindConfigString.push_back(kind_to_char(givenVarNamesAndKinds[originalName]));
        }
        else{
          givenKindConfigString.push_back('x');
        }
      }

      return givenKindConfigString;
    }
}; // end ApplyConfigurationAction


/*
* Input: a pointer to an InitializedName Node or a VariableDeclaration Node
* Output: if underlying base type is a constant, a pointer to that SgInitializedName object. Otherwise, NULL
*/
SgInitializedName* is_constant(SgNode* node){
  assert( (isSgVariableDeclaration(node)) || (isSgInitializedName(node)) );

  SgInitializedName* var;
  SgDeclarationStatement* declStmt;

  // if it is a variableDeclaration node, get the type of the first initializedName in the declaration
  if ( SgVariableDeclaration* varDecl = isSgVariableDeclaration(node) ){
    var = SageInterface::getFirstInitializedName(varDecl);
    declStmt = isSgDeclarationStatement(varDecl);
  }
  // otherwise, just get associated type of given initializedName
  else {
    var = isSgInitializedName(node);
    declStmt = var->get_declaration();
  }

  if ( declStmt && declStmt->get_declarationModifier().get_typeModifier().get_constVolatileModifier().isConst() ) {
    return var;
  }
  else {
    return NULL;
  }
}

/*
* Input: a pointer to an InitializedName Node or a VariableDeclaration Node
* Output: if underlying base type is real, a pointer to that type. Otherwise, NULL
*/
SgTypeFloat* is_underlying_real_type(SgNode* node){

  assert( (isSgVariableDeclaration(node)) || (isSgInitializedName(node)) );

  SgType* variableType;

  // if it is a variableDeclaration node, get the type of the first initializedName in the declaration
  if ( SgVariableDeclaration* varDecl = isSgVariableDeclaration(node) ){
    variableType = SageInterface::getFirstInitializedName(varDecl)->get_type();
  }
  // otherwise, just get associated type of given initializedName
  else if ( SgInitializedName* var = isSgInitializedName(node) ){
    variableType = var->get_type();
  }

  // get base type if it is an array or pointer
  if ( SgArrayType* arrayType = isSgArrayType(variableType) ){
    variableType = arrayType->findBaseType();
  }
  else if ( SgPointerType* pointerType = isSgPointerType(variableType) ){
    variableType = pointerType->findBaseType();
  }

  return isSgTypeFloat(variableType);
}


/*
* Input: a pointer to an InitializedName Node and the replacement type
* Output: 0 for success, non-zero for failure
*
* Sets type of scalars, arrays, and pointers
* with the type specified by replacementType.
*/
void set_real_kind( SgInitializedName* targetNode, SgType* replacementType ){

  SgType* targetNodeType = targetNode->get_type();

  // if the type of the target node does not contain internal types,
  // we can simply set the type of the targetNode
  if ( !(targetNodeType->containsInternalTypes()) ){
    targetNode->set_type(replacementType);
  }

  // otherwise, we wrap the base SgTypeFloat in deep copies of each internal type
  else{
    auto internalTypes = targetNodeType->getInternalTypes();
    assert( isSgTypeFloat(internalTypes.back()) );

    for ( auto it = internalTypes.rbegin() + 1; it != internalTypes.rend(); ++it ){
      SgType* clonedType = isSgType(SageInterface::deepCopyNode(*it));
      if ( SgPointerType* pointerType = isSgPointerType(clonedType) ){
        pointerType->set_base_type(replacementType);
      }
      else if ( SgArrayType* arrayType = isSgArrayType(clonedType) ){
        arrayType->set_base_type(replacementType);
      }
      replacementType = clonedType;
    }
    targetNode->set_type(replacementType);
  }

  // special case: if this is an implicit return value node(i.e. no accompanying variable declaration)
  if ( SgProcedureHeaderStatement* func = isSgProcedureHeaderStatement(targetNode->get_parent()) ){
    if ( SgFunctionType* funcType = isSgFunctionType(func->get_type()) ){
      funcType->set_return_type(targetNode->get_type());
      funcType->set_orig_return_type(targetNode->get_type());
    }
  }

  return;
}


string get_scoped_name( SgNode* node ){

  string scope_str = "";

  if ( isSgInitializedName(node) ){
    SgInitializedName* var = isSgInitializedName(node);

    auto* parent = var->get_parent();

    // special case for function returns
    if ( SgProcedureHeaderStatement* temp = isSgProcedureHeaderStatement(parent) ){
      if ( temp->get_result_name() == var ){
        parent = var->get_scope();
      }
    }

    // special case for imported globals
    else if ( SgProcedureHeaderStatement* procHeader = isSgProcedureHeaderStatement(var->get_declaration()) ){
      for ( const auto& x : NodeQuery::querySubTree(SageInterface::getEnclosingSourceFile(procHeader), V_SgUseStatement) ){
        if ( SgUseStatement* useStatement = isSgUseStatement(x) ){
          for ( const auto& y : useStatement->get_rename_list() ){
            if ( SgRenamePair* renamePair = isSgRenamePair(y) ){
              if ( var->get_name() == renamePair->get_use_name() ){
                return "::" + useStatement->get_name() + "::" + var->get_name();
              }
            }
          }
        }
      }
    }

    do{
      if (SgFunctionDefinition * func = isSgFunctionDefinition(parent)){
        scope_str = func->get_declaration()->get_qualified_name()+"::"+scope_str;
      }
    }
    while (parent = parent->get_parent());

    scope_str = scope_str+(var->get_qualified_name());
  }
  else if ( SgFunctionDefinition* funcDef = isSgFunctionDefinition(node) ){

    SgNode* parent = funcDef->get_parent();

    do{
      if ( SgFunctionDefinition* enclosingFunc = isSgFunctionDefinition(parent) ){
        scope_str = enclosingFunc->get_declaration()->get_qualified_name() + "::" + scope_str;
      }
    }
    while (parent = parent->get_parent());

    scope_str = scope_str+(funcDef->get_declaration()->get_qualified_name());
  }
  else if ( SgClassDefinition* classDef = isSgClassDefinition(node) ){

    SgNode* parent = classDef->get_parent();

    do{
      if ( SgFunctionDefinition* enclosingFunc = isSgFunctionDefinition(parent) ){
        scope_str = enclosingFunc->get_declaration()->get_qualified_name() + "::" + scope_str;
      }
    }
    while (parent = parent->get_parent());

    scope_str = scope_str+(classDef->get_declaration()->get_qualified_name());
  }
  else if ( SgFunctionCallExp* procedureCall = isSgFunctionCallExp(node) ){

    SgNode* parent = SageInterface::getEnclosingScope(node);
    while ( !(isSgFunctionDefinition(parent) || isSgClassDefinition(parent) ) ){
      parent = SageInterface::getEnclosingScope(parent);
    }

    scope_str = get_scoped_name(parent) + procedureCall->unparseToString();
  }
  else {
    assert(false);
  }

  scope_str.erase(remove(scope_str.begin(), scope_str.end(), ' '), scope_str.end());
  boost::algorithm::to_lower(scope_str);

  return scope_str;
}


/*
* Input: a pointer to the root project node, the path to the target file, the scopeString CL arg, and the scopePtr
* Output: the new scopePtr
*
* Sets the scopePtr to the specified tuning scope
*/
SgScopeStatement* get_specified_scope(SgProject* n, string targetScopeString ) {

  SgScopeStatement* scopePtr = isSgScopeStatement(n);   // initialize to NULL

  // split targetScopeString into vector of individual scope names
  vector <string> scopeNodeNames;
  targetScopeString.replace(targetScopeString.find("::"), sizeof("::") - 1, ":");
  boost::split(scopeNodeNames, targetScopeString, boost::is_any_of(":"));

  // get all files included in the project and find the specified scope
  auto fileList = n->get_fileList();
  for ( const auto& x : fileList ){
    if ( SgSourceFile* file = isSgSourceFile(x) ){
      scopePtr = file->get_globalScope();

      // traverse each nested scope to reach specified top-level scope
      for ( const auto& name : scopeNodeNames ){
        if ( name.empty() ){continue;}

        // look for specified nested scope within current scope
        // if the nested scope being searched for isn't found, error
        SgSymbol* nextScope = scopePtr->lookup_symbol(SgName(name));
        if ( nextScope == NULL ){
          scopePtr = NULL;
          break;
        }

        // update scopePtr to current scope
        if ( SgClassSymbol* temp1 = isSgClassSymbol(nextScope) ){
          if ( SgModuleStatement* temp2 = isSgModuleStatement(temp1->get_declaration()->get_definingDeclaration()) ){
            scopePtr = temp2->get_definition();
          }
          else if ( SgDerivedTypeStatement* temp2 = isSgDerivedTypeStatement(temp1->get_declaration()->get_definingDeclaration())){
            scopePtr = temp2->get_definition();
          }
        }
        else if ( SgFunctionSymbol* temp1 = isSgFunctionSymbol(nextScope) ){
          if ( SgProcedureHeaderStatement* temp2 = isSgProcedureHeaderStatement(temp1->get_declaration()->get_definingDeclaration()) ){
            scopePtr = temp2->get_definition();
          }
          else if ( SgProgramHeaderStatement* temp2 = isSgProgramHeaderStatement(temp1->get_declaration()->get_definingDeclaration())){
            scopePtr = temp2->get_definition();
          }
        }
      }

      if ( scopePtr != NULL ){
        return scopePtr;
      }
      else{
        continue;
      }
    }
  }

  // assert(scopePtr != NULL);
  return scopePtr;
}


/*
 *  This function iterates through the argument lists of a procedure call and a procedure definition, performing type comparisons; the optional arguments determine
 *  the function's purpose:
 *        1) if just a mismatchPtr is provided, check_args is supposed to determine whether or not the specified function call corresponds to the given procedure definition.
 *            it should return as soon as a mismatch is encountered
 *        2) if an argKindConfigsPtr is provided, check_args is supposed to check all of the arguments and log the kinds of any floating point variables it encounters to make note of any differences
 *        3) if a accumulatedBoundVariablesPtr is provided, check_args is supposed to check all of the arguments and log the bindings between scoped variable names of floating point variables
 *        4) if an accumulatedConstantListPtr is provided, check_args is supposed to check all of the arguments and log the full scoped names of floating point named constants.
 */
void check_args(SgFunctionCallExp* functionCall,
                SgProcedureHeaderStatement* calleeProcHeader,
                bool* mismatchPtr,
                map<string,pair<int, int>>* argKindConfigsPtr,
                vector<variable_binding>* accumulatedBoundVariablesPtr,
                vector<string>* accumulatedConstantListPtr ) {

  #if VERBOSE >= 11
  cout << "\n[BEGIN] checking args for " << functionCall->unparseToString() << " against expected arguments for " << calleeProcHeader->get_name() << endl;
  #endif

  bool mismatch = false;
  map<string,pair<int,int>> argKindConfigs;
  vector<variable_binding> boundVariables;
  vector<string> constantList;

  vector<SgInitializedName*> expectedArguments = calleeProcHeader->get_parameterList()->get_args();
  vector<SgExpression*> givenArguments = functionCall->get_args()->get_expressions();

  try{
    if ( givenArguments.size() > expectedArguments.size() ){
      throw NOT_THE_SAME_PROCEDURE;
    }

    for ( int i = 0; i < givenArguments.size(); ++i ){

      #if VERBOSE >= 11
      cout << "\t checking " << givenArguments[i]->unparseToString() << endl;
      #endif

      //////////////////////////////////////////////////// BEGIN GET GIVEN ARGUMENT TYPE
      SgType* givenType = get_expression_type(givenArguments[i]);
      //////////////////////////////////////////////////// END GET GIVEN ARGUMENT TYPE

      //////////////////////////////////////////////////// BEGIN GET CORRESPONDING EXPECTED ARGUMENT TYPE

      SgType* expectedType;
      string expectedName;
      SgInitializedName* matchingExpectedArg;

      // if we've encountered a keyword arg, find corresponding expected arg with matching name
      // otherwise, just take the name and type of the expected arg in the corresponding position
      if ( SgActualArgumentExpression* keywordAssignment = isSgActualArgumentExpression(givenArguments[i]) ){

        // get the given argument's name
        string givenName = keywordAssignment->get_argument_name().str();
        boost::algorithm::to_lower(givenName);

        // search for a dummy variable with a matching name
        // abort if none found
        bool found = false;
        vector<SgInitializedName*>::iterator ptrToExpectedArgItr;
        for ( ptrToExpectedArgItr = expectedArguments.begin(); ptrToExpectedArgItr < expectedArguments.end() ; ++ptrToExpectedArgItr){

          matchingExpectedArg = *ptrToExpectedArgItr;
          expectedName = matchingExpectedArg->get_name().str();
          boost::algorithm::to_lower(expectedName);
          expectedType = matchingExpectedArg->get_type();

          if ( expectedName == givenName ){
            found = true;
            break;
          }
        }
        if ( !found ){
          throw NOT_THE_SAME_PROCEDURE;
        }

        expectedArguments.erase(ptrToExpectedArgItr);

        #if VERBOSE >= 11
        cout << "\t\t comparing against corresponding named argument in function definition" << endl;
        #endif
      }
      else{

        #if VERBOSE >= 11
        cout << "\t\t comparing against corresponding positional argument in function definition" << endl;
        #endif

        matchingExpectedArg = expectedArguments.front();
        expectedName = matchingExpectedArg->get_name().str();
        boost::algorithm::to_lower(expectedName);
        expectedType = matchingExpectedArg->get_type();
        expectedArguments.erase(expectedArguments.begin());
      }
      //////////////////////////////////////////////////// END GET CORRESPONDING EXPECTED ARGUMENT TYPE

      #if VERBOSE >= 11
      cout << "\t\t\t Given: " << givenType->class_name() << endl;
      cout << "\t\t\t Expected: " << expectedType->class_name() << endl;
      #endif

      //////////////////////////////////////////////////// BEGIN COARSE-GRAINED TOP-LEVEL CHECK FOR TYPE EQUIVALENCE

      // strip away pointer types
      while ( SgPointerType* pointerType = isSgPointerType(givenType) ){
        givenType = pointerType->get_base_type();
      }
      while ( SgPointerType* pointerType = isSgPointerType(expectedType) ){
        expectedType = pointerType->get_base_type();
      }

      #if VERBOSE >= 11
      cout << "\t\t\t Stripping away pointer types..." << endl;
      cout << "\t\t\t\t Given: " << givenType->class_name() << endl;
      cout << "\t\t\t\t Expected: " << expectedType->class_name() << endl;
      #endif


      // abort if the type names (once pointers are stripped away) don't match
      if ( givenType->class_name() != expectedType->class_name() ){

        // handle bug in ROSE in which mismatching case in derived
        // type names result in the expectedType being an integer; the
        // chances of this being a genuine mismatch seem small...
        if ( (givenType->class_name() == "SgClassType") && (expectedType->class_name() == "SgTypeInt") ){
          continue;
        }
        // handle bug in ROSE in which the result of a binary
        // operation expression between imported constants is typed as
        // an SgFunctionType; again, the chances of this being a
        // genuine mismatch feel small
        else if ( (givenType->class_name() == "SgFunctionType") && (expectedType->class_name() == "SgTypeInt") ){
          continue;
        }
        else if ( (givenType->class_name() == "SgFunctionType") && (expectedType->class_name() == "SgArrayType") ){
          continue;
        }
        else if ( (calleeProcHeader->get_functionModifier().isElemental()) && (givenType->class_name() == "SgArrayType") && ((expectedType->class_name() == "SgTypeFloat") || (expectedType->class_name() == "SgTypeInt")) ){
          // handle ELEMENTAL functions
        }
        else if ( (expectedType->class_name() == "SgArrayType") && ((givenType->class_name() == "SgTypeFloat") || (givenType->class_name() == "SgTypeInt")) && ( NodeQuery::querySubTree(givenArguments[i], V_SgPntrArrRefExp).size() > 0 )){
          // handles this case: https://stackoverflow.com/questions/41630361/passing-scalars-and-array-elements-to-a-procedure-expecting-an-array
        }
        else{
          throw NOT_THE_SAME_PROCEDURE;
        }
      }

      SgArrayType* expectedArrayType = isSgArrayType(expectedType);
      if ( (expectedArrayType != NULL) ){

        #if VERBOSE >= 11
        cout << "\t\t making sure array ranks match" << endl;
        #endif

        SgArrayType* givenArrayType = isSgArrayType(givenType);
        if ( givenArrayType != NULL ){

          int givenRank = get_array_rank(givenArguments[i], givenArrayType);

          #if VERBOSE >= 11
          cout << "\t\t\t Comparing array ranks:" << endl;
          cout << "\t\t\t\t Given: " << givenRank << endl;
          cout << "\t\t\t\t Expected: " << expectedArrayType->get_rank() << endl;
          #endif

          if ( givenRank != expectedArrayType->get_rank() ){
            throw NOT_THE_SAME_PROCEDURE;
          }
        }
      }

      //////////////////////////////////////////////////// END COARSE-GRAINED TOP-LEVEL CHECK FOR TYPE EQUIVALENCE

      //////////////////////////////////////////////////// BEGIN FINER-GRAINED BASE-LEVEL CHECK FOR TYPE EQUIVALENCE

      givenType = givenType->findBaseType();
      expectedType = expectedType->findBaseType();

      #if VERBOSE >= 11
      cout << "\t\t\t\t Getting Base types..." << endl;
      cout << "\t\t\t\t\t Given: " << givenType->class_name() << endl;
      cout << "\t\t\t\t\t Expected: " << expectedType->class_name() << endl;
      #endif

      // abort if the base type names don't match
      if ( givenType->class_name() != expectedType->class_name() ){

        // handle bug in ROSE in which mismatching case in derived
        // type names result in the expectedType being an integer; the
        // chances of this being a genuine mismatch seem small...
        if ( (givenType->class_name() == "SgClassType") && (expectedType->class_name() == "SgTypeInt") ){
          continue;
        }
        // handle bug in ROSE in which the result of a binary
        // operation expression between imported constants is typed as
        // an SgFunctionType; again, the chances of this being a
        // genuine mismatch feel small
        else if ( (givenType->class_name() == "SgFunctionType") && (expectedType->class_name() == "SgTypeInt") ){
          continue;
        }
        else{
          throw NOT_THE_SAME_PROCEDURE;
        }
      }

      //////////////////////////////////////////////////// END FINER-GRAINED BASE-LEVEL CHECK FOR TYPE EQUIVALENCE

      //////////////////////////////////////////////////// BEGIN FINEST-GRAINED KIND-LEVEL CHECK FOR TYPE EQUIVALENCE

      SgTypeFloat* givenFloat = isSgTypeFloat(givenType);
      SgTypeFloat* expectedFloat = isSgTypeFloat(expectedType);

      // if they are not reals, take note and move on to next given argument
      if ( (givenFloat == NULL) && (expectedFloat == NULL) ){
        argKindConfigs[expectedName] = make_pair(NON_FLOAT_TYPE, NON_FLOAT_TYPE);
        continue;
      }

      // extract kinds or use default and record them
      int givenKind = DEFAULT_KIND;
      int expectedKind = DEFAULT_KIND;
      if ( SgIntVal* temp1 = isSgIntVal(givenFloat->get_type_kind()) ){
        givenKind = temp1->get_value();
      }
      if ( SgIntVal* temp2 = isSgIntVal(expectedFloat->get_type_kind()) ){
        expectedKind = temp2->get_value();
      }
      argKindConfigs[expectedName] = make_pair(givenKind, expectedKind);

      #if VERBOSE >= 11
      cout << "\t\t\t\t\t Getting Kinds..." << endl;
      cout << "\t\t\t\t\t\t Given: " << givenKind << endl;
      cout << "\t\t\t\t\t\t Expected: " << expectedKind << endl;
      #endif

      // record if there is a mismatch in kinds
      if ( givenKind != expectedKind ){
        mismatch = true;
      }

      // bind expected variables and all variables that appear in the corresponding given expression
      if ( is_underlying_real_type(matchingExpectedArg) ){
        vector<pair<string,int>> binding;
        binding.push_back(make_pair(get_scoped_name(matchingExpectedArg), expectedKind));
        vector<SgVariableSymbol*> temp = SageInterface::getSymbolsUsedInExpression(givenArguments[i]);
        for ( const auto& symbol : temp ){
          SgInitializedName* var = symbol->get_declaration();
          if ( is_underlying_real_type(var) ){
            binding.push_back(make_pair(get_scoped_name(var), givenKind));
          }
          if ( is_constant(var) != NULL ) {
            constantList.push_back(get_scoped_name(var));
          }
        }
        if ( (binding.size() > 1) && (accumulatedBoundVariablesPtr) ){
          variable_binding temp;
          temp.binding = binding;

          // retrieve profiling info
          double weight = get_edgeWeight_from_profiling_info(functionCall);

          // flag weights representing interprocedural floating point flow by making them negative
          temp.weight = -1 * weight;
          boundVariables.push_back(temp);
        }
      }

      //////////////////////////////////////////////////// END FINEST-GRAINED KIND-LEVEL CHECK FOR TYPE EQUIVALENCE
    } // end iterating through given arguments


    //////////////////////////////////////////////////// BEGIN OPTIONAL-ARG CHECKING
    // iterate through any expected arguments that haven't been processed yet;
    // ensure that they are optional since they were not provided in the function call
    for ( int i = 0; i < expectedArguments.size(); ++i ){

      #if VERBOSE >= 11
      cout << "\t checking if missing argument " << expectedArguments[i]->get_name() << " is optional" << endl;
      #endif

      // abort if a non-optional argument was not present
      if ( !(expectedArguments[i]->get_declaration()->get_declarationModifier().get_typeModifier().isOptional()) ){
        throw NOT_THE_SAME_PROCEDURE;
      }

      // record the config pair
      string expectedName = expectedArguments[i]->get_name().str();
      boost::algorithm::to_lower(expectedName);
      SgType* expectedType = expectedArguments[i]->get_type()->findBaseType();
      if ( SgTypeFloat* expectedFloat = isSgTypeFloat(expectedType) ){

        int expectedKind = DEFAULT_KIND;
        if ( SgIntVal* temp = isSgIntVal(expectedFloat->get_type_kind()) ){
          expectedKind = temp->get_value();
        }

        argKindConfigs[expectedName] = make_pair(OPTIONAL_ARG_NOT_PROVIDED, expectedKind);
      }
      else{
        argKindConfigs[expectedName] = make_pair(OPTIONAL_ARG_NOT_PROVIDED, NON_FLOAT_TYPE);
      }
    }
    //////////////////////////////////////////////////// END OPTIONAL-ARG CHECKING

  }
  catch (int e){
    // catch statement if we've discovered via argument checking that these are not the same procedure

    mismatch = true;
    argKindConfigs.clear();
    boundVariables.clear();
    constantList.clear();
  }

  if ( mismatchPtr ){
    *mismatchPtr = mismatch;
  }
  if ( argKindConfigsPtr ){
    *argKindConfigsPtr = argKindConfigs;
  }
  if ( accumulatedBoundVariablesPtr ){
    (*accumulatedBoundVariablesPtr).insert((*accumulatedBoundVariablesPtr).end(), boundVariables.begin(), boundVariables.end());
  }
  if ( accumulatedConstantListPtr ) {
    (*accumulatedConstantListPtr).insert((*accumulatedConstantListPtr).end(), constantList.begin(), constantList.end());
  }

  #if VERBOSE >= 11
  if ( !argKindConfigs.empty() ){
    cout << "\t constructed the following argKindConfig:" << endl;
    for (const auto& x : argKindConfigs ){
      cout << "\t\t " << x.first << " : kind" << x.second.first << " -> kind" << x.second.second << endl;
    }
  }

  if ( mismatch ){
      cout << "[END] (MISMATCH) checking args for " << functionCall->unparseToString() << " against expected arguments for " << calleeProcHeader->get_name() << endl;
  }
  else{
      cout << "[END] (MATCH) checking args for " << functionCall->unparseToString() << " against expected arguments for " << calleeProcHeader->get_name() << endl;
  }
  cout << endl;
  #endif

  return;
} // end check_args


vector<SgInitializedName*> get_dummy_var_decl_order(SgProcedureHeaderStatement* procHeader){
  vector<SgInitializedName*> dummyVars_declListOrder;
  vector<SgInitializedName*> dummyVars_paramListOrder = procHeader->get_parameterList()->get_args();
  vector<SgStatement*> statementList = procHeader->get_definition()->get_body()->get_statements();
  for (const auto& statement : statementList ){
    if ( SgVariableDeclaration* varDecl = isSgVariableDeclaration(statement) ){
      for (const auto& var : varDecl->get_variables()){
        if ( find(dummyVars_paramListOrder.begin(), dummyVars_paramListOrder.end(), var) != dummyVars_paramListOrder.end() ){
          dummyVars_declListOrder.push_back(var);
        }
      }
    }
  }
  return dummyVars_declListOrder;
}


SgInitializedName* find_imported_variable_declaration( SgInitializedName* var , SgSourceFile* originalFile){

  for ( const auto& x : NodeQuery::querySubTree(originalFile, V_SgUseStatement) ){
    if ( SgUseStatement* useStatement = isSgUseStatement(x) ){
      if ( SgModuleStatement* moduleStatement = useStatement->get_module() ){

        // because of the way that new files are parsed into the AST, the above moduleStatement points to the stale
        // module statement from the rmod file; this is not updated when we apply a configuration, a step that precedes any calls to this
        // function in the context of the ApplyConfiguration plugin action. This leads to faulty reasoning about kind mismatches. 
        // here, we reference the project's file list (which is updated, having removed the reference to the stale rmod file
        // in the parse_source_files_into_AST() function) and find the module with the same name. This appears to work!
        for ( const auto& x : SageInterface::getProject()->get_fileList() ){
          for ( const auto& y : NodeQuery::querySubTree(x, V_SgModuleStatement) ){
            SgModuleStatement* otherModule = isSgModuleStatement(y);
            if ( (moduleStatement->get_name() == otherModule->get_name()) && (moduleStatement != otherModule) ){
              moduleStatement = otherModule;
            }
          }
        }

        // grab the possibly updated definition of this variable from the module in which it is declared
        if ( SgClassDefinition* externalModuleDef = moduleStatement->get_definition() ){
          if ( SgSymbolTable* symbolTable = externalModuleDef->get_symbol_table() ){
            if ( SgVariableSymbol* varSymbol = symbolTable->find_variable(var->get_name()) ){
              return isSgInitializedName(varSymbol->get_declaration());
            }
            else{

              // search through any derived type definitions that may contain the variable that is sought
              for ( const auto& x : NodeQuery::querySubTree(externalModuleDef, V_SgClassDefinition) ){
                SgClassDefinition* derivedTypeDefinition = isSgClassDefinition(x);
                if ( SgSymbolTable* symbolTable = derivedTypeDefinition->get_symbol_table() ){
                  if ( SgVariableSymbol* varSymbol = symbolTable->find_variable(var->get_name()) ){
                    return isSgInitializedName(varSymbol->get_declaration());
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return NULL;
}


void preserve_aliasing( vector<SgNode*> functionCalls ){

  if ( functionCalls.size() == 0 ){
    return;
  }

  // populate map of possible function aliases
  map<string,string> possibleAliases;
  SgSourceFile* sourceFile = SageInterface::getEnclosingSourceFile(functionCalls.front());

  // start by gathering all renamePairs from this module's use statements which are the sites of possible aliasing
  map<string,SgRenamePair*> renamePairsMap;
  for ( const auto& x : NodeQuery::querySubTree(sourceFile, V_SgRenamePair) ){
    if ( SgRenamePair* renamePair = isSgRenamePair(x) ){
      renamePairsMap[renamePair->get_use_name()] = renamePair;
    }
  }

  // for each rename pair
  for ( const auto& x : renamePairsMap ){
    SgRenamePair* localRenamePair = x.second;

    #if VERBOSE >= 11
    cout << "\t Searching symbol tables for symbols corresponding to the imported symbol " << localRenamePair->get_local_name() << endl;
    #endif

    // query symbol table for the function symbol corresponding to the imported symbol
    // and then check to see if there exists a corresponding alias symbol
    SgFunctionSymbol* functionSymbol = NULL;
    SgSymbolTable* symbolTable = NULL;
    SgAliasSymbol* localAliasSymbol = NULL;
    SgScopeStatement* enclosingScope = SageInterface::getEnclosingScope(functionCalls.front());
    while ( enclosingScope != NULL ){
      symbolTable = enclosingScope->get_symbol_table();
      if ( symbolTable != NULL ){
        functionSymbol = symbolTable->find_function(localRenamePair->get_local_name());
        if ( functionSymbol != NULL ){
          localAliasSymbol = symbolTable->find_aliased_symbol(localRenamePair->get_local_name(), functionSymbol);
          if ( localAliasSymbol != NULL){
            break;
          }
        }
      }
      if ( isSgGlobal(enclosingScope) ){
        break;
      }
      enclosingScope = SageInterface::getEnclosingScope(enclosingScope);
    }

    // if the symbol table returned an alias symbol that is also marked as renamed, then aliasing
    // for the imported function symbol occurred in this module; save the alias!
    if ( localAliasSymbol != NULL ) {
      if ( localAliasSymbol->get_isRenamed() ){

        #if VERBOSE >= 11
        cout << "\t\t FOUND LOCAL ALIAS SYMBOL! " << localAliasSymbol->get_new_name() << endl;
        #endif

        possibleAliases[functionSymbol->get_name()] = localAliasSymbol->get_new_name();
      }
      // otherwise, we need to start following a chain of modules
      else{

        SgRenamePair* externalRenamePair = localRenamePair;
        while ( externalRenamePair ){

          // traverse to classDefinition of module from which this symbol was imported
          // perform same check as above for the aliasSymbol
          if ( SgUseStatement* useStatement = isSgUseStatement(externalRenamePair->get_parent()) ){
            if ( SgModuleStatement* moduleStatement = useStatement->get_module() ){
              if ( SgClassDefinition* externalModuleDef = moduleStatement->get_definition() ){
                if ( SgSymbolTable* symbolTable = externalModuleDef->get_symbol_table() ){

                  #if VERBOSE >= 11
                  cout << "\t\t NO ALIAS FOUND; now searching " << moduleStatement->get_name() << endl;
                  #endif

                  SgFunctionSymbol* functionSymbol = symbolTable->find_function(localRenamePair->get_local_name());

                  SgAliasSymbol* externalAliasSymbol = symbolTable->find_aliased_symbol(localRenamePair->get_local_name(), functionSymbol);
                  if ( externalAliasSymbol != NULL ){
                    if ( externalAliasSymbol->get_isRenamed() ){

                      #if VERBOSE >= 11
                      cout << "\t\t\t FOUND EXTERNAL ALIAS SYMBOL! " << externalAliasSymbol->get_name() << endl;
                      #endif

                      // update the local alias symbol with the corret name and mark it as renamed
                      localAliasSymbol->set_new_name(externalAliasSymbol->get_name());
                      localAliasSymbol->set_isRenamed(true);

                      // save the alias
                      possibleAliases[functionSymbol->get_name()] = localAliasSymbol->get_new_name();
                      break;
                    }
                    else{
                      auto renamePairs = NodeQuery::querySubTree(externalModuleDef, V_SgRenamePair);
                      SgRenamePair* oldRenamePair = externalRenamePair;
                      if ( renamePairs.size() > 0 ){
                        for ( const auto& x : renamePairs ){
                          if ( SgRenamePair* temp = isSgRenamePair(x) ){
                            if ( temp->get_local_name() == localRenamePair->get_local_name() ){
                              externalRenamePair = temp;
                              continue;
                            }
                          }
                        }
                        if ( oldRenamePair == externalRenamePair ){

                          #if VERBOSE >= 11
                          cout << "\t\t avoiding loop; " << functionSymbol->get_name() << " determined to not have an alias" << endl;
                          #endif

                          break;
                        }
                      }
                      else{

                      #if VERBOSE >= 11
                      cout << "\t\t no more renamePairs; " << functionSymbol->get_name() << " determined to not have an alias" << endl;
                      #endif

                      break;
                      }
                    }
                  }
                  else{

                    #if VERBOSE >= 11
                    cout << "\t\t no more alias symbols; " << functionSymbol->get_name() << " determined to not have an alias" << endl;
                    #endif

                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  for ( const auto& x : functionCalls ){
    if ( SgFunctionCallExp* functionCall = isSgFunctionCallExp(x) ){
      if ( SgFunctionRefExp* procedureReference = isSgFunctionRefExp(functionCall->get_function()) ){
        if ( SgFunctionSymbol* thisFunctionSymbol = isSgFunctionSymbol(procedureReference->get_symbol()) ){

          string functionName = thisFunctionSymbol->get_name();
          if ( possibleAliases.count(functionName) > 0 ){

            SgFunctionSymbol* functionSymbol = NULL;
            SgSymbolTable* symbolTable = NULL;
            SgAliasSymbol* aliasSymbol = NULL;
            SgScopeStatement* enclosingScope = SageInterface::getEnclosingScope(x);
            while ( enclosingScope != NULL ){
              symbolTable = enclosingScope->get_symbol_table();
              if ( symbolTable != NULL ){
                functionSymbol = symbolTable->find_function(possibleAliases[functionName]);
                if ( functionSymbol != NULL ){
                  aliasSymbol = symbolTable->find_aliased_symbol(possibleAliases[functionName], functionSymbol);
                  if ( aliasSymbol != NULL){
                    break;
                  }
                }
              }
              if ( isSgGlobal(enclosingScope) ){
                break;
              }
              enclosingScope = SageInterface::getEnclosingScope(enclosingScope);
            }

            if ( aliasSymbol != NULL ){
              if ( aliasSymbol->get_isRenamed() ){

                #if VERBOSE >= 11
                cout << "\t aliasing " << thisFunctionSymbol->get_name() << " as " << aliasSymbol->get_name() << endl;
                #endif

                SgRenameSymbol* renameSymbol = new SgRenameSymbol(thisFunctionSymbol->get_declaration(), thisFunctionSymbol, aliasSymbol->get_name());
                renameSymbol->set_parent(SageInterface::getEnclosingClassDefinition(functionCall)->get_symbol_table());
                procedureReference->set_symbol(renameSymbol);
              }
            }
          }
        }
      }
    }
  }
  return;
}


int get_array_rank( SgExpression* expression, SgArrayType* arrayExpressionType ){

  int rank = arrayExpressionType->get_rank();

  auto arrayAccessExpressions = NodeQuery::querySubTree(expression, V_SgPntrArrRefExp);
  if ( arrayAccessExpressions.size() > 0 ){
    for ( const auto& arrayAccess : arrayAccessExpressions ){

      // query for slice operations
      auto subscriptExpressions = NodeQuery::querySubTree(arrayAccess, V_SgSubscriptExpression);
      if ( subscriptExpressions.size() > 0 ){


        // the number of subscript expressions (i.e. those that include one or more colons) is equivalent to the rank of the argument, i.e.
        // array(1,1) has rank 0
        // array(:,1) has rank 1
        // array(1,1:4:2) has rank 1
        rank = subscriptExpressions.size();
      }
    }
  }

  return rank;
}


SgType* get_expression_type( SgExpression* expression ){

  if ( SgActualArgumentExpression* keywordArg = isSgActualArgumentExpression(expression) ){
    expression = keywordArg->get_expression();
  }

  // these operators will always result in a boolean.
  if (isSgOrOp(expression) || isSgAndOp(expression) || isSgNotOp(expression) ||
      isSgEqualityOp(expression) || isSgGreaterOrEqualOp(expression) || isSgGreaterThanOp(expression) ||
      isSgNotEqualOp(expression) || isSgLessOrEqualOp(expression) || isSgLessThanOp(expression)
  ) {
    return SgTypeBool::createType();
  }
  
  if ( SgBinaryOp* binaryOp = isSgBinaryOp(expression) ){

    // if this is a field-access operand for a derived type (%), ensure
    // the resulting type is taken to be the type of the rightmost
    // expression rather than SgClassType for the derived type
    if ( SgDotExp* dotExp = isSgDotExp(binaryOp) ){
      return get_expression_type(dotExp->get_rhs_operand());
    }
    else if ( !isSgPntrArrRefExp(expression) ){

      SgType* lhsType = get_expression_type(binaryOp->get_lhs_operand());
      SgType* rhsType = get_expression_type(binaryOp->get_rhs_operand());

      if ( SgTypeFloat* lhsFloat = isSgTypeFloat(lhsType) ){
        if ( SgTypeInt* rhsInt = isSgTypeInt(rhsType) ){
          return lhsType;
        }
        else if ( SgTypeFloat* rhsFloat = isSgTypeFloat(rhsType) ){

          if ( isSgValueExp(binaryOp->get_lhs_operand()) ){
            return rhsType;
          }
          else if ( isSgValueExp(binaryOp->get_rhs_operand()) ){
            return lhsType;
          }
          else{
            int lhsKind = DEFAULT_KIND;
            int rhsKind = DEFAULT_KIND;

            if ( SgIntVal* temp = isSgIntVal(lhsFloat->get_type_kind()) ){
              lhsKind = temp->get_value();
            }
            if ( SgIntVal* temp = isSgIntVal(rhsFloat->get_type_kind()) ){
              rhsKind = temp->get_value();
            }

            if ( lhsKind > rhsKind ){
              return lhsType;
            }
            else{
              return rhsType;
            }
          }
        }
      }
      else if ( SgTypeFloat* rhsFloat = isSgTypeFloat(rhsType) ){
        if ( SgTypeInt* lhsInt = isSgTypeInt(lhsType) ){
          return rhsType;
        }
      }
      else{
        return rhsType;
      }
    }
  }

  SgType* expressionType = expression->get_type();

  // if there is an array slice operation, make sure that the resulting type
  // is noted to be arrayType; the default get_type operation does not do this.
  auto arrayAccessExpressions = NodeQuery::querySubTree(expression, V_SgPntrArrRefExp);
  if ( arrayAccessExpressions.size() > 0 ){
    for ( const auto& arrayAccess : arrayAccessExpressions ){

      // query for slice operations
      auto subscriptExpressions = NodeQuery::querySubTree(arrayAccess, V_SgSubscriptExpression);
      if ( subscriptExpressions.size() > 0 ){
        SgExpression* expressionContainingArray = isSgPntrArrRefExp(arrayAccess)->get_lhs_operand();
        assert(expressionContainingArray != NULL);
        SgVariableSymbol* firstSymbol = SageInterface::getSymbolsUsedInExpression(expressionContainingArray).front(); // this makes an assumption that multiple symbols won't be involved here
        expressionType = firstSymbol->get_type();

        return expressionType;
      }
    }
  }

  // reason about kinds of imported float vars
  auto variableSymbols = SageInterface::getSymbolsUsedInExpression(expression);
  if ( (variableSymbols.size() == 1) || (arrayAccessExpressions.size() > 0) ){
    SgVariableSymbol* variableSymbol = variableSymbols[0];
    if ( SgInitializedName* var = isSgInitializedName(variableSymbol->get_declaration()) ){
      if (is_underlying_real_type(var)){
        if ( boost::algorithm::ends_with(SageInterface::getEnclosingSourceFile(var->get_declaration())->getFileName(), ".rmod") ||
            ( isSgProcedureHeaderStatement(var->get_declaration()) ) ){
          var = find_imported_variable_declaration(var, SageInterface::getEnclosingSourceFile(expression));
          assert ( var != NULL );
          expressionType = var->get_type();
          if (arrayAccessExpressions.size() > 0){
            expressionType = expressionType->findBaseType();
          }
          return expressionType;
        }
      }
    }
  }

  // handle array access expressions
  if (arrayAccessExpressions.size() > 0){
    return expressionType->findBaseType(); 
  }

  // if the expression is a literal array, reason about the type accordingly
  if ( SgTypeDefault* t = isSgTypeDefault(expression->get_type()) ){
    string argText = expression->unparseToString();
    if ( boost::algorithm::starts_with(argText, "(/") ){

      // split text of the literal array
      vector<string> elementsText;
      boost::split(elementsText, argText, boost::is_any_of(","));
      if ( elementsText.size() > 0 ){

        // assuming array literals are of rank 1
        vector<SgExpression*> temp = {SageBuilder::buildIntVal(elementsText.size())};
        SgExprListExp* dimInfo = SageBuilder::buildExprListExp(temp);

        SgType* arrayElementType;

        // check for variable references
        string firstElementText = elementsText.front().substr(2,elementsText.front().size()-1);
        SgVariableSymbol* varSymbol = isSgVariableSymbol(SageInterface::lookupVariableSymbolInParentScopes(firstElementText, SageInterface::getEnclosingScope(expression)));
        if ( varSymbol ){
          arrayElementType = varSymbol->get_declaration()->get_type();
        }

        // check for integer literals
        else if ( (elementsText.front()).find('.') == string::npos ){
          arrayElementType = SageBuilder::buildIntType();
        }

        // otherwise, they are float literals. Check their kind
        else{

          int kind = DEFAULT_KIND;

          // if the literal is tagged with a kind, check the text
          string valueString = *(elementsText.begin());
          if ( valueString.find_last_of("_") != string::npos ){
            string kindString = valueString.substr(valueString.find_last_of("_") + 1, 2);
            if ( kindString.length() > 0 ){
              int i = 0;
              do{
                if ( !isdigit(kindString[i]) ){
                  break;
                }
              }
              while (++i < kindString.length());

              if ( i == kindString.length() ){
                kind = stoi(kindString);
              }
            }
          }
          arrayElementType = SgTypeFloat::createType(SageBuilder::buildIntVal(kind));
        }

        expressionType = SageBuilder::buildArrayType(arrayElementType, dimInfo);
        return expressionType;
      }
    }
  }

  // special case for the "real" intrinsic function
  if ( SgFunctionCallExp* procCall = isSgFunctionCallExp(expression) ){
    std::string procCallName = procCall->getAssociatedFunctionSymbol()->get_name();
    boost::algorithm::to_lower(procCallName);
    if ( procCallName == "real" ){
      int kind = DEFAULT_KIND;
      auto argExpressions = procCall->get_args()->get_expressions();
      if ( argExpressions.size() > 1 ){
        if ( SgIntVal* intVal = isSgIntVal(argExpressions.back()) ){
          kind = intVal->get_value();
        }
      }
      if ( SgArrayType* arrType = isSgArrayType(get_expression_type(argExpressions.front())) ){
        expressionType = SageBuilder::buildArrayType(SgTypeFloat::createType(SageBuilder::buildIntVal(kind)), arrType->get_dim_info());
      }
      else{
        expressionType = SgTypeFloat::createType(SageBuilder::buildIntVal(kind));
      }
      return expressionType;
    }
    else if ( procCallName == "dble" ){
      int kind = DEFAULT_KIND;
      auto argExpressions = procCall->get_args()->get_expressions();
      if ( SgArrayType* arrType = isSgArrayType(get_expression_type(argExpressions.front())) ){
        expressionType = SageBuilder::buildArrayType(SgTypeFloat::createType(SageBuilder::buildIntVal(kind)), arrType->get_dim_info());
      }
      else{
        expressionType = SgTypeFloat::createType(SageBuilder::buildIntVal(kind));
      }
      return expressionType;
    }
    else if ( procCallName == "sum" ){
      int kind = DEFAULT_KIND;
      auto argExpressions = procCall->get_args()->get_expressions();
      if ( argExpressions.size() > 1 ){
        if ( SgIntVal* intVal = isSgIntVal(argExpressions.front()) ){
          kind = intVal->get_value();
        }
      }
      expressionType = SgTypeFloat::createType(SageBuilder::buildIntVal(kind));
      return expressionType;
    }
    else if ( procCallName== "associated" ){
      return SgTypeBool::createType();
    }
    else if ( procCallName == "size" ){
      return SageBuilder::buildIntType();
    }
    else if ( procCallName == "trim" ){
      auto argExpressions = procCall->get_args()->get_expressions();
      if ( argExpressions.size() == 1) {
        return SageBuilder::buildStringType(SageBuilder::buildIntVal(1));
      }
      else {
        assert(false);
      }
    }
    else if ( procCallName == "sqrt" ||
              procCallName == "abs" ) {
      auto argExpressions = procCall->get_args()->get_expressions();
      if ( argExpressions.size() == 1 ){
        expressionType = get_expression_type(argExpressions[0]);
        return expressionType;
      }
    }
  }
  return expressionType;
}

void intraprocedural_variable_reasoning( vector<SgNode*> binaryOps, vector<variable_binding> & intraboundVariables ){

  if ( binaryOps.size() == 0 ){
    return;
  }

  // this code operates on the assumption that the
  // NodeQuery::querySubTree method queries in a DFS manner and
  // returns results in that order.
  // SgNode* latestProcessed = binaryOps.front();
  for ( const auto& x : binaryOps ){

    // if ( SageInterface::isAncestor(latestProcessed, x) ){
    //   continue;
    // }
    // else{
    //   latestProcessed = x;
    // }

    if ( SgBinaryOp* binaryOp = isSgBinaryOp(x) ){

      bool pointerAssignment = false;
      SgExpression* lhs = binaryOp->get_lhs_operand();
      SgExpression* rhs = binaryOp->get_rhs_operand();

      auto lhsVarRefs = NodeQuery::querySubTree(lhs, V_SgVarRefExp);
      vector<pair<string,int>> lhsVars;
      for ( const auto& x : lhsVarRefs ){
        if ( SgVarRefExp* varRef = isSgVarRefExp(x) ){
          if ( SgTypeFloat* floatType = isSgTypeFloat(varRef->get_type()->findBaseType() )){
            SgInitializedName* var = varRef->get_symbol()->get_declaration();
            string scopedVarName = get_scoped_name(var);

            int kind = DEFAULT_KIND;
            if ( SgIntVal* intVal = isSgIntVal(floatType->get_type_kind()) ){
              kind = intVal->get_value();
            }

            lhsVars.push_back(make_pair(scopedVarName, kind));

            if ( isSgPointerType(varRef->get_symbol()->get_type()) ){
              pointerAssignment = true;
            }
          }
        }
      }

      auto rhsVarRefs = NodeQuery::querySubTree(rhs, V_SgVarRefExp);
      vector<pair<string,int>> rhsVars;
      for ( const auto& x : rhsVarRefs ){
        if ( SgVarRefExp* varRef = isSgVarRefExp(x) ){
          if ( SgTypeFloat* floatType = isSgTypeFloat(varRef->get_type()->findBaseType() )){
            SgInitializedName* var = varRef->get_symbol()->get_declaration();
            string scopedVarName = get_scoped_name(var);

            int kind = DEFAULT_KIND;
            if ( SgIntVal* intVal = isSgIntVal(floatType->get_type_kind()) ){
              kind = intVal->get_value();
            }

            rhsVars.push_back(make_pair(scopedVarName, kind));
          }
        }
      }

      if ( lhsVars.size() + rhsVars.size() > 1 ){

        // retrieve profiling info if it is available
        double weight = get_edgeWeight_from_profiling_info(binaryOp);    
        
        // if ( pointerAssignment ){
        //   // since pointer assignments do not cast automatically,
        //   // instead resulting in compilation errors, we need to weigh such bindings heavily when performing clustering.
        //   // flag it by making it negative
        //   weight = -1 * weight;
        // }

        // add a single binding between the each lhsVar (there
        // should be only 1?) and each rhsVar.
        // TODO: refine this; how can we better weigh
        // intraprocedural dataflow between floating-point
        // variables?
        for ( const auto& l : lhsVars ){
          for ( const auto& r : rhsVars ){
            variable_binding temp;
            temp.binding = {l,r};
            temp.weight = weight;
            intraboundVariables.push_back(temp);
          }
        }
      }
    }
  }
}

void intraprocedural_variable_reasoning2( vector<SgNode*> binaryOps, vector<variable_binding> & intraboundVariables ){

  if ( binaryOps.size() == 0 ){
    return;
  }

  // this code operates on the assumption that the
  // NodeQuery::querySubTree method queries in a DFS manner and
  // returns results in that order.
  for ( const auto& x : binaryOps ){
    if ( SgBinaryOp* binaryOp = isSgBinaryOp(x) ){

      bool pointerAssignment = false;
      SgExpression* lhs = binaryOp->get_lhs_operand();
      SgExpression* rhs = binaryOp->get_rhs_operand();

      auto childBinaryOps = NodeQuery::querySubTree(lhs, V_SgBinaryOp);
      if ( childBinaryOps.size() > 0 ){
        continue;
      }
      childBinaryOps = NodeQuery::querySubTree(rhs, V_SgBinaryOp);
      if ( childBinaryOps.size() > 0 ){
        continue;
      }

      auto lhsVarRefs = NodeQuery::querySubTree(lhs, V_SgVarRefExp);
      vector<pair<string,int>> lhsVars;
      for ( const auto& x : lhsVarRefs ){
        if ( SgVarRefExp* varRef = isSgVarRefExp(x) ){
          if ( SgTypeFloat* floatType = isSgTypeFloat(varRef->get_type()->findBaseType() )){
            SgInitializedName* var = varRef->get_symbol()->get_declaration();
            string scopedVarName = get_scoped_name(var);

            int kind = DEFAULT_KIND;
            if ( SgIntVal* intVal = isSgIntVal(floatType->get_type_kind()) ){
              kind = intVal->get_value();
            }

            lhsVars.push_back(make_pair(scopedVarName, kind));

            if ( isSgPointerType(varRef->get_symbol()->get_type()) ){
              pointerAssignment = true;
            }
          }
        }
      }

      auto rhsVarRefs = NodeQuery::querySubTree(rhs, V_SgVarRefExp);
      vector<pair<string,int>> rhsVars;
      for ( const auto& x : rhsVarRefs ){
        if ( SgVarRefExp* varRef = isSgVarRefExp(x) ){
          if ( SgTypeFloat* floatType = isSgTypeFloat(varRef->get_type()->findBaseType() )){
            SgInitializedName* var = varRef->get_symbol()->get_declaration();
            string scopedVarName = get_scoped_name(var);

            int kind = DEFAULT_KIND;
            if ( SgIntVal* intVal = isSgIntVal(floatType->get_type_kind()) ){
              kind = intVal->get_value();
            }

            rhsVars.push_back(make_pair(scopedVarName, kind));
          }
        }
      }

      if ( lhsVars.size() + rhsVars.size() > 1 ){

        // retrieve profiling info if it is available
        double weight = get_edgeWeight_from_profiling_info(binaryOp);
        // if ( pointerAssignment ){
        //   // since pointer assignments do not cast automatically,
        //   // instead resulting in compilation errors, we need to weigh such bindings heavily when performing clustering.
        //   // flag it by making it negative
        //   weight = -1 * weight;
        // }

        // add a single binding between the each lhsVar (there
        // should be only 1?) and each rhsVar.
        // TODO: refine this; how can we better weigh
        // intraprocedural dataflow between floating-point
        // variables?
        for ( const auto& l : lhsVars ){
          for ( const auto& r : rhsVars ){
            variable_binding temp;
            temp.binding = {l,r};
            temp.weight = weight;
            intraboundVariables.push_back(temp);
          }
        }
      }
    }
  }
}


/*
 * For all given file paths, parse the source into the AST and set the
 * unparse location while removing any lingering references to
 * placeholder rmod files.
 */
unordered_set<SgSourceFile*> parse_source_files_into_AST( SgProject* n, vector<string> filePaths, string unparseLocation ){

  unordered_set<SgSourceFile*> sourceFilesToBeUnparsed;
  for ( string filePath : filePaths ){

    // construct different flavors of filename
    string fileNameWithExtLower = boost::algorithm::to_lower_copy(filePath.substr(filePath.find_last_of("/") + 1, string::npos));
    string rmodFileNameLower = boost::algorithm::to_lower_copy(fileNameWithExtLower.substr(0, fileNameWithExtLower.find_last_of(".")) + ".rmod");
    string fullUnparseName = unparseLocation + "/" + filePath.substr(filePath.find_last_of("/") + 1, string::npos);

    vector<SgFile*> fileList = n->get_fileList();

    // search for corresponding SgFile in the project
    bool already_parsed = false;
    auto it = fileList.begin();
    while ( it != fileList.end() ){

      // if it's an rmod file, clear the global scope's symbol table of symbols before parsing in full source
      if (rmodFileNameLower == boost::algorithm::to_lower_copy(isSgSourceFile(*it)->get_sourceFileNameWithoutPath()) ){        
        while ( isSgSourceFile(*it)->get_globalScope()->symbol_table_size() > 0 ){
          isSgSourceFile(*it)->get_globalScope()->remove_symbol(isSgSourceFile(*it)->get_globalScope()->first_any_symbol());
        }
        break;
      }
      else if ( fileNameWithExtLower == boost::algorithm::to_lower_copy(isSgSourceFile(*it)->get_sourceFileNameWithoutPath())){
        isSgSourceFile(*it)->set_unparse_output_filename(fullUnparseName);
        sourceFilesToBeUnparsed.insert(isSgSourceFile(*it));
        already_parsed = true;
        break;
      }

      ++it;
    }

    // parse full source code and add it to project if necessary
    if ( !already_parsed ){

      SgSourceFile* newFile = isSgSourceFile(SageBuilder::buildFile( filePath, fullUnparseName, n));
      
      // either add the brand new file or replace the existing rmod file
      if ( it == fileList.end() ){
        fileList.push_back(newFile);
      }
      else{
        replace(fileList.begin(), fileList.end(), *it, isSgFile(newFile) );
      }
      n->set_fileList(fileList);
      sourceFilesToBeUnparsed.insert(newFile);
    }
  }
  return sourceFilesToBeUnparsed;
}


void init_edgeWeight_info( SgSourceFile* sourceFile ) {
  string filePath = sourceFile->get_sourceFileNameWithPath();

  boost::to_upper(filePath);
  boost::replace_all(filePath,"/","_");
  boost::replace_all(filePath,".","_");

  ifstream inFile;
  inFile.open(WORKING_DIR + "prose_workspace/__profiling/code_coverage/" + filePath + ".bcov");
  if ( inFile.is_open() ){

    string inputBuffer;
    vector<string> splitBuffer;

    while ( getline(inFile, inputBuffer) ){

      boost::split(splitBuffer, inputBuffer, boost::is_any_of(":") );

      int lineNo = stoi(boost::algorithm::trim_copy(splitBuffer[0]));
      double edgeWeight = stod(boost::algorithm::trim_copy(splitBuffer[1]));
      string codeText = splitBuffer[2];
      for ( int i = min(3, int(splitBuffer.size()));  i < splitBuffer.size(); ++i ){
        codeText = codeText + ":" + splitBuffer[i];
      }

      if (lineNo >= EXECUTION_COUNTS.size()) {
        EXECUTION_COUNTS.resize(lineNo + 1);
        ORIGINAL_CODE_TEXT.resize(lineNo + 1);
      }

      EXECUTION_COUNTS[lineNo] = edgeWeight;
      ORIGINAL_CODE_TEXT[lineNo] = boost::algorithm::erase_all_copy(codeText, " ");
    }
    inFile.close();
  }
}

double get_edgeWeight_from_profiling_info( SgNode* targetNode ){

  double weight = -1.0;

  if ( EXECUTION_COUNTS.size() > 0 ){
    SgNode* currentNode = targetNode;
    string modCode = "";

    // declare vars to control statement search termination
    int maxDepth = 50;
    int depth = 0;

    // if this is an array declaration, get the node of the containing scope
    if ( isSgDeclarationStatement(targetNode) ){
      currentNode = SageInterface::getEnclosingFunctionDeclaration(targetNode);
      if ( !isSgFunctionDeclaration(currentNode) ){
        currentNode = SageInterface::getEnclosingClassDeclaration(targetNode);
      }
      modCode = boost::algorithm::trim_copy(targetNode->unparseToString());
    }

    // otherwise, get the node of the statement by making sure it's parent is a basic block
    else{
      while ((currentNode->get_parent() != NULL) && (depth < maxDepth)){
        depth++;
        currentNode = currentNode->get_parent();
        if ( isSgBasicBlock(currentNode->get_parent()) ){
          assert(isSgStatement(currentNode));
          break;
        }
        // handles a case revealed by micro_mg_utils.F90; module variables declared with binary ops of constants in the initialization
        // the compiler should optimize this out
        else if ( isSgVariableDeclaration(currentNode) ){
          return 0.0;
        } 
      } 
      if ( isSgBasicBlock(currentNode->get_parent()) ){
        modCode = boost::algorithm::trim_copy(currentNode->unparseToString());
      }
    }

    assert ( modCode != "" );
    boost::algorithm::to_lower(modCode);
    boost::algorithm::erase_all(modCode, " ");

    int startingLineNo;
    if ( (isSgClassDeclaration(currentNode)) && !(isSgDerivedTypeStatement(currentNode)) ){
      startingLineNo = isSgLocatedNode(currentNode)->get_file_info()->get_line();
    }
    else{
      // look for attached comment above the function call with the line number information from the profiling run
      AttachedPreprocessingInfoType* comments = isSgLocatedNode(currentNode)->getAttachedPreprocessingInfo();
      depth = 0;
      while ((!comments) && (depth < maxDepth)){
        depth++;
        currentNode = SageInterface::getPreviousStatement(isSgStatement(currentNode));
        comments = isSgLocatedNode(currentNode)->getAttachedPreprocessingInfo();
      }
      assert(comments != NULL);
      assert(boost::algorithm::starts_with((*prev(comments->end()))->getString(), "!PROSE_"));
      string comment = (*prev(comments->end()))->getString();
      startingLineNo = atoi(comment.c_str() + 7); // len("!PROSE_")
    }

    // starting from the line number in the attached comment, search for the most similar source code line within the window of the lookahead
    int targetLineNo;
    int lookAhead = 100;
    
    double maxSimilarity = -1;
    double similarity;
    for ( int i = startingLineNo; i < min(startingLineNo + lookAhead,int(ORIGINAL_CODE_TEXT.size())); ++i ){

      if ( boost::algorithm::starts_with(boost::algorithm::trim_copy(ORIGINAL_CODE_TEXT[i]), "!") ){
        ++lookAhead;
        continue;
      }

      similarity = calc_similarity( modCode, ORIGINAL_CODE_TEXT[i] );
      if ( similarity == 1 ){
        targetLineNo = i;
        break;
      }
      else if ( similarity > maxSimilarity ){
        maxSimilarity = similarity;
        targetLineNo = i;
      }
    }

    while (weight < 0.0) {
      if (targetLineNo >= EXECUTION_COUNTS.size()){
        break;
      }
      weight = EXECUTION_COUNTS[targetLineNo];

      if (weight >= 0.0){
        weight = abs(weight);
        cout << "(line no " << targetLineNo << ") " << isSgLocatedNode(targetNode)->unparseToString() << ", value = " << EXECUTION_COUNTS[targetLineNo] <<endl;
        cout << "\tExpected code: " << modCode << endl;
        cout << "\t Matched code: " << ORIGINAL_CODE_TEXT[targetLineNo] << endl;
        break;
      }
      else{
        targetLineNo++;
      }
    }
  }

  if (weight < 0.0){
    //cout << "**COULDN'T FIND INFO FOR " << isSgLocatedNode(targetNode)->unparseToString() << endl;
    weight = 1.0;
  }

  return weight;
}


//---------------------------------------------------------------------------
// Similarity based on Sørensen–Dice index
double calc_similarity( string s1, string s2 )
{

    // Check banal cases
    if( s1.empty() || s2.empty() )
       {// Empty string is never similar to another
        return 0.0;
       }
    else if( s1==s2 )
       {// Perfectly equal
        return 1.0;
       }
    else if( s1.length()==1 || s2.length()==1 )
       {// Single (not equal) characters have zero similarity
        return 0.0;
       }

    /////////////////////////////////////////////////////////////////////////
    // Represents a pair of adjacent characters
    class charpair_t final
    {
     public:
         charpair_t(const char a, const char b) noexcept : c1(a), c2(b) {}
         [[nodiscard]] bool operator==(const charpair_t& other) const noexcept { return c1==other.c1 && c2==other.c2; }
     private:
        char c1, c2;
    };

    /////////////////////////////////////////////////////////////////////////
    // Collects and access a sequence of adjacent characters (skipping spaces)
    class charpairs_t final
    {
     public:
         charpairs_t(string s)
            {
             assert( !s.empty() );
             const std::size_t i_last = s.size()-1;
             std::size_t i = 0;
             chpairs.reserve(i_last);
             while( i<i_last )
               {
                // Accepting also single-character words (the second is a space)
                //if( !std::isspace(s[i]) ) chpairs.emplace_back( std::tolower(s[i]), std::tolower(s[i+1]) );
                // Skipping single-character words (as in the original article)
                if( std::isspace(s[i]) ) ; // Skip
                else if( std::isspace(s[i+1]) ) ++i; // Skip also next
                else chpairs.emplace_back( std::tolower(s[i]), std::tolower(s[i+1]) );
                ++i;
               }
            }

         [[nodiscard]] auto size() const noexcept { return chpairs.size(); }
         [[nodiscard]] auto cbegin() const noexcept { return chpairs.cbegin(); }
         [[nodiscard]] auto cend() const noexcept { return chpairs.cend(); }
         auto erase(std::vector<charpair_t>::const_iterator i) { return chpairs.erase(i); }

     private:
        std::vector<charpair_t> chpairs;
    };

    charpairs_t chpairs1{s1},
                chpairs2{s2};
    const double orig_avg_bigrams_count = 0.5 * static_cast<double>(chpairs1.size() + chpairs2.size());
    std::size_t matching_bigrams_count = 0;
    for( auto ib1=chpairs1.cbegin(); ib1!=chpairs1.cend(); ++ib1 )
       {
        for( auto ib2=chpairs2.cbegin(); ib2!=chpairs2.cend(); )
           {
            if( *ib1==*ib2 )
               {
                ++matching_bigrams_count;
                ib2 = chpairs2.erase(ib2); // Avoid to match the same character pair multiple times
                break;
               }
            else ++ib2;
           }
       }
    return static_cast<double>(matching_bigrams_count) / orig_avg_bigrams_count;
}


static Rose::PluginRegistry::Add<GenerateGraph>  uniquePluginName0("prose-generate-graph", "build FP dataflow graph");
static Rose::PluginRegistry::Add<LinkGraph> uniquePluginName1("prose-link-graph", "link FP dataflow graph");
static Rose::PluginRegistry::Add<ApplyConfiguration>  uniquePluginName4("prose-apply-configuration", "applies kind configuration specified by command line arguments");

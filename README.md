# PaMILO
## Parametric Mixed Integer Linear Optimization

Implementation of the dual Benson algorithm to computes the breakpoints (and facets) of a parametric MILP or MIQP.

These problems are also known as Multi-Objective Mixed Integer Linear Programs (MOMILPs) and
Multi-Objective Mixed Integer Quadratic Programs (MOMIQPs).

The only assumption on the input instances is that an ideal point exists.

The implementation is based on the paper:

Bökler, Fritz and Mutzel, Petra: Output-Sensitive Algorithms for Enumerating the Extreme Nondominated Points of Multiobjective Combinatorial Optimization Problems, Algorithms-ESA, LNCS 9294, Springer, DOI 10.1007/978-3-662-48350-3_25, 2015.

For support please contact Mirko H. Wagner (mirwagner@uos.de), Levin Nemesch (lnemesch@uos.de),
 or Fritz Bökler (fboekler@uos.de).

# Installation
## UNIX/Linux
```
unzip pamilo.zip
cd pamilo
cmake .
make
```

## Windows
As none of the authors of this program are experts in Windows whatsoever,
please feel free to let us know if the configuration and compilation process
can be made more eﬀicient or the descriptions more precise.


 - Unzip pamilo.zip and go into the directory it is unzipped into.
 - Run `<path to cmake>\cmake.exe .` . If you have a standard installation of CPLEX this should suﬀice, but sometimes it does not. Refer to
the Gurobi section below.

 - Start Visual Studio and build pamilo_cli in Release mode.
 - Go to `<path to pamilo>\Release` and run `pamilo_cli.exe` (for detailed
usage see below).

## CMake
Version 3.18 or higher is needed.

## Solver
PaMILO needs either CPLEX or Gurobi for LP/MILP/MIQP solving.
Capabilities of PaMILO depend on the underlying solver.
At the moment (11.2022), only Gurobi offers the option of nonconvex quadratic problem solving for multiobjective optimization.

You can switch between the solvers in cmake, default is Gurobi.

### Gurobi
Version 8.0 or higher is needed, as we rely on the multi-objective optimization
features introduced in 8.0.

If the environment variable `GUROBI_HOME` is set (and correctly points to your Gurobi installation), Gurobi should be found automatically.
Otherwise you have to set `GUROBI_CXX_LIB` and `GUROBI_LIB` to the libraries `libgurobi_c++.a` and `libgurobi95.so` respectively and `GUROBI_INCLUDE_DIR` to the include directory.

### CPLEX
Version 12.9 or higher is needed, as we rely on the multiobjective optimization
features introduced in in 12.9.

If CPLEX is installed in a default location (`/opt/ibm/ILOG/CPLEX_Studio12<9+>` or `C:\Program Files\IBM\ILOG\CPLEX_Studio12<9+>`), it should be found automatically. If it is installed in such a way that the directories `<path to cplex>/concert` and `<path to cplex>/cplex` are present you can set CPLEX to the `<path to cplex>`. Otherwise, you have to set `CONCERT_INCLUDE_DIR` and `CPLEX_INCLUDE_DIR` to the directories containing `ilconcert/iloenv.h` and `ilcplex/ilocplex.h` respectively and `CONCERT_LIB`, `CPLEX_LIB`, and `ILOCPLEX_LIB` to the corresponding static libraries. Furthermore, on UNIX CPLEX needs the libraries `dl` and `pthread`.

## Vertex Enumeration

PaMILO has an onboard vertex enumeration.
As an alternative, PaMILO is also compatible with [cdd](https://people.inf.ethz.ch/fukudak/cdd_home/).
For numerical challenging programs, cdd might be more stable.

Todo Fukuda citen und anschreiben

To use cdd, install it on your system and then simply enable it in cmake.
If installed in a default location, it should be found automatically.
Otherwise, set `CDD_INCLUDE_PATH` to `<path to folder of cdd.>` and `CDD_LIB` as `<path to cdd>/libcdd.so`.

## Usage
`<pamilo_cli> [<parameters>] <instance>`


with ` <pamilo_cli> = ./pamilo_cli` on Linux/UNIX and `<pamilo_cli> = Release\pamilo_cli.exe` on Windows.

Important optional parameters are:

`-o <output>` The basename for all output files. This defaults to `<instance>`
on UNIX and an empty string on Windows.


`-e <epsilon>` Epsilon to be used in floating point calculations.


`<pamilo_cli> -h` shows all available parameters.


For a list of all parameters with short explanations, use `<pamilo_cli> -h`.

## Instance format
The format of an instance needs to be supported by Gurobi. ([Supported formats](https://www.gurobi.com/documentation/9.5/refman/model_file_formats.html)).
All multi-objective parameters (Priority, Weight, AbsTol, and RelTol) are overwritten by PaMILO.
In some formats (for example `.lp` for Gurobi) the parser still needs those fields, just insert some dummy values in this case.
If the demand is there, this could be changed so that tolerances are taken into account.


An easily human readable format is the `.lp` format. An in depth description of the `.lp` format can be found for [Gurobi](https://www.gurobi.com/documentation/9.5/refman/lp_format.html) and in the [IBM Knowledge Center](https://www.ibm.com/support/knowledgecenter/SSSA5P_20.1.0/ilog.odms.cplex.help/CPLEX/FileFormats/topics/LP.html).


## Output
A list of all extreme points is printed to the command line. Every line is one
point, where the coordinates are separated by spaces. The last line is the total number of extreme points found.
The file `<output>_[gurobi|cplex]` contains all logs from the respective solver.
`<output>_log` logs the running time, the vertex enumeration time and the time
the MILP solver took. If you want anything else to be logged, let us know.
`<output>_sol` stores all extreme points with their corresponding variable assignment.
The format is json, every solution is stored as an object in an array, with one solution looking like:
```
{ "values" : [ <obj1>, <obj2>, ... ], "variables : { "<var1>":<val1>, "<var2>":<val2>, ... } }
```
The `variables` object stores only variables with a value other than 0.
Additionally, a single field `solutionCount` at the end shows the number of extreme points found.

If you prefer any other solution format let us know.

## Notes on MIQP

* Solvers only allow quadratic terms in constraints, but not in objectives.
  You can still solve problems with quadratic terms in objectives the following way:

  1. Add a new continuous variable $z_i$ for each objective $f_i(x)$
  2. Add a new constraint $f_i(x)=z_i$ for each objective
  3. Replace each objective function $f_i(x)$ by $z_i$

* If you attempt to solve a non-convex QP in Gurobi, you need to add the flag `--non-con`.
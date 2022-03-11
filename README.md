# PILP-Solver
Implementation of the dual Benson algorithm to computes the breakpoints (and facets) of a parametric MILP in incremental polynomial time.

The only assumption on the input instances is that an ideal point exists.

The implementation is based on the paper:

Bökler, Fritz and Mutzel, Petra: Output-Sensitive Algorithms for Enumerating the Extreme Nondominated Points of Multiobjective Combinatorial Optimization Problems, Algorithms-ESA, LNCS 9294, Springer, DOI 10.1007/978-3-662-48350-3_25, 2015.

For support please contact Mirko H. Wagner (mirwagner@uos.de).

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
the CMake section below. I have not testet CDD on Windows.

 - Start Visual Studio and build pamilo_cli in Release mode.
 - Go to `<path to pamilo>\Release` and run `pamilo_cli.exe` (for detailed
usage see below).


## CMake
Version 3.18 or higher is needed.

## CPLEX
Version 12.9 or higher is needed, as we rely on the multiobjective optimization
features introduced in in 12.9.


If CPLEX is installed in a default location (`/opt/ibm/ILOG/CPLEX_Studio12<9+>` or `C:\Program Files\IBM\ILOG\CPLEX_Studio12<9+>`), it should be found automaticly. If it is installed in such a way that the directories `<path to cplex>/concert` and `<path to cplex>/cplex` are present you can set CPLEX to the `<path to cplex>`. Otherwise you have to set `CONCERT_INCLUDE_DIR` and `CPLEX_INCLUDE_DIR` to the directories containing `ilconcert/iloenv.h` and `ilcplex/ilocplex.h` respectively and `CONCERT_LIB`, `CPLEX_LIB`, and `ILOCPLEX_LIB` to the corresponding static libraries. Furthermore, on UNIX `CPLEX` needs the libraries `dl` and `pthread`.

## CDD
For four or more objectives the vertex enumeration can be sped up by using the `cddlib`. This is done by activating the CMake flag `USE_CDD` and (if they are not found automatically) providing `libcdd.a` in `CDD_LIB` and the directory containing `cdd.h` in `CDD_INCLUDE_PATH`. The cddlib can be found at [github.com/cddlib/cddlib](https://github.com/cddlib/cddlib).

## Usage
`<pamilo_cli> [<parameters>] <instance>`


with ` <pamilo_cli> = ./pamilo_cli` on Linux/UNIX and `<pamilo_cli> = Release\pamilo_cli.exe` on Windows.


Important parameters are:


`-o <output>` The basename for all output files. This defaults to `<instance>`
on UNIX and an empty string on Windows.


`-e <epsilon>` Epsilon to be used in floating point calculations.


`<pamilo_cli> -h` shows all available parameters.

## Instance format
An instance is to be in CPLEX `.lp` format with multiple objective functions.
All multiobjective parameters (Priority, Weight, AbsTol, and RelTol) are ignored.
If the demand is there this could be changed for the tolerances.


An indepth description of the CPLEX `.lp` format can be found in the [IBM Knowledge Center](https://www.ibm.com/support/knowledgecenter/SSSA5P_20.1.0/ilog.odms.cplex.help/CPLEX/FileFormats/topics/LP.html).



## Output
A list of all extreme points is printed to the command line. Every line is one
point, where the coordinates are separated by spaces.
`<output>_log` logs the running time, the vertex enumeration time and the time
that CPLEX took. If you want anything else to be logged, let us know.
`<output>_sol` stores all extreme points with their corresponding variable assignment. 
The format is json, every solution is stored as an object in an array, with one solution looking like:
```
{ "values" : [ <obj1>, <obj2>, ... ], "variables : { "<var1>":<val1>, "<var2>":<val2>, ... } }
```
The "variables" object stores only variables with a value other than o. If you prefer any other solution format let us know.


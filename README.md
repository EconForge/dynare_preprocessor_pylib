<a name="logo"/>
<div align="center">
<a href="https://www.dynare.org/" target="_blank">
<img src="https://www.dynare.org/assets/images/logo/dlogo.svg" alt="Dynare Logo"></img>
</a>
</div>

# Dynare Preprocessor

The Dynare Preprocessor defines the Dynare model language. It takes in a `.mod`
file, computes the derivatives of the model represented therein, and produces
MATLAB/Octave, Julia, or JSON output.

# Dynare Preprocessor Pylib

This repository contains a lightly modified version of the Dynare preprocessor that makes it usable as a Python library. This was done by replacing most `exit` calls in the code by exceptions, and by adding a wrapper class named `DynareModel` with `pybind11` python bindings. See below for details about both. A meson option has also been added to choose between building the Dynare preprocessor as a binary or as a library. The value of this option can be modified by editing `meson.options` file.

## Installation
`dynare-preprocessor-pylib` is available in conda-forge as well as emscripten-forge.
It can be installed with [micromamba](https://github.com/mamba-org/micromamba-releases) for instance as follows:
```
micromamba install -c https://repo.prefix.dev/conda-forge dynare-preprocessor-pylib
```

## Details about the python bindings

### Interface of DynareModel

#### Initialization
An instance of the `DynareModel` class is initialized by passing a mod file in string format, as well as two optional arguments, the first is the derivation order for variables (defaults to 1) and the second is the derivation order for parameters (defaults to 0).

#### Attributes

- `json_string` (_str_): contains the preprocessor's json output
- `endogenous`, `exogenous`, `exogenous_det` and `parameters` (_list[str]_): lists of symbols
- `equations` (_list[str]_): list of model equations in string form
- `context` (_dict[str,float]_): dictionary of values indexed by the model symbols, filled using steady state then initval blocks, uninitialized values are set to 0
- `covariances` (_dict[tuple[str,str],float_): contains the variances and covariances of the exogenous variables declared in Shocks block
- `trajectories` (_dict[str,list[tuple[int,int,float]]]_): contains the trajectories of exogenous variables declared in surpise shocks block, encoded as lists of elements of the form (p1, p2, v) signifying that the variable in question takes on the value v in periods p1 to p2.
- `symbol_info` (_list[tuple[SymbolType,int,int]]_): maps derivation IDs (currently only appear in the output of the `derivatives` function) to symbol type, (type specific) symbol number and lag. `SymbolType` is the same enum type as in the C++ code base, it most notably contains the values `endogenous`, `exogenous`, `exogenous_det` and `parameter`

#### Methods
- `residuals`:
  - Input:
    - `endo_future` (_list[float]_): list of values of endogenous variables at time t+1, given in the same order as the `endogenous` attribute
    - `endo_present` (_list[float]_): list of values of endogenous variables at time t, given in the same order as the `endogenous` attribute
    - `endo_past` (_list[float]_): list of values of endogenous variables at time t-1, given in the same order as the `endogenous` attribute
    - `exo` (_list[float]_): list of values of exogenous variables at time t, given in the same order as the `exogenous` field
    - `exo_det` (_list[float]_): list of values of deterministic exogenous variables at time t, given in the same order as the `exogenous_det` attribute
    - `params` (_list[float]_): list of values of parameters, given in the same order as the `parameters` attribute
  - Output (_list[float]_): list of residuals at time t


- `jacobians`:
  - Input:
    - `endo_future` (_list[float]_): list of values of endogenous variables at time t+1, given in the same order as the `endogenous` attribute
    - `endo_present` (_list[float]_): list of values of endogenous variables at time t, given in the same order as the `endogenous` attribute
    - `endo_past` (_list[float]_): list of values of endogenous variables at time t-1, given in the same order as the `endogenous` attribute
    - `exo` (_list[float]_): list of values of exogenous variables at time t, given in the same order as the `exogenous` field
    - `exo_det` (_list[float]_): list of values of deterministic exogenous variables at time t, given in the same order as the `exogenous_det` attribute
    - `params` (_list[float]_): list of values of parameters, given in the same order as the `parameters` attribute
  - Output (_list[dict[tuple[int,int], float]]_): list of jacobians of equations with regards to endo_future, endo_present, endo_past, exo, exo_det and params vectors respectively evaluated at time t, each jacobian has the type dict[tuple[int,int], float]. The first index is the equation number and the second is the symbol number. For example, `parameter_jacobian[(0,0)]` is supposed to be the partial derivative of the first residual w.r.t. the first parameter and it exists as long as the partial derivative in equation is not identically equal to zero.

- `derivatives`:
  - Input:
    - `endo_future` (_list[float]_): list of values of endogenous variables at time t+1, given in the same order as the `endogenous` attribute
    - `endo_present` (_list[float]_): list of values of endogenous variables at time t, given in the same order as the `endogenous` attribute
    - `endo_past` (_list[float]_): list of values of endogenous variables at time t-1, given in the same order as the `endogenous` attribute
    - `exo` (_list[float]_): list of values of exogenous variables at time t, given in the same order as the `exogenous` field
    - `exo_det` (_list[float]_): list of values of deterministic exogenous variables at time t, given in the same order as the `exogenous_det` attribute
    - `params` (_list[float]_): list of values of parameters, given in the same order as the `parameters` attribute
  - Output (_list[list[tuple[list[int], float]]]_): list where the kth element is a sparse matrix in COO format containing all k-th order derivatives of the model equations. Each such matrix is a list of pairs: the first element of each element of the list contains indices and the second contains the value of the derivatives at these indices. Each list of indices starts with the equation number as a first index, followed by the derivation ids of the variables with regards to which the derivation was done in non-decreasing order. In order to get information about a variable using its derivation id, refer to the `symbol_info` list.

### Exception handling
Four types of exceptions are exposed to the Python side:
- `PreprocessorException` is the base exception. It is used virtually anywhere an exit was used (except `DynareMain` and `Configuration` which shouldn't be needed if the preprocessor is used as a library) and is intended as a temporary catchall until the tedious task of defining more specific Exceptions and handling them is undertaken.
- `ParserException`: inherits from `PreprocessorException`. It's used for errors stemming from the parsing process.
- `UnsupportedFeatureException`: inherits from `PreprocessorException`. It's used for features that exist in the dynare-preprocessor binary but that do not work as of yet in the python library.
- `EvalException`: inherits from `PreprocessorException`. It's used for errors arising from failed evaluations (Divide by zero, uninitialized variable, etc…)

# License

Most of the source files are covered by the GNU General Public Licence version
3 or later. There are some exceptions, see the respective file headers.

<a name="logo"/>
<div align="center">
<a href="https://www.dynare.org/" target="_blank">
<img src="https://www.dynare.org/assets/images/logo/dlogo.svg" alt="Dynare Logo"></img>
</a>
</div>

# Dynare Preprocessor

The **Dynare Preprocessor** defines and parses the Dynare model specification language. It parses `.mod` files, runs macroprocessing directives (`@#define`, `@#if`, `@#for`, etc.), constructs the symbolic model representation, computes analytical derivatives (Jacobian, Hessian, and higher-order tensors), and produces simulation/estimation drivers and evaluators for **MATLAB/Octave**, **Julia**, **Python**, or structured **JSON AST** dumps.

In addition to the standalone CLI binary, it provides a high-performance native Python library (`dynare_preprocessor`) powered by [nanobind](https://github.com/wjakob/nanobind). The Python library allows you to parse `.mod` files (from disk or in-memory strings), introspect model equations and symbols, and evaluate residuals, Jacobians, and higher-order derivatives directly via NumPy arrays without intermediate code generation.

> [!NOTE]
> **Project Status & Upstream Integration**
> This repository is **not** an independent fork of the Dynare Preprocessor. It is developed with the blessing and active collaboration of the Dynare team, and all enhancements and packaging features will be integrated directly into the official [Dynare project](https://git.dynare.org/Dynare/preprocessor).
>
> Please note that this is an **early preview**: the codebase is in active development, and the Python library API is subject to change as upstream integration proceeds.

---

## Table of Contents

1. [Quick Start with Pixi](#1-quick-start-with-pixi)
2. [Python Library Guide & API Reference](#2-python-library-guide--api-reference)
   - [Basic Example](#basic-example)
   - [DynareModel Class](#dynaremodel-class)
   - [JacobianBlocks & Derivatives](#jacobianblocks--derivatives)
   - [Exception Handling & Source Locations](#exception-handling--source-locations)
3. [Installation & Packaging](#3-installation--packaging)
   - [Conda / Prefix.dev](#conda--prefixdev)
   - [Pip / Distributable Wheels](#pip--distributable-wheels)
4. [Standalone C++ CLI & Meson Build](#4-standalone-c-cli--meson-build)
5. [WebAssembly / Emscripten Packaging](#5-webassembly--emscripten-packaging)
6. [License](#license)

---

## 1. Quick Start with Pixi

The repository uses [Pixi](https://pixi.sh) to manage reproducible cross-platform development environments across Linux, macOS (Apple Silicon and Intel), and Windows.

### Full Development Environment (CLI + Python Library)

```bash
# Install dependencies into the default development environment
pixi install

# Compile the preprocessor CLI and nanobind Python extension
pixi run compile

# Run Python pytest test suite
pixi run test

# Run all test suites via Meson (C++ CLI tests + Python tests)
pixi run test-all
```

### Standalone C++ CLI Only (`cli` environment)

To build and test the standalone `dynare-preprocessor` CLI without any Python or nanobind dependencies:

```bash
# Configure the minimal CLI build directory
pixi run -e cli setup-cli

# Compile the CLI binary
pixi run -e cli compile-cli

# Run CLI tests
pixi run -e cli test-cli
```

---

## 2. Python Library Guide & API Reference

The Python package `dynare_preprocessor` provides zero-copy C++ AST evaluation using [nanobind](https://github.com/wjakob/nanobind) and NumPy arrays.

### Basic Example

```python
import numpy as np
import dynare_preprocessor as dp

# Load model from a .mod file on disk or an in-memory string
mod_text = """
var c k y;
varexo e;
parameters alpha beta delta;

alpha = 0.33;
beta  = 0.99;
delta = 0.025;

model;
  c + k = y + (1 - delta) * k(-1);
  y = k(-1)^alpha;
  1/c = beta * (1/c(+1)) * (alpha * y(+1)/k + 1 - delta);
end;

initval;
  k = 10.0;
  c = 1.0;
  y = 1.2;
  e = 0.0;
end;
"""

model = dp.DynareModel(mod_text, derivs_order=1)

# Inspect model declarations
print("Endogenous:", model.endogenous)    # ['c', 'k', 'y']
print("Exogenous:", model.exogenous)      # ['e']
print("Parameters:", model.parameters)    # ['alpha', 'beta', 'delta']
print("Equations count:", len(model.equations))

# Extract initial values from the model context
y_ss = np.array([model.context[v] for v in model.endogenous])
e_ss = np.array([model.context[v] for v in model.exogenous])
ed_ss = np.array([model.context[v] for v in model.exogenous_det])
p_ss = np.array([model.context[v] for v in model.parameters])

# 1. Evaluate dynamic model residuals F(y_{t+1}, y_t, y_{t-1}, e_t, params)
res = model.residuals(y_ss, y_ss, y_ss, e_ss, ed_ss, p_ss)
print("Residuals shape:", res.shape)  # (3,)

# 2. Evaluate structured Jacobian blocks (lead, current, lag, exo, exo_det, params)
blocks = model.jacobian_blocks(y_ss, y_ss, y_ss, e_ss, ed_ss, p_ss)
print("dF/dy_{t+1} (lead):\n", blocks.lead)
print("dF/dy_t     (curr):\n", blocks.curr)
print("dF/dy_{t-1} (lag):\n",  blocks.lag)
print("dF/de_t     (exo):\n",  blocks.exo)

# Blocks can also be unpacked directly
lead, curr, lag, exo, exo_det, params = blocks

# 3. Evaluate static (steady-state) residuals and Jacobian
s_res = model.static_residuals(y_ss, e_ss, ed_ss, p_ss)
s_jac = model.static_jacobian(y_ss, e_ss, ed_ss, p_ss)
```

### DynareModel Class

#### Constructor

```python
model = dp.DynareModel(
    modfile_content_or_path: str,
    derivs_order: int = 1,
    params_derivs_order: int = 0,
    strict: bool = False
)
```

- `modfile_content_or_path`: Either a filesystem path to a `.mod` file or a raw string containing the model definition.
- `derivs_order`: Maximum derivation order with respect to variables (default `1`).
- `params_derivs_order`: Maximum derivation order with respect to parameters (default `0`).
- `strict`: If `True`, treat undeclared variables as immediate fatal parsing errors (default `False`).

#### Attributes

| Attribute | Type | Description |
|---|---|---|
| `endogenous` | `list[str]` | List of endogenous variable names $y$. |
| `exogenous` | `list[str]` | List of exogenous shock names $\epsilon$. |
| `exogenous_det` | `list[str]` | List of deterministic exogenous variable names. |
| `parameters` | `list[str]` | List of parameter names $\theta$. |
| `equations` | `list[str]` | Mathematical equations formatted as strings. |
| `context` | `dict[str, float]` | Numerical values populated from assignments, `initval`, or `steady_state_model` blocks. |
| `covariances` | `dict[tuple[str, str], float]` | Variances and covariances declared in `shocks` blocks. |
| `trajectories` | `dict[str, list[tuple[int, int, float]]]` | Deterministic shock paths from surprise/shocks blocks. |
| `lead_lag_incidence` | `list[list[int]]` | Dynare lead-lag incidence matrix across lags $(-1)$, current $(0)$, and leads $(+1)$. |
| `max_endo_lag` | `int` | Maximum lag depth across endogenous variables. |
| `max_endo_lead` | `int` | Maximum lead depth across endogenous variables. |
| `symbol_info` | `list[tuple[SymbolType, int, int]]` | Maps derivation IDs to `(SymbolType, symbol_id, lag)`. |
| `json_string` | `str` | Complete JSON AST representation of the parsed model. |

#### Methods

All evaluator methods accept inputs as either contiguous 1D NumPy arrays (`np.ndarray`) or standard Python sequences of floats:

- `residuals(endo_future, endo_present, endo_past, exo, exo_det, params) -> np.ndarray`
  Evaluates dynamic model residuals $F(y_{t+1}, y_t, y_{t-1}, \epsilon_t, \epsilon_{det,t}, \theta)$ as a 1D NumPy array of shape `(n_eq,)`.
- `jacobian_blocks(endo_future, endo_present, endo_past, exo, exo_det, params) -> JacobianBlocks`
  Returns structured analytical Jacobian blocks as 2D NumPy arrays:
  - `lead` ($n_{eq} \times n_{endo}$): $\partial F / \partial y_{t+1}$
  - `curr` ($n_{eq} \times n_{endo}$): $\partial F / \partial y_t$
  - `lag` ($n_{eq} \times n_{endo}$): $\partial F / \partial y_{t-1}$
  - `exo` ($n_{eq} \times n_{exo}$): $\partial F / \partial \epsilon_t$
  - `exo_det` ($n_{eq} \times n_{exo\_det}$): $\partial F / \partial \epsilon_{det,t}$
  - `params` ($n_{eq} \times n_{params}$): $\partial F / \partial \theta$
- `jacobian(endo_future, endo_present, endo_past, exo, exo_det, params) -> np.ndarray`
  Evaluates the combined dynamic Jacobian as a 2D NumPy array across dynamic incidence columns and shocks.
- `static_residuals(endo, exo, exo_det, params) -> np.ndarray`
  Evaluates static (steady-state) residuals $F(y, \epsilon, \epsilon_{det}, \theta)$ as a 1D NumPy array of shape `(n_eq,)`.
- `static_jacobian(endo, exo, exo_det, params) -> np.ndarray`
  Evaluates static (steady-state) Jacobian $\partial F / \partial y$ as a 2D NumPy array of shape `(n_eq, n_endo)`.
- `derivatives(endo_future, endo_present, endo_past, exo, exo_det, params) -> list[list[tuple[list[int], float]]]`
  Evaluates higher-order analytical derivatives in sparse coordinate (COO) format. The $k$-th element contains the list of non-zero entries for the $k$-th order tensor: `([eq, var_1, ..., var_k], value)`.

### Exception Handling & Source Locations

The library provides structured C++ exceptions translated natively into Python exceptions:

```python
import dynare_preprocessor as dp

try:
    model = dp.DynareModel("var x; model; x = ; end;")
except dp.ParserException as e:
    print(f"Parse error: {e.message}")
    if e.location:
        print(f"Location: {e.location.filename}:{e.location.line}:{e.location.column}")
```

#### Exception Hierarchy

- `DynareException`: Base class for all preprocessor exceptions.
  - `SourceFileException`: Base class for exceptions tied to source code positions.
    - `ParserException`: Syntax errors, grammar violations, or undeclared symbols.
      - Properties: `location` (`SourceLocation`), `line`, `column`, `filename`, `undeclared_variables`.
    - `MacroException`: Macroprocessor expansion errors (`@#error`, invalid directives).
      - Properties: `location` (`SourceLocation`), `line`, `column`, `backtrace`.
  - `ModelSemanticException`: Mathematical or semantic errors in the model definition (e.g. equation-variable count mismatch, invalid leads/lags).
    - Properties: `equation_number`, `equation_lineno`, `equation_tag`, `symbol_name`, `location`, `line`.
    - `EquationException`: Specific equation error.
  - `StatementException`: Incompatible statements or options (e.g. mixing perfect foresight with stochastic simulation).
    - Properties: `statement_name`, `option_name`.
  - `FileIOException`: File read/write failures.
    - Properties: `path`, `action`.
  - `EvaluationException`: Numerical evaluation errors during steady state evaluation (e.g. division by zero, domain errors).

---

## 3. Installation & Packaging

### Conda / Prefix.dev

Prebuilt packages for Linux (`linux-64`), macOS (`osx-arm64`, `osx-64`), Windows (`win-64`), and WebAssembly (`emscripten-wasm32`) are distributed via `conda-forge` and `prefix.dev/econforge`:

```bash
# Install via conda or mamba
conda install -c https://repo.prefix.dev/econforge -c conda-forge dynare-preprocessor-pylib

# Or add to an existing Pixi project
pixi add --channel https://repo.prefix.dev/econforge dynare-preprocessor-pylib
```

### Pip / Distributable Wheels

The Python package is PEP 517 compliant using [meson-python](https://meson-python.readthedocs.io/):

```bash
# Local editable installation
pip install --no-build-isolation -e .

# Or build standalone wheels (.whl) for distribution
pixi run wheel
# The wheel will be generated in dist/
pip install dist/dynare_preprocessor-*.whl
```

---

## 4. Standalone C++ CLI & Meson Build

If you are not using Pixi and have system dependencies installed (`meson >= 1.3.0`, `ninja`, modern C++20 compiler, `boost`, `flex`, `bison`):

### Building the Standalone CLI Only

```bash
meson setup build -Dbuild_cli=enabled -Dbuild_library=disabled -Dbuild_doc=false
meson compile -C build

# Run the preprocessor on a .mod file
./build/src/dynare-preprocessor example.mod

# Dump JSON AST to stdout
./build/src/dynare-preprocessor example.mod json=parse jsonstdout
```

### Building Both CLI and Python Extension

```bash
meson setup build -Dbuild_cli=enabled -Dbuild_library=enabled -Dbuild_doc=false
meson compile -C build
PYTHONPATH=build/src pytest tests
```

---

## 5. WebAssembly / Emscripten Packaging

The Python library can be compiled to WebAssembly for browser runtimes such as [Pyodide](https://pyodide.org), [JupyterLite](https://jupyterlite.readthedocs.io/), and [Stlite](https://github.com/whitphx/stlite).

### Building & Testing Packages

```bash
# WebAssembly: Build emscripten-wasm32 package locally via rattler-build
pixi run build-wasm

# WebAssembly: Run tests in a headless browser (Chromium via Playwright + pytester)
pixi run test-wasm

# WebAssembly: Upload package to prefix.dev/econforge
pixi run upload-wasm

# Linux: Build packages across Python 3.11, 3.12, 3.13 variants
pixi run build-linux

# Linux: Upload packages to prefix.dev/econforge
pixi run upload-linux
```

### Using in WebAssembly Runtimes (JupyterLite / Pyodide / Pixi)

Include `econforge` and `emscripten-forge-4x` channels in your `pixi.toml` or environment configuration:

```toml
[workspace]
channels = [
    "https://repo.prefix.dev/econforge",
    "https://repo.prefix.dev/emscripten-forge-4x",
    "conda-forge"
]
platforms = ["emscripten-wasm32"]

[dependencies]
dynare-preprocessor-pylib = ">=0.0.1.dev0"
```

Or create an environment using standard Conda:

```bash
conda create -n wasm-env \
    --platform=emscripten-wasm32 \
    -c https://repo.prefix.dev/econforge \
    -c https://repo.prefix.dev/emscripten-forge-4x \
    -c conda-forge \
    dynare-preprocessor-pylib
```

---

## License

Most of the source files are covered by the GNU General Public License version 3 or later. There are some exceptions; see the respective file headers.

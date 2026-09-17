<a name="logo"/>
<div align="center">
<a href="https://www.dynare.org/" target="_blank">
<img src="https://www.dynare.org/assets/images/logo/dlogo.svg" alt="Dynare Logo"></img>
</a>
</div>

# Dynare Preprocessor

The Dynare Preprocessor defines the Dynare model language. It parses `.mod`
files, performs macroprocessing, constructs the symbolic model representation,
computes analytical derivatives (Jacobian, Hessian, and higher order), and produces
simulation/estimation drivers and evaluators for MATLAB/Octave, Julia, Python, or JSON AST dumps.

It also provides a native Python library (`dynare_preprocessor`) powered by [nanobind](https://github.com/wjakob/nanobind)
for programmatically loading `.mod` files and evaluating model residuals, Jacobians, and higher-order derivatives directly in Python.

---

## 1. Quick Start with Pixi

The repository is configured with [Pixi](https://pixi.sh) to provide reproducible cross-platform development environments across Linux, macOS, and Windows.

### Python Library & Full Development (Default)

```bash
# Install dependencies into default environment
pixi install

# Compile the preprocessor CLI and Python extension module
pixi run compile

# Run the full test suite (Pytest)
pixi run test

# Run all tests via Meson (CLI exception tests + Python tests)
pixi run test-all
```

### Standalone C++ CLI Only (`cli` environment)

To build and test the standalone `dynare-preprocessor` CLI binary without any Python or nanobind dependencies:

```bash
# Setup the minimal C++ build directory
pixi run -e cli setup-cli

# Compile the CLI executable
pixi run -e cli compile-cli

# Run CLI tests
pixi run -e cli test-cli
```

---

## 2. Python Package & Wheel Installation

The Python package is standard PEP 517 compliant using [meson-python](https://meson-python.readthedocs.io/).

### Editable / Local Installation

Inside your Python or Pixi environment:

```bash
# Editable install (in-place development)
pip install --no-build-isolation -e .

# Or using the Pixi task:
pixi run install-editable
```

### Building Distributable Wheels (`.whl`)

To build binary wheels for distribution (e.g. for PyPI or local installation):

```bash
# Using Pixi (handles platform flags automatically):
pixi run wheel

# Or directly with python -m build:
python -m build --wheel --no-isolation --skip-dependency-check
```

The resulting wheel is saved in `dist/` (e.g. `dist/dynare_preprocessor-0.0.1.dev0-cp312-cp312-linux_x86_64.whl`) and can be installed via pip:

```bash
pip install dist/dynare_preprocessor-*.whl
```

---

## 3. Standalone Meson Build (Manual)

If you are not using Pixi and have system dependencies installed (C++20 compiler, Boost, Flex, Bison, Meson, Ninja):

### Building the CLI only:
```bash
meson setup build -Dbuild_cli=enabled -Dbuild_library=disabled -Dbuild_doc=false
meson compile -C build
./build/src/dynare-preprocessor example.mod
```

### Building both CLI and Python library:
```bash
meson setup build -Dbuild_cli=enabled -Dbuild_library=enabled -Dbuild_doc=false
meson compile -C build
PYTHONPATH=build/src pytest tests
```

---

## 4. WebAssembly / Emscripten Build & Testing

### WebAssembly Packaging (`emscripten-wasm32`)

```bash
# Build the emscripten-wasm32 package locally (~1 min)
pixi run build-wasm

# Run the test suite inside a headless browser (Chromium via pytester)
pixi run test-wasm

# Upload the package to the econforge channel on prefix.dev
pixi run upload-wasm
```

The resulting package will be generated under `output/emscripten-wasm32/dynare-preprocessor-pylib-*.conda`.

### Linux Packaging (`linux-64`)

```bash
# Build packages for Python 3.11, 3.12, 3.13 variants
pixi run build-linux

# Upload Linux packages to the econforge channel on prefix.dev
pixi run upload-linux
```

The resulting packages will be generated under `output/linux-64/dynare-preprocessor-pylib-*.conda`.

### Consuming the Package (JupyterLite / Pyodide / Pixi)

To install the WebAssembly package from `econforge`, include the `emscripten-forge-4x` channel for WebAssembly runtime dependencies (`emscripten-abi`, `python`, `numpy`):

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

Or via `micromamba`:
```bash
micromamba create -n wasm-env \
    --platform=emscripten-wasm32 \
    -c https://repo.prefix.dev/econforge \
    -c https://repo.prefix.dev/emscripten-forge-4x \
    -c conda-forge \
    dynare-preprocessor-pylib
```

---

## License

Most of the source files are covered by the GNU General Public License version
3 or later. There are some exceptions; see the respective file headers.

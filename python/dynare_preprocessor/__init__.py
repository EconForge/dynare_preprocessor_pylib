"""
Dynare Preprocessor Python Library.

This package provides Python bindings to the Dynare Preprocessor, allowing
programmatic parsing of Dynare .mod files, evaluation of dynamic and static
model residuals, Jacobians, and higher-order analytical derivatives.
"""

try:
    from ._dynare_preprocessor import (
        DynareException,
        DynareModel,
        EquationException,
        EvalException,
        EvaluationException,
        FileIOException,
        InternalCompilerException,
        JacobianBlocks,
        MacroException,
        ModelSemanticException,
        ParserException,
        PreprocessorException,
        SourceFileException,
        StatementException,
        SymbolType,
        UnsupportedFeatureException,
    )
except ModuleNotFoundError:
    # Fallback for uninstalled in-tree builds where the shared object resides in build/src
    from _dynare_preprocessor import (  # type: ignore[no-redef]
        DynareException,
        DynareModel,
        EquationException,
        EvalException,
        EvaluationException,
        FileIOException,
        InternalCompilerException,
        JacobianBlocks,
        MacroException,
        ModelSemanticException,
        ParserException,
        PreprocessorException,
        SourceFileException,
        StatementException,
        SymbolType,
        UnsupportedFeatureException,
    )

__version__ = "8.0.0.dev0"

__all__ = [
    "DynareModel",
    "JacobianBlocks",
    "SymbolType",
    "DynareException",
    "SourceFileException",
    "ParserException",
    "MacroException",
    "ModelSemanticException",
    "EquationException",
    "StatementException",
    "FileIOException",
    "EvaluationException",
    "InternalCompilerException",
    "PreprocessorException",
    "UnsupportedFeatureException",
    "EvalException",
    "__version__",
]

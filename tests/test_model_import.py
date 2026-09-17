import os
import sys
import pytest

build_src = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "build", "src"))
if build_src not in sys.path:
    sys.path.insert(0, build_src)

from dynare_preprocessor import (
    DynareModel,
    DynareException,
    ParserException,
    ModelSemanticException,
    StatementException,
    EvaluationException,
)

files = [
    "modfile.mod",
    "ramst.mod",
    "NK_baseline.mod",
    "fs2000.mod",
    "Ramsey_steady_file.mod",
    "example1.mod",
    "model_KR2000_IRF.mod",
    "example1_reporting.mod",
    "model_KR2000_STAT.mod",
    "bkk.mod",  # Contains macro directives (@#) -> succeeds now!
    "RBC.mod",
    "fs2000_nonstationary.mod",
    "Ramsey_Example.mod",  # Contains macro directives (@#) -> succeeds now!
    "Gali_2015.mod",
    "example3.mod",
    "agtrend.mod",  # Contains macro directives (@#) -> succeeds now!
    "example2.mod",
    "Occbin_example.mod",
]

# Modfiles expected to raise specific exceptions
expected_exceptions = {
    "modfile.mod": ParserException,  # Symbol declared twice with different types
    "example1.mod": EvaluationException,  # Shock block evaluation error
    "example1_reporting.mod": EvaluationException,  # Shock block evaluation error
    "example3.mod": ModelSemanticException,  # External function in steady state
    "ramst.mod": StatementException,  # Conflict: perfect foresight and stochastic
}


@pytest.mark.parametrize("filename", files)
def test_modfile_import(filename):
    path = os.path.join(os.path.dirname(__file__), "modfiles", filename)
    with open(path, encoding="utf-8") as f:
        mod_string = f.read()

    if filename in expected_exceptions:
        exc_type = expected_exceptions[filename]
        with pytest.raises(exc_type):
            DynareModel(mod_string)
    else:
        model = DynareModel(mod_string)
        assert len(model.endogenous) > 0
        assert len(model.equations) > 0


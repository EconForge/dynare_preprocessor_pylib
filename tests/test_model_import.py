from dynare_preprocessor import DynareModel, PreprocessorException, ParserException, UnsupportedFeatureException
import pytest

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
    "bkk.mod",
    "RBC.mod",
    "fs2000_nonstationary.mod",
    "Ramsey_Example.mod",
    "Gali_2015.mod",
    "example3.mod",
    "agtrend.mod",
    "example2.mod",
    "Occbin_example.mod",
]

parser_exception = [
    "modfile.mod", # symbol declared twice with different types
    "bkk.mod", # character unrecognized by lexer `@`
    "Ramsey_Example.mod", # character unrecognized by lexer `@`
    "agtrend.mod", # character unrecognized by lexer `@`
]

unsupported_exception = [
    "example1.mod", # unsupported native statement
    "example1_reporting.mod", # unsupported native statement
    "example3.mod", # external steady state helper
]

generic_exception = [
    "ramst.mod", # mixed perfect foresight context with stochastic context
]

@pytest.mark.parametrize("filename", files)
def test_modfile_import(filename):
    try:
        f = filename
        filename = "tests/modfiles/" + f
        mod_string = open(filename).read()
        model = DynareModel(mod_string)
    except ParserException:
        assert(f in parser_exception)
    except UnsupportedFeatureException:
        assert(f in unsupported_exception)
    except PreprocessorException:
        assert(f in generic_exception)
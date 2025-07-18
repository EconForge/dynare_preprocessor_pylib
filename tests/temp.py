#!/usr/bin/env python
from dynare_preprocessor import DynareModel
f = "tests/modfiles/Occbin_example.mod"
mod_string = open(f).read()
model = DynareModel(mod_string)

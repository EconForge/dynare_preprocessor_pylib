# from ctypes import cdll
# mylib = cdll.LoadLibrary("./libdynare_preprocessor_lib.so")
# mylib.just_try("am a your father.")


import dynare_preprocessor

txt = open('/home/work/dyno.py/examples/modfiles/example1.mod').read()

from time import time
import json

for mode, modename in enumerate(["nothing", "parsing", "checking", "transforming", "computing"]):
    t0 = time()
    res = dynare_preprocessor.preprocess(txt, mode)
    data = json.loads(res)
    t1 = time()
    print(f"Steps up to {modename} take {t1-t0} seconds")

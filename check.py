# from ctypes import cdll
# mylib = cdll.LoadLibrary("./libdynare_preprocessor_lib.so")
# mylib.just_try("am a your father.")


import dynare_preprocessor

txt = open('/home/pablo/Source/dynare/examples/example1.mod').read()

import time
res = dynare_preprocessor.preprocess(txt)

with open("out.json", "w") as f:
    f.write(res)

import json



t1 = time.time()
dynare_preprocessor.preprocess(txt)
data = json.loads(res)
t2 = time.time()

print("Elapsed : ", t2-t1)
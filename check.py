# from ctypes import cdll
# mylib = cdll.LoadLibrary("./libdynare_preprocessor_lib.so")
# mylib.just_try("am a your father.")


import dynare_preprocessor

txt = open('/home/work/dyno.py/examples/modfiles/example1.mod').read()

res = dynare_preprocessor.preprocess(txt, 4)
import json
data = json.loads(res)
# print(data)

with open("out.json", "w") as f:
    f.write(res)

print("Last line of python program executed")
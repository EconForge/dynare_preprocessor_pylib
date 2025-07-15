import numpy as np

def reexecute_if_unbuffered():
    """Ensures that output is immediately flushed (e.g. for segfaults).
    ONLY use this at your entrypoint. Otherwise, you may have code be
    re-executed that will clutter your console."""
    import os
    import shlex
    import sys
    if os.environ.get("PYTHONUNBUFFERED") in (None, ""):
        os.environ["PYTHONUNBUFFERED"] = "1"
        argv = list(sys.argv)
        if argv[0] != sys.executable:
            argv.insert(0, sys.executable)
        cmd = " ".join([shlex.quote(arg) for arg in argv])
        sys.stdout.flush()
        os.execv(argv[0], argv)


def traced(func, ignoredirs=None):
    """Decorates func such that its execution is traced, but filters out any
     Python code outside of the system prefix."""
    import functools
    import sys
    import trace
    if ignoredirs is None:
        ignoredirs = ["/usr", sys.prefix]
    tracer = trace.Trace(trace=1, count=0, ignoredirs=ignoredirs)

    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        return tracer.runfunc(func, *args, **kwargs)

    return wrapped


# NOTE: You don't have to trace all of your code. If you can identify a
# single function, then you can just decorate it with this. If you're
# decorating a class method, then be sure to declare these functions above
# it.
# @traced
def main():
    from dynare_preprocessor import DynareModel
    f = "/home/work/dyno.py/examples/modfiles/example1.mod"
    mod_string = open(f).read()
    from time import time
    t0 = time()
    dyn = DynareModel(mod_string)
    t1 = time()
    print("Symbols:")
    print(f"""
    endogenous: {dyn.endogenous}
    exogenous: {dyn.exogenous}
    exogenous deterministic: {dyn.exogenous_det}
    parameters: {dyn.parameters}
    """)
    print("Equations:\n")
    for eq in dyn.equations:
        print(f"    {eq}")
    print("\nCalibration:\n")
    for k,v in dyn.calibration.items():
        print(f"    {k} = {v}")
    print("\nResiduals:\n")
    endo = [dyn.calibration[x] for x in dyn.endogenous]
    exo = [dyn.calibration[x] for x in dyn.exogenous]
    exo_det = [dyn.calibration[x] for x in dyn.exogenous_det]
    params = [dyn.calibration[x] for x in dyn.parameters]
    print(f"    {dyn.dynamic_function(endo,endo,endo,exo,exo_det,params)}")
    print("\nCovariances:\n")
    for [x,y],v in dyn.covariances.items():
        print(f"    Covar({x}, {y}) = {v}")
    print("\nTrajectories:\n")
    for v, traj in dyn.trajectories.items():
        print(f"    Trajectory of exogenous variable {v}")
        for p1,p2,val in traj:
            print(f"        {val} from {p1} to {p2}")
    print("\nJacobians:\n")
    jacobians = dyn.jacobians(endo,endo,endo,exo,exo_det,params)
    n = len(dyn.equations)
    symbtypes = ["future endogenous",
                "present endogenous",
                "past endogenous",
                "exogenous",
                "deterministic exogenous",
                "parameters"]
    symblists = [dyn.endogenous,
                dyn.endogenous,
                dyn.endogenous,
                dyn.exogenous,
                dyn.exogenous_det,
                dyn.parameters]
    for ind, symblist in enumerate(symblists):
        m = len(symblist)
        matrix = np.zeros((n,m))
        for (i,j), v in jacobians[ind].items():
            matrix[i,j] = v
        print(f"    {symbtypes[ind]}\n")
        print(f"    {matrix}\n")
    t2 = time()
    print(f"Elapsed: {t1-t0}")






if __name__ == "__main__":
    reexecute_if_unbuffered()
    main()

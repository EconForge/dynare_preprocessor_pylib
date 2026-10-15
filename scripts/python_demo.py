"""Parse a small RBC model with `dynare_preprocessor` and exercise the Python API.

Loads an RBC model written in Dynare syntax (as an in-memory string), evaluates
the dynamic residuals at the steady state, and validates the analytic Jacobian
blocks against finite-difference approximations.

Usage: pixi run python_demo
"""

import contextlib
import ctypes
import os
import sys
import tempfile
import time

import numpy as np
from rich.console import Console, Group
from rich.panel import Panel
from rich.rule import Rule
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text

from dynare_preprocessor import DynareModel

RBC = """
var c, r, w, n, k, i, y, a;
varexo epsilon;
parameters beta, delta, khi, eta, alpha, rho, nss;

beta=0.985;
delta=0.025;
nss=0.33;
eta=1.0;
alpha=.33;
rho=.95;
khi=(1-alpha)*(1-nss)^eta/nss*(1/beta-1+delta)/(1/beta-1+delta-delta*alpha);

model;
1/c = beta*(r(1)+1-delta)/c(1);
w = khi*c/(1-n)^eta;
k = (1-delta)*k(-1)+i;
y = a*k(-1)^alpha*n^(1-alpha);
log(a) = rho*log(a(-1))+epsilon;
w = (1-alpha)*y/n;
r = alpha*y/k(-1);
y = c+i;
end;

steady_state_model;
a = 1;
r = 1/beta-1+delta;
n = nss;
k = (alpha/r)^(1/(1-alpha))*n;
y = k^alpha*n^(1-alpha);
w = (1-alpha)*y/n;
i = delta*k;
c = y-i;
end;

shocks;
var epsilon; stderr .009;
end;
"""

console = Console()
libc = ctypes.CDLL(None)


@contextlib.contextmanager
def capture_fd_stdout():
    """Capture what the C++ preprocessor writes to file descriptor 1."""
    out = []
    sys.stdout.flush()
    saved = os.dup(1)
    with tempfile.TemporaryFile(mode="w+b") as tmp:
        os.dup2(tmp.fileno(), 1)
        try:
            yield out
        finally:
            sys.stdout.flush()
            libc.fflush(None)  # flush the C stdio buffer before restoring fd 1
            os.dup2(saved, 1)
            os.close(saved)
            tmp.seek(0)
            out.append(tmp.read().decode(errors="replace"))


def show_code(title, code):
    """Print the instructions of a step, syntax-highlighted, before running them."""
    console.print(Rule(f"[bold]{title}[/bold]", style="blue"))
    console.print(Syntax(code.strip(), "python", theme="monokai", background_color="default"))
    console.print()


def steady_state_values(model):
    """Collect the steady-state (context) values of each symbol type."""
    y_ss = np.array([model.context[v] for v in model.endogenous])
    e_ss = np.array([model.context[v] for v in model.exogenous])
    ed_ss = np.array([model.context[v] for v in model.exogenous_det])
    p_ss = np.array([model.context[v] for v in model.parameters])
    return y_ss, e_ss, ed_ss, p_ss


def finite_difference_jacobian(func, x0, delta=1e-7):
    """Compute a numerical Jacobian using central differences."""
    f0 = func(x0)
    m, n = len(f0), len(x0)
    jac = np.zeros((m, n))
    for j in range(n):
        dx = np.zeros(n)
        dx[j] = delta
        jac[:, j] = (func(x0 + dx) - func(x0 - dx)) / (2.0 * delta)
    return jac


def main():
    console.print(
        Panel(
            "[bold red]⚠  Warning: the Python API is [u]not[/u] definitive "
            "and may change.[/bold red]",
            border_style="red",
            expand=False,
        )
    )
    console.print()
    show_code(
        "1. Parse the model",
        """
from dynare_preprocessor import DynareModel

model = DynareModel(RBC, stochastic=True, params_derivs_order=1)   # RBC: .mod source (or a path to a .mod file)
""",
    )
    with capture_fd_stdout() as log:
        t0 = time.perf_counter()
        model = DynareModel(RBC, stochastic=True, params_derivs_order=1)
        elapsed = time.perf_counter() - t0

    summary = Text.assemble(
        ("RBC model", "bold"),
        f": {len(model.endogenous)} endogenous, {len(model.exogenous)} shocks, "
        f"{len(model.parameters)} parameters — parsed in ",
        (f"{1e3 * elapsed:.2f} ms", "bold green"),
    )
    console.print(
        Panel(
            Group(summary, Text(""), Text(log[0].strip(), style="dim")),
            title="[bold]Parsing[/bold]",
            border_style="blue",
            expand=False,
        )
    )

    y_ss, e_ss, ed_ss, p_ss = steady_state_values(model)

    console.print()
    show_code(
        "2. Evaluate residuals at the steady state",
        """
res = model.residuals(y_ss, y_ss, y_ss, e_ss, ed_ss, p_ss)   # F(y_{t+1}, y_t, y_{t-1}, e_t, params)
res.shape[0]   # one residual per equation
""",
    )
    with capture_fd_stdout():
        t0 = time.perf_counter()
        res = model.residuals(y_ss, y_ss, y_ss, e_ss, ed_ss, p_ss)
        elapsed = time.perf_counter() - t0

    console.print(
        Panel(
            Group(
                Text(f"max |residual| = {np.abs(res).max():.3e}"),
                Text(""),
                Text(f"computed in {1e3 * elapsed:.3f} ms", style="dim"),
            ),
            title="[bold]Steady state[/bold]",
            border_style="blue",
            expand=False,
        )
    )

    console.print()
    show_code(
        "3. Validate analytic Jacobians against finite differences",
        """
jb = model.jacobian_blocks(y_ss, y_ss, y_ss, e_ss, ed_ss, p_ss)
# jb.lead, jb.curr, jb.lag, jb.exo, jb.params are 2D NumPy arrays
""",
    )
    with capture_fd_stdout():
        t0 = time.perf_counter()
        jb = model.jacobian_blocks(y_ss, y_ss, y_ss, e_ss, ed_ss, p_ss)
        elapsed = time.perf_counter() - t0

    checks = [
        ("lead", jb.lead, lambda u: model.residuals(u, y_ss, y_ss, e_ss, ed_ss, p_ss), y_ss),
        ("curr", jb.curr, lambda u: model.residuals(y_ss, u, y_ss, e_ss, ed_ss, p_ss), y_ss),
        ("lag", jb.lag, lambda u: model.residuals(y_ss, y_ss, u, e_ss, ed_ss, p_ss), y_ss),
        ("exo", jb.exo, lambda u: model.residuals(y_ss, y_ss, y_ss, u, ed_ss, p_ss), e_ss),
        ("params", jb.params, lambda u: model.residuals(y_ss, y_ss, y_ss, e_ss, ed_ss, u), p_ss),
    ]

    table = Table(title="JacobianBlocks", title_style="bold")
    table.add_column("block", justify="left", style="cyan")
    table.add_column("shape", justify="right")
    table.add_column("max |analytic − fd|", justify="right")
    for name, block, func, x0 in checks:
        fd = finite_difference_jacobian(func, x0)
        err = float(np.abs(block - fd).max())
        err_style = "green" if err < 1e-5 else "bold red"
        table.add_row(name, str(block.shape), f"[{err_style}]{err:.2e}[/{err_style}]")
    console.print(table)
    console.print(
        Text(f"all blocks computed in {1e3 * elapsed:.2f} ms", style="dim"), justify="right"
    )


if __name__ == "__main__":
    main()

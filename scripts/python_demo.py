"""Solve a small RBC model with TaylorPerturbation at increasing orders and print timings.

Usage: pixi run python_demo [max_order]   (default max_order = 6)
"""

import contextlib
import ctypes
import os
import sys
import tempfile
import time

from rich.console import Console, Group
from rich.panel import Panel
from rich.rule import Rule
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text

from dynare_preprocessor import DynareModel, TaylorPerturbation

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


def main(max_order=6):
    console.print(
        Panel(
            "[bold red]⚠  Warning: the TaylorPerturbation API is [u]not[/u] definitive "
            "and may change.[/bold red]",
            border_style="red",
            expand=False,
        )
    )
    console.print()
    show_code(
        "1. Parse the model",
        """
from dynare_preprocessor import DynareModel, TaylorPerturbation

model = DynareModel(RBC, stochastic=True)   # RBC: .mod source (or a path to a .mod file)
""",
    )
    with capture_fd_stdout() as log:
        t0 = time.perf_counter()
        model = DynareModel(RBC, stochastic=True)
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

    console.print()
    show_code(
        f"2. Solve at orders 1 to {max_order}",
        f"""
for order in range(1, {max_order + 1}):
    tp = TaylorPerturbation(model, order).solve()
    tp.exponents.shape[0]   # number of monomials in w = (xhat, u, sigma)
    tp.max_residual         # max |degree-k residual|
    tp.coef(s=2)            # risk correction (sigma^2 term)
""",
    )
    table = Table(title="TaylorPerturbation", title_style="bold")
    table.add_column("order", justify="right", style="cyan")
    table.add_column("monomials", justify="right")
    table.add_column("time (ms)", justify="right", style="green")
    table.add_column("max residual", justify="right")
    table.add_column("risk corr. c", justify="right", style="magenta")

    ic = model.endogenous.index("c")
    for order in range(1, max_order + 1):
        with capture_fd_stdout():
            t0 = time.perf_counter()
            tp = TaylorPerturbation(model, order).solve()
            elapsed = time.perf_counter() - t0
        risk = tp.coef(s=2)[ic] if order >= 2 else 0.0
        res_style = "green" if tp.max_residual < 1e-10 else "bold red"
        table.add_row(
            str(order),
            str(tp.exponents.shape[0]),
            f"{1e3 * elapsed:.2f}",
            f"[{res_style}]{tp.max_residual:.2e}[/{res_style}]",
            f"{risk:.4e}",
        )
    console.print(table)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

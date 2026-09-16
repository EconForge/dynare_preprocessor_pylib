#!/usr/bin/env python3
"""
Test suite for Dynare Preprocessor exception handling and error reporting.
Verifies line numbers, column numbers/ranges, filenames, and semantic error messages.
"""

import os
import re
import subprocess
import sys
import tempfile
import unittest

# Locate dynare-preprocessor binary
PREPROCESSOR_BIN = None
if len(sys.argv) > 1 and os.path.isfile(sys.argv[1]) and os.access(sys.argv[1], os.X_OK):
    PREPROCESSOR_BIN = os.path.abspath(sys.argv.pop(1))
elif "DYNARE_PREPROCESSOR" in os.environ:
    PREPROCESSOR_BIN = os.path.abspath(os.environ["DYNARE_PREPROCESSOR"])
else:
    candidates = [
        os.path.join(os.getcwd(), "src", "dynare-preprocessor"),
        os.path.join(os.getcwd(), "build", "src", "dynare-preprocessor"),
        os.path.join(os.path.dirname(__file__), "..", "build", "src", "dynare-preprocessor"),
        os.path.join(os.path.dirname(__file__), "src", "dynare-preprocessor"),
    ]
    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            PREPROCESSOR_BIN = os.path.abspath(c)
            break


class PreprocessorExceptionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not PREPROCESSOR_BIN or not os.path.isfile(PREPROCESSOR_BIN):
            raise unittest.SkipTest(
                f"dynare-preprocessor binary not found (searched: {PREPROCESSOR_BIN})"
            )

    def run_mod(self, mod_content, args=None):
        """Helper to create a temporary .mod file and run the preprocessor."""
        if args is None:
            args = []
        with tempfile.TemporaryDirectory() as tmpdir:
            modfile = os.path.join(tmpdir, "model.mod")
            with open(modfile, "w") as f:
                f.write(mod_content)

            proc = subprocess.run(
                [PREPROCESSOR_BIN, modfile] + args,
                cwd=tmpdir,
                capture_output=True,
                text=True,
            )
            combined_output = proc.stdout + proc.stderr
            return proc.returncode, combined_output, modfile

    def parse_error_location(self, output):
        """
        Parses location info from preprocessor output lines formatted like:
        ERROR: <file>: line <line>, col <col>: <msg>
        ERROR: <file>: line <line>, cols <col1>-<col2>: <msg>
        ERROR: <file>: line <line1>, col <col1> - line <line2>, col <col2>: <msg>
        """
        # Pattern 1: Single line, single column
        p1 = re.search(
            r"ERROR:\s*(?:(?P<file>.+?):\s*)?line\s+(?P<line>\d+),\s*col\s+(?P<col>\d+):\s*(?P<msg>.+)",
            output,
        )
        if p1:
            return {
                "file": p1.group("file"),
                "begin_line": int(p1.group("line")),
                "end_line": int(p1.group("line")),
                "begin_col": int(p1.group("col")),
                "end_col": int(p1.group("col")),
                "message": p1.group("msg").strip(),
            }

        # Pattern 2: Single line, column range
        p2 = re.search(
            r"ERROR:\s*(?:(?P<file>.+?):\s*)?line\s+(?P<line>\d+),\s*cols\s+(?P<col1>\d+)-(?P<col2>\d+):\s*(?P<msg>.+)",
            output,
        )
        if p2:
            return {
                "file": p2.group("file"),
                "begin_line": int(p2.group("line")),
                "end_line": int(p2.group("line")),
                "begin_col": int(p2.group("col1")),
                "end_col": int(p2.group("col2")),
                "message": p2.group("msg").strip(),
            }

        # Pattern 3: Multi-line span
        p3 = re.search(
            r"ERROR:\s*(?:(?P<file>.+?):\s*)?line\s+(?P<line1>\d+),\s*col\s+(?P<col1>\d+)\s*-\s*line\s+(?P<line2>\d+),\s*col\s+(?P<col2>\d+):\s*(?P<msg>.+)",
            output,
        )
        if p3:
            return {
                "file": p3.group("file"),
                "begin_line": int(p3.group("line1")),
                "end_line": int(p3.group("line2")),
                "begin_col": int(p3.group("col1")),
                "end_col": int(p3.group("col2")),
                "message": p3.group("msg").strip(),
            }

        return None

    def test_valid_model_succeeds(self):
        """Test that a valid model exits cleanly with status 0."""
        mod = """var y c;
varexo e;
parameters alpha;
alpha = 0.5;
model;
c = y;
y = alpha * c + e;
end;
initval;
y = 0;
c = 0;
end;
steady;
"""
        code, out, _ = self.run_mod(mod)
        self.assertEqual(code, 0, f"Expected 0 exit code, got {code}.\nOutput: {out}")
        self.assertIn("Preprocessing completed.", out)

    def test_syntax_error_single_column(self):
        """Test syntax error on a single token reporting line and exact column."""
        mod = """var y;
model;
y = +;
end;
"""
        code, out, modfile = self.run_mod(mod)
        self.assertEqual(code, 1)

        loc = self.parse_error_location(out)
        self.assertIsNotNone(loc, f"Could not parse location from output:\n{out}")
        self.assertEqual(loc["file"], modfile)
        self.assertEqual(loc["begin_line"], 3)
        self.assertEqual(loc["begin_col"], 6)
        self.assertIn("syntax error", loc["message"].lower())
        self.assertIn(";", loc["message"])

    def test_syntax_error_unexpected_token(self):
        """Test syntax error with unexpected binary operator reporting correct line and column."""
        mod = """var y;
model;
y = 1
    + * 2;
end;
"""
        code, out, _ = self.run_mod(mod)
        self.assertEqual(code, 1)

        loc = self.parse_error_location(out)
        self.assertIsNotNone(loc, f"Could not parse location from output:\n{out}")
        self.assertEqual(loc["begin_line"], 4)
        self.assertEqual(loc["begin_col"], 7)
        self.assertIn("unexpected TIMES", loc["message"])

    def test_undeclared_symbol_column_range(self):
        """Test undeclared model variable reports correct line and column range."""
        mod = """var y;
varexo e;
model;
y = 1 +
    2 +
    undeclared_z;
end;
"""
        code, out, _ = self.run_mod(mod)
        self.assertEqual(code, 1)

        loc = self.parse_error_location(out)
        self.assertIsNotNone(loc, f"Could not parse location from output:\n{out}")
        self.assertEqual(loc["begin_line"], 6)
        # undeclared_z begins at column 5 and ends at column 16
        self.assertEqual(loc["begin_col"], 5)
        self.assertEqual(loc["end_col"], 16)
        self.assertIn("Unknown symbol: undeclared_z", loc["message"])

    def test_duplicate_symbol_multiline_span(self):
        """Test symbol declared twice with different types reports multi-line range."""
        mod = """var y;
parameters
    y;
"""
        code, out, _ = self.run_mod(mod)
        self.assertEqual(code, 1)

        loc = self.parse_error_location(out)
        self.assertIsNotNone(loc, f"Could not parse location from output:\n{out}")
        self.assertEqual(loc["begin_line"], 2)
        self.assertEqual(loc["end_line"], 3)
        self.assertEqual(loc["begin_col"], 1)
        self.assertEqual(loc["end_col"], 6)
        self.assertIn("Symbol y declared twice with different types", loc["message"])

    def test_forbidden_lead_lag_in_planner_objective(self):
        """Test lead/lag in planner_objective reports line and column span."""
        mod = """var y;
planner_objective y(+1);
"""
        code, out, _ = self.run_mod(mod)
        self.assertEqual(code, 1)

        loc = self.parse_error_location(out)
        self.assertIsNotNone(loc, f"Could not parse location from output:\n{out}")
        self.assertEqual(loc["begin_line"], 2)
        self.assertEqual(loc["begin_col"], 19)
        self.assertEqual(loc["end_col"], 23)
        self.assertIn("Leads and lags on variables are forbidden in 'planner_objective'", loc["message"])

    def test_equation_count_mismatch(self):
        """Test semantic error for equation count vs. endogenous variables count mismatch."""
        mod = """var y c;
varexo e;
model;
y = c + e;
end;
"""
        code, out, _ = self.run_mod(mod)
        self.assertEqual(code, 1)
        self.assertIn("ERROR: There are 1 equations but 2 endogenous variables!", out)

    def test_empty_model(self):
        """Test semantic error for empty model when computing task requires equations."""
        mod = """var y;
varexo e;
check;
"""
        code, out, _ = self.run_mod(mod, ["nostrict"])
        self.assertEqual(code, 1)
        self.assertIn("ERROR: At least one model equation must be declared!", out)

    def test_incompatible_statements(self):
        """Test statement conflict error (e.g. discretionary_policy with ramsey_model)."""
        mod = """var y;
varexo e;
model;
y = e;
end;
ramsey_model(planner_discount = 0.99);
discretionary_policy;
"""
        code, out, _ = self.run_mod(mod)
        self.assertEqual(code, 1)
        self.assertTrue(
            "discretionary_policy" in out and ("ramsey_model" in out or "instruments" in out)
        )

    def test_symbol_redeclaration(self):
        """Test semantic error on redeclaring an existing symbol."""
        mod = """var y;
varexo y;
model;
y = 0;
end;
"""
        code, out, modfile = self.run_mod(mod)
        self.assertEqual(code, 1)
        self.assertIn("Symbol y declared twice with different types!", out)
        loc = self.parse_error_location(out)
        self.assertIsNotNone(loc)
        self.assertEqual(loc["file"], modfile)
        self.assertEqual(loc["begin_line"], 2)

    def test_macro_error_directive(self):
        """Test macro processor exception and source location reporting."""
        mod = """@#define foo = 1
@#if foo == 1
@#error "Explicit macro error triggered"
@#endif
var y;
model;
y = 0;
end;
"""
        code, out, modfile = self.run_mod(mod)
        self.assertEqual(code, 1)
        self.assertIn("Explicit macro error triggered", out)
        loc = self.parse_error_location(out)
        self.assertIsNotNone(loc)
        self.assertEqual(loc["file"], modfile)
        self.assertEqual(loc["begin_line"], 3)

    def test_shocks_invalid_variable(self):
        """Test error in shocks block on undeclared/invalid symbol."""
        mod = """var y;
parameters alpha;
model;
y = 0;
end;
shocks;
var alpha; stderr 0.1;
end;
"""
        code, out, _ = self.run_mod(mod)
        self.assertEqual(code, 1)
        self.assertTrue("shocks" in out and ("alpha" in out or "exogenous" in out))

    def test_linear_model_nonlinear_op(self):
        """Test error when model(linear) contains a nonlinear function."""
        mod = """var y;
varexo e;
model(linear);
y = abs(y(-1)) + e;
end;
"""
        code, out, _ = self.run_mod(mod)
        self.assertEqual(code, 1)
        self.assertIn("declared your model 'linear'", out)

    def test_perfect_foresight_stochastic_conflict(self):
        """Test error when mixing perfect foresight solver with stochastic simulation."""
        mod = """var y;
varexo e;
model;
y = e;
end;
perfect_foresight_setup;
perfect_foresight_solver;
stoch_simul;
"""
        code, out, _ = self.run_mod(mod)
        self.assertEqual(code, 1)
        self.assertIn("cannot contain both one of {perfect_foresight_solver", out)


if __name__ == "__main__":
    unittest.main()

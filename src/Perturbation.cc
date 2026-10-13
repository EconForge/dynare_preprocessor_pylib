/*
 * Copyright © 2026 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "Perturbation.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <numbers>
#include <numeric>
#include <unordered_map>

#include <Eigen/Eigenvalues>

#include "DynLib.hh"

namespace perturbation
{

namespace
{

double
seconds_since(std::chrono::steady_clock::time_point t0)
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

} // namespace

/*
 * Basis
 */

Basis::Basis(int nvars_arg, int order_arg, bool products) : nvars {nvars_arg}, order {order_arg}
{
  // Same order as itertools.combinations_with_replacement, degree by degree
  std::vector<int> e(nvars, 0);
  offsets.push_back(0);
  auto gen = [&](auto&& self, int start, int remaining) -> void {
    if (remaining == 0)
      {
        exps.insert(exps.end(), e.begin(), e.end());
        return;
      }
    for (int i = start; i < nvars; i++)
      {
        e[i]++;
        self(self, i, remaining - 1);
        e[i]--;
      }
  };
  for (int d = 0; d <= order; d++)
    {
      gen(gen, 0, d);
      offsets.push_back(static_cast<int>(exps.size() / std::max(nvars, 1)));
    }
  M = offsets.back();
  for (int d = 0; d <= order; d++)
    for (int j = offsets[d]; j < offsets[d + 1]; j++)
      degree.push_back(d);
  for (int j = 0; j < M; j++)
    index_map.emplace(std::vector<int>(exponents(j), exponents(j) + nvars), j);

  parent.assign(M, 0);
  var.assign(M, 0);
  for (int j = 1; j < M; j++)
    {
      std::vector<int> pe(exponents(j), exponents(j) + nvars);
      int v = nvars - 1;
      while (pe[v] == 0)
        v--;
      pe[v]--;
      var[j] = v;
      parent[j] = index(pe);
    }

  const int low = offsets[order];
  lin.resize(static_cast<size_t>(low) * nvars);
  for (int j = 0; j < low; j++)
    {
      std::vector<int> se(exponents(j), exponents(j) + nvars);
      for (int i = 0; i < nvars; i++)
        {
          se[i]++;
          lin[static_cast<size_t>(j) * nvars + i] = index(se);
          se[i]--;
        }
    }

  // Product table, built along the parent tree: e_a + e_b = (e_a + e_parent(b)) + e_var(b)
  if (!products)
    return;
  pair_start.assign(M + 1, 0);
  for (int a = 0; a < M; a++)
    pair_start[a + 1] = pair_start[a] + offsets[order - degree[a] + 1];
  ic.resize(pair_start[M]);
  for (int a = 0; a < M; a++)
    {
      int* row = &ic[pair_start[a]];
      const int nb = offsets[order - degree[a] + 1];
      row[0] = a;
      for (int b = 1; b < nb; b++)
        row[b] = lin[static_cast<size_t>(row[parent[b]]) * nvars + var[b]];
    }
}

int
Basis::index(const std::vector<int>& e) const
{
  auto it = index_map.find(e);
  if (it == index_map.end())
    throw std::out_of_range("Monomial not in basis");
  return it->second;
}

/*
 * Algebra
 */

Algebra::Algebra(const Basis& basis_arg, int k_arg) :
    basis {basis_arg}, k {k_arg}, M {basis_arg.size(k_arg)}
{
}

std::vector<double>
Algebra::constant(double v) const
{
  std::vector<double> c(M, 0.0);
  c[0] = v;
  return c;
}

std::vector<double>
Algebra::variable(int i) const
{
  std::vector<double> c(M, 0.0);
  if (k >= 1)
    c[1 + i] = 1.0;
  return c;
}

std::vector<double>
Algebra::mul(const std::vector<double>& a, const std::vector<double>& b) const
{
  // Iterate over the nonzero terms of the sparser factor
  auto nnz = [this](const std::vector<double>& x) {
    return std::count_if(x.begin(), x.begin() + M, [](double v) { return v != 0.0; });
  };
  const std::vector<double>& outer = nnz(a) <= nnz(b) ? a : b;
  const std::vector<double>& inner = &outer == &a ? b : a;
  std::vector<double> out(M, 0.0);
  for (int ja = 0; ja < M; ja++)
    {
      const double av = outer[ja];
      if (av == 0.0)
        continue;
      const int nb = basis.offsets[k - basis.degree[ja] + 1];
      const int* row = &basis.ic[basis.pair_start[ja]];
      for (int jb = 0; jb < nb; jb++)
        out[row[jb]] += av * inner[jb];
    }
  return out;
}

std::vector<double>
Algebra::shift(const std::vector<double>& a, int i) const
{
  std::vector<double> out(M, 0.0);
  const int nv = basis.nvars;
  for (int j = 0; j < basis.offsets[k]; j++)
    out[basis.lin[static_cast<size_t>(j) * nv + i]] = a[j];
  return out;
}

std::vector<double>
Algebra::mul_linear(const std::vector<double>& a, const double* lf, int deg) const
{
  std::vector<double> out(M, 0.0);
  if (deg >= k)
    return out;
  const int nv = basis.nvars;
  for (int j = basis.offsets[deg]; j < basis.offsets[deg + 1]; j++)
    {
      const double av = a[j];
      if (av == 0.0)
        continue;
      const int* row = &basis.lin[static_cast<size_t>(j) * nv];
      for (int i = 0; i < nv; i++)
        if (lf[i] != 0.0)
          out[row[i]] += av * lf[i];
    }
  return out;
}

std::vector<double>
Algebra::series(const std::vector<double>& x, const std::vector<double>& d) const
{
  // Horner in h = x - x0
  std::vector<double> h(x.begin(), x.begin() + M);
  h[0] = 0.0;
  std::vector<double> r = constant(d[k]);
  for (int j = k - 1; j >= 0; j--)
    {
      r = mul(r, h);
      r[0] += d[j];
    }
  return r;
}

std::vector<double>
Algebra::exp(const std::vector<double>& x) const
{
  std::vector<double> d(k + 1);
  d[0] = std::exp(x[0]);
  for (int j = 1; j <= k; j++)
    d[j] = d[j - 1] / j;
  return series(x, d);
}

std::vector<double>
Algebra::log(const std::vector<double>& x) const
{
  const double x0 = x[0];
  std::vector<double> d(k + 1);
  d[0] = std::log(x0);
  for (int j = 1; j <= k; j++)
    d[j] = (j % 2 == 1 ? 1.0 : -1.0) / (j * std::pow(x0, j));
  return series(x, d);
}

std::vector<double>
Algebra::pow(const std::vector<double>& x, double a) const
{
  if (a == std::floor(a) && a >= 0)
    {
      std::vector<double> r = constant(1.0);
      for (int i = 0; i < static_cast<int>(a); i++)
        r = mul(r, x);
      return r;
    }
  // (x0 + h)^a = x0^a sum_j binom(a, j) (h / x0)^j
  const double x0 = x[0];
  std::vector<double> d(k + 1);
  d[0] = std::pow(x0, a);
  for (int j = 1; j <= k; j++)
    d[j] = d[j - 1] * (a - j + 1) / j / x0;
  return series(x, d);
}

namespace
{

// Taylor coefficients f^(j)(x0) / j! for a function whose derivatives cycle
std::vector<double>
cyclic_coefs(int k, const std::vector<double>& cycle)
{
  std::vector<double> d(k + 1);
  double fact = 1.0;
  for (int j = 0; j <= k; j++)
    {
      if (j > 0)
        fact *= j;
      d[j] = cycle[j % cycle.size()] / fact;
    }
  return d;
}

} // namespace

std::vector<double>
Algebra::sin(const std::vector<double>& x) const
{
  const double s0 = std::sin(x[0]), c0 = std::cos(x[0]);
  return series(x, cyclic_coefs(k, {s0, c0, -s0, -c0}));
}

std::vector<double>
Algebra::cos(const std::vector<double>& x) const
{
  const double s0 = std::sin(x[0]), c0 = std::cos(x[0]);
  return series(x, cyclic_coefs(k, {c0, -s0, -c0, s0}));
}

std::vector<double>
Algebra::sinh(const std::vector<double>& x) const
{
  return series(x, cyclic_coefs(k, {std::sinh(x[0]), std::cosh(x[0])}));
}

std::vector<double>
Algebra::cosh(const std::vector<double>& x) const
{
  return series(x, cyclic_coefs(k, {std::cosh(x[0]), std::sinh(x[0])}));
}

/*
 * Tape
 */

Tape::Tape(const DynareModel& model, std::span<const double> params) :
    n_endo {static_cast<int>(model.endogenous.size())},
    n_exo {static_cast<int>(model.exogenous.size())}
{
  const ModFile& mod_file = model.get_mod_file();
  const SymbolTable& table = mod_file.symbol_table;
  std::unordered_map<expr_t, int> memo;

  auto constant = [&](double v) {
    nodes.push_back({Kind::constant, v});
    return static_cast<int>(nodes.size()) - 1;
  };
  auto is_constant = [&](int i) { return nodes[i].kind == Kind::constant; };

  auto compile = [&](auto&& self, expr_t e) -> int {
    if (auto it = memo.find(e); it != memo.end())
      return it->second;
    int id;
    try
      {
        switch (expression_type(e))
          {
          case ExprNodeType::NumConstNode:
            id = constant(mod_file.num_constants.getDouble(static_cast<NumConstNode*>(e)->id));
            break;
          case ExprNodeType::VariableNode:
            {
              auto* v = static_cast<VariableNode*>(e);
              const int sid = table.getTypeSpecificID(v->symb_id);
              switch (table.getType(v->symb_id))
                {
                case SymbolType::endogenous:
                  if (v->lag < -1 || v->lag > 1)
                    throw ModelSemanticException(
                        "Perturbation: only lags/leads -1, 0 and 1 are supported");
                  nodes.push_back({Kind::argument, 0.0, (1 - v->lag) * n_endo + sid});
                  id = static_cast<int>(nodes.size()) - 1;
                  break;
                case SymbolType::exogenous:
                  if (v->lag != 0)
                    throw ModelSemanticException(
                        "Perturbation requires exogenous variables at date t only");
                  nodes.push_back({Kind::argument, 0.0, 3 * n_endo + sid});
                  id = static_cast<int>(nodes.size()) - 1;
                  break;
                case SymbolType::parameter:
                  id = constant(params[sid]);
                  break;
                default:
                  throw ModelSemanticException(
                      "Perturbation: unsupported variable type in equation ("
                      + table.getName(v->symb_id) + ")");
                }
              break;
            }
          case ExprNodeType::UnaryOpNode:
            {
              auto* u = static_cast<UnaryOpNode*>(e);
              const int a = self(self, u->arg);
              if (u->op_code == UnaryOpcode::steadyState)
                {
                  // STEADY_STATE(x): x evaluated at the expansion point
                  if (is_constant(a))
                    id = a;
                  else
                    {
                      nodes.push_back({Kind::steady_state, 0.0, 0, a});
                      id = static_cast<int>(nodes.size()) - 1;
                      has_steady_state = true;
                    }
                }
              else if (is_constant(a))
                id = constant(UnaryOpNode::eval_opcode(u->op_code, nodes[a].value));
              else
                {
                  nodes.push_back({Kind::unary, 0.0, static_cast<int>(u->op_code), a});
                  id = static_cast<int>(nodes.size()) - 1;
                }
              break;
            }
          case ExprNodeType::BinaryOpNode:
            {
              auto* bn = static_cast<BinaryOpNode*>(e);
              const int a = self(self, bn->arg1);
              const int b = self(self, bn->arg2);
              BinaryOpcode op = bn->op_code;
              if (op == BinaryOpcode::equal)
                op = BinaryOpcode::minus;
              if (is_constant(a) && is_constant(b))
                id = constant(BinaryOpNode::eval_opcode(nodes[a].value, op, nodes[b].value,
                                                        bn->powerDerivOrder));
              else
                {
                  nodes.push_back({Kind::binary, 0.0, static_cast<int>(op), a, b});
                  id = static_cast<int>(nodes.size()) - 1;
                }
              break;
            }
          default:
            throw ModelSemanticException("Perturbation: unsupported expression node");
          }
      }
    catch (const ExprNode::EvalException&)
      {
        throw EvaluationException("Perturbation: evaluation error in constant subexpression");
      }
    memo.emplace(e, id);
    return id;
  };

  for (expr_t eq : model.dynamic_equations())
    outputs.push_back(compile(compile, eq));
}

namespace
{

// Arithmetic on Values. With alg == nullptr all values are scalars.
struct ValueOps
{
  const Algebra* alg;

  [[nodiscard]] std::vector<double>
  as_series(const Value& v) const
  {
    return v.scalar ? alg->constant(v.s) : v.c;
  }
  [[nodiscard]] static Value
  scalar(double s)
  {
    return {true, s, {}};
  }
  [[nodiscard]] static Value
  series(std::vector<double> c)
  {
    return {false, 0.0, std::move(c)};
  }
  [[nodiscard]] static double
  head(const Value& v)
  {
    return v.scalar ? v.s : v.c[0];
  }

  [[nodiscard]] Value
  add(const Value& a, const Value& b, double sign) const
  {
    if (a.scalar && b.scalar)
      return scalar(a.s + sign * b.s);
    if (a.scalar)
      {
        Value r = series(b.c);
        for (double& x : r.c)
          x *= sign;
        r.c[0] += a.s;
        return r;
      }
    Value r = series(a.c);
    if (b.scalar)
      r.c[0] += sign * b.s;
    else
      for (size_t i = 0; i < r.c.size(); i++)
        r.c[i] += sign * b.c[i];
    return r;
  }

  [[nodiscard]] static Value
  scale(const Value& a, double s)
  {
    Value r = a;
    if (r.scalar)
      r.s *= s;
    else
      for (double& x : r.c)
        x *= s;
    return r;
  }

  [[nodiscard]] Value
  times(const Value& a, const Value& b) const
  {
    if (a.scalar)
      return scale(b, a.s);
    if (b.scalar)
      return scale(a, b.s);
    return series(alg->mul(a.c, b.c));
  }

  [[nodiscard]] Value
  divide(const Value& a, const Value& b) const
  {
    if (b.scalar)
      return scale(a, 1.0 / b.s);
    return times(a, series(alg->pow(b.c, -1.0)));
  }

  [[nodiscard]] Value
  power(const Value& a, const Value& b) const
  {
    if (a.scalar && b.scalar)
      return scalar(std::pow(a.s, b.s));
    if (b.scalar)
      return series(alg->pow(a.c, b.s));
    // a^b = exp(b log a)
    Value la = a.scalar ? scalar(std::log(a.s)) : series(alg->log(a.c));
    return series(alg->exp(times(la, b).c));
  }

  [[nodiscard]] Value
  unary(UnaryOpcode op, const Value& a) const
  {
    if (a.scalar)
      return scalar(UnaryOpNode::eval_opcode(op, a.s));
    const auto& x = a.c;
    switch (op)
      {
      case UnaryOpcode::uminus:
        return scale(a, -1.0);
      case UnaryOpcode::exp:
        return series(alg->exp(x));
      case UnaryOpcode::log:
        return series(alg->log(x));
      case UnaryOpcode::log10:
        return scale(series(alg->log(x)), 1.0 / std::numbers::ln10);
      case UnaryOpcode::sqrt:
        return series(alg->pow(x, 0.5));
      case UnaryOpcode::cbrt:
        return series(alg->pow(x, 1.0 / 3.0));
      case UnaryOpcode::sin:
        return series(alg->sin(x));
      case UnaryOpcode::cos:
        return series(alg->cos(x));
      case UnaryOpcode::tan:
        return divide(series(alg->sin(x)), series(alg->cos(x)));
      case UnaryOpcode::sinh:
        return series(alg->sinh(x));
      case UnaryOpcode::cosh:
        return series(alg->cosh(x));
      case UnaryOpcode::tanh:
        return divide(series(alg->sinh(x)), series(alg->cosh(x)));
      case UnaryOpcode::abs:
        return x[0] < 0 ? scale(a, -1.0) : a;
      case UnaryOpcode::sign:
        return scalar(x[0] > 0 ? 1.0 : (x[0] < 0 ? -1.0 : 0.0));
      default:
        throw ModelSemanticException("Perturbation: unary operator (opcode "
                                     + std::to_string(static_cast<int>(op))
                                     + ") not supported on power series");
      }
  }

  [[nodiscard]] Value
  binary(BinaryOpcode op, const Value& a, const Value& b) const
  {
    if (a.scalar && b.scalar)
      return scalar(BinaryOpNode::eval_opcode(a.s, op, b.s, 0));
    switch (op)
      {
      case BinaryOpcode::plus:
        return add(a, b, 1.0);
      case BinaryOpcode::minus:
        return add(a, b, -1.0);
      case BinaryOpcode::times:
        return times(a, b);
      case BinaryOpcode::divide:
        return divide(a, b);
      case BinaryOpcode::power:
        return power(a, b);
      // Piecewise operators: the branch is selected at the expansion point
      case BinaryOpcode::max:
        return head(a) >= head(b) ? a : b;
      case BinaryOpcode::min:
        return head(a) <= head(b) ? a : b;
      case BinaryOpcode::less:
      case BinaryOpcode::greater:
      case BinaryOpcode::lessEqual:
      case BinaryOpcode::greaterEqual:
      case BinaryOpcode::equalEqual:
      case BinaryOpcode::different:
        return scalar(BinaryOpNode::eval_opcode(head(a), op, head(b), 0));
      default:
        throw ModelSemanticException("Perturbation: binary operator (opcode "
                                     + std::to_string(static_cast<int>(op))
                                     + ") not supported on power series");
      }
  }
};

} // namespace

std::vector<Value>
Tape::evaluate(const Algebra* alg, const std::vector<Value>& args) const
{
  const std::vector<Value> vals = run(alg, args, false);
  std::vector<Value> out;
  out.reserve(outputs.size());
  for (int o : outputs)
    out.push_back(vals[o]);
  return out;
}

std::vector<Value>
Tape::run(const Algebra* alg, const std::vector<Value>& args, bool at_steady_state) const
{
  // Values of all nodes at the steady state, for STEADY_STATE(...) operators
  std::vector<Value> ss;
  if (has_steady_state && !at_steady_state)
    {
      if (steady_state.size() != args.size())
        throw EvaluationException("Perturbation: steady state not set for STEADY_STATE operator");
      std::vector<Value> ss_args;
      for (double v : steady_state)
        ss_args.push_back({true, v, {}});
      ss = run(nullptr, ss_args, true);
    }
  ValueOps ops {alg};
  std::vector<Value> vals(nodes.size());
  try
    {
      for (size_t i = 0; i < nodes.size(); i++)
        {
          const Node& nd = nodes[i];
          switch (nd.kind)
            {
            case Kind::constant:
              vals[i] = ValueOps::scalar(nd.value);
              break;
            case Kind::argument:
              vals[i] = args[nd.op];
              break;
            case Kind::steady_state:
              vals[i] = at_steady_state ? vals[nd.a] : ss[nd.a];
              break;
            case Kind::unary:
              vals[i] = ops.unary(static_cast<UnaryOpcode>(nd.op), vals[nd.a]);
              break;
            case Kind::binary:
              vals[i] = ops.binary(static_cast<BinaryOpcode>(nd.op), vals[nd.a], vals[nd.b]);
              break;
            }
        }
    }
  catch (const ExprNode::EvalException&)
    {
      throw EvaluationException("Perturbation: evaluation error in model equations");
    }
  return vals;
}

/*
 * Expectation
 */

namespace
{

// E[eps^alpha] for eps ~ N(0, Sigma) (Isserlis recursion)
double
gaussian_moment(const std::vector<int>& alpha, const Matrix& Sigma,
                std::map<std::vector<int>, double>& memo)
{
  const int total = std::accumulate(alpha.begin(), alpha.end(), 0);
  if (total == 0)
    return 1.0;
  if (total % 2 == 1)
    return 0.0;
  if (auto it = memo.find(alpha); it != memo.end())
    return it->second;
  const int i = static_cast<int>(
      std::distance(alpha.begin(), std::find_if(alpha.begin(), alpha.end(), [](int a) { return a > 0; })));
  std::vector<int> a1 = alpha;
  a1[i]--;
  double s = 0.0;
  for (size_t j = 0; j < a1.size(); j++)
    if (a1[j] > 0 && Sigma(i, j) != 0.0)
      {
        std::vector<int> a2 = a1;
        a2[j]--;
        s += Sigma(i, j) * a1[j] * gaussian_moment(a2, Sigma, memo);
      }
  memo.emplace(alpha, s);
  return s;
}

ExpectationOp
expectation_operator(const Basis& b, int i_sigma, const std::vector<int>& i_up, const Matrix& Sigma)
{
  ExpectationOp E;
  std::map<std::vector<int>, double> memo;
  std::vector<int> alpha(i_up.size());
  for (int j = 0; j < b.M; j++)
    {
      const int* e = b.exponents(j);
      int total = 0;
      for (size_t p = 0; p < i_up.size(); p++)
        total += (alpha[p] = e[i_up[p]]);
      const double mom = gaussian_moment(alpha, Sigma, memo);
      if (mom == 0.0)
        continue;
      std::vector<int> e2(e, e + b.nvars);
      for (int i : i_up)
        e2[i] = 0;
      e2[i_sigma] += total;
      E.src.push_back(j);
      E.tgt.push_back(b.index(e2));
      E.val.push_back(mom);
    }
  return E;
}

} // namespace

RowMatrix
ExpectationOp::apply(const RowMatrix& F) const
{
  const auto Mk = static_cast<int>(F.cols());
  RowMatrix out = RowMatrix::Zero(F.rows(), Mk);
  for (Eigen::Index r = 0; r < F.rows(); r++)
    for (size_t p = 0; p < src.size() && src[p] < Mk; p++)
      out(r, tgt[p]) += val[p] * F(r, src[p]);
  return out;
}

/*
 * Linear solvers
 */

namespace
{

// Reorder a complex Schur form A = U T U^H so that the eigenvalues satisfying `select`
// come first, by swapping adjacent diagonal entries with Givens rotations (as in
// LAPACK's ztrexc). Returns the number of selected eigenvalues.
template<typename Select>
Eigen::Index
reorder_schur(Eigen::MatrixXcd& T, Eigen::MatrixXcd& U, Select select)
{
  const Eigen::Index N = T.rows();
  Eigen::Index ks = 0;
  for (Eigen::Index i = 0; i < N; i++)
    {
      if (!select(T(i, i)))
        continue;
      for (Eigen::Index k = i - 1; k >= ks; k--)
        {
          // Unitary Q whose first column is the eigenvector of the 2x2 block for T(k+1, k+1)
          const std::complex<double> t11 = T(k, k), t22 = T(k + 1, k + 1);
          const std::complex<double> v1 = T(k, k + 1), v2 = t22 - t11;
          const double nv = std::hypot(std::abs(v1), std::abs(v2));
          if (nv == 0.0)
            continue; // equal eigenvalues, nothing to swap
          Eigen::Matrix2cd Q;
          Q << v1 / nv, -std::conj(v2) / nv, v2 / nv, std::conj(v1) / nv;
          T.middleRows(k, 2) = (Q.adjoint() * T.middleRows(k, 2)).eval();
          T.middleCols(k, 2) = (T.middleCols(k, 2) * Q).eval();
          U.middleCols(k, 2) = (U.middleCols(k, 2) * Q).eval();
          T(k + 1, k) = 0.0;
          T(k, k) = t22;
          T(k + 1, k + 1) = t11;
        }
      ks++;
    }
  return ks;
}

} // namespace

Matrix
solve_quadratic(const Matrix& A0, const Matrix& A1, const Matrix& A2)
{
  /* The stable solvent spans the stable invariant subspace of the companion form.
     A2 is typically singular (infinite eigenvalues), so we first substitute
     x = (w - tau) / (1 - tau w), |tau| < 1: multiplied by (1 - tau w)^2, the polynomial
     becomes B0 + B1 w + B2 w^2 with B2 = tau^2 phi(-1/tau), nonsingular for a generic
     tau. This Möbius map preserves the unit disk, so the stable solvent W of the new
     polynomial gives X = (W - tau I)(I - tau W)^{-1}; infinite eigenvalues become 1/tau,
     outside the unit disk. A few Newton steps on the original equation then remove the
     error due to the conditioning of B2. */
  using CMatrix = Eigen::MatrixXcd;
  const Eigen::Index n = A0.rows();
  const Matrix I = Matrix::Identity(n, n);
  auto norm1 = [](const Matrix& A) { return A.cwiseAbs().colwise().sum().maxCoeff(); };
  const double scale = std::max({norm1(A0), norm1(A1), norm1(A2), 1.0});

  // Pick tau giving the best-conditioned B2 (tau = 0 would leave B2 = A2, usually
  // singular). Eigen's rcond() estimate is unreliable for exactly singular matrices,
  // so it is capped by the pivot ratio.
  double tau = 0.0, best = -1.0;
  for (double t : {0.31, -0.37, 0.53, -0.59, 0.71, -0.77})
    {
      const Eigen::PartialPivLU<Matrix> lu(t * t * A0 - t * A1 + A2);
      const auto pivots = lu.matrixLU().diagonal().cwiseAbs();
      double rc = std::min(lu.rcond(), pivots.minCoeff() / pivots.maxCoeff());
      if (!std::isfinite(rc))
        rc = 0.0;
      if (rc > best)
        {
          best = rc;
          tau = t;
        }
    }
  if (!(best > 1e-14))
    throw EvaluationException("First-order solution: singular matrix polynomial");
  const Matrix B0 = A0 - tau * A1 + tau * tau * A2;
  const Matrix B1 = -2 * tau * A0 + (1 + tau * tau) * A1 - 2 * tau * A2;
  const Matrix B2 = tau * tau * A0 - tau * A1 + A2;

  // Companion matrix of w^2 + B2^{-1} B1 w + B2^{-1} B0
  Eigen::PartialPivLU<Matrix> lu2(B2);
  Matrix L = Matrix::Zero(2 * n, 2 * n);
  L.topRightCorner(n, n) = I;
  L.bottomLeftCorner(n, n) = -lu2.solve(B0);
  L.bottomRightCorner(n, n) = -lu2.solve(B1);
  Eigen::ComplexSchur<CMatrix> schur(L.cast<std::complex<double>>());
  if (schur.info() != Eigen::Success)
    throw EvaluationException("First-order solution: Schur decomposition failed");
  CMatrix T = schur.matrixT(), U = schur.matrixU();
  // Same criterion as Dynare (qz_criterium = 1 + 1e-6): unit roots count as stable
  const Eigen::Index nstable
      = reorder_schur(T, U, [](std::complex<double> w) { return std::abs(w) < 1.0 + 1e-6; });
  if (nstable != n)
    throw EvaluationException("Blanchard-Kahn conditions are not satisfied: "
                              + std::to_string(nstable) + " stable eigenvalues for "
                              + std::to_string(n) + " variables");
  const CMatrix U11 = U.topLeftCorner(n, n), U21 = U.bottomLeftCorner(n, n);
  const Matrix W = U11.transpose().partialPivLu().solve(U21.transpose()).transpose().real();
  Matrix X = Eigen::PartialPivLU<Matrix>(I - tau * W).solve(W - tau * I);

  // Newton steps: (A1 + A2 X) H + A2 H X = -(A0 + A1 X + A2 X^2)
  double res = (A0 + A1 * X + A2 * X * X).cwiseAbs().maxCoeff();
  for (int it = 0; it < 3 && res > 0.0; it++)
    {
      const Matrix H = solve_gen_sylvester(A1 + A2 * X, A2, X, A0 + A1 * X + A2 * X * X);
      const Matrix Xn = X + H;
      const double rn = (A0 + A1 * Xn + A2 * Xn * Xn).cwiseAbs().maxCoeff();
      if (!(rn < res))
        break;
      X = Xn;
      res = rn;
    }
  if (!X.allFinite() || res > 1e-8 * scale)
    throw EvaluationException("First-order solution: inaccurate solution (residual "
                              + std::to_string(res) + ")");
  return X;
}

namespace
{

using CMatrix = Eigen::MatrixXcd;

// Solve A D + B D T = -R for D, T being upper triangular, by back-substitution over the
// columns. The factorizations of A + t B are shared between equal diagonal entries t.
CMatrix
solve_triangular_sylvester(const Matrix& A, const Matrix& B, const CMatrix& T, const CMatrix& R)
{
  const Eigen::Index n = R.rows(), N = R.cols();
  const CMatrix Ac = A.cast<std::complex<double>>(), Bc = B.cast<std::complex<double>>();
  std::map<std::pair<double, double>, Eigen::PartialPivLU<CMatrix>> lus;
  CMatrix D = CMatrix::Zero(n, N);
  for (Eigen::Index j = 0; j < N; j++)
    {
      Eigen::VectorXcd rhs = -R.col(j);
      if (j > 0)
        rhs -= Bc * (D.leftCols(j) * T.col(j).head(j));
      const std::complex<double> t = T(j, j);
      auto it = lus.find({t.real(), t.imag()});
      if (it == lus.end())
        it = lus.emplace(std::pair {t.real(), t.imag()}, (Ac + t * Bc).partialPivLu()).first;
      D.col(j) = it->second.solve(rhs);
    }
  return D;
}

/* Matrix of the linear substitution x -> M x on the homogeneous polynomials of degree j
   in b.nvars variables: row p holds the coefficients of m_p(M x) on the degree-j
   monomials of b (in basis order). Sym(M1 M2) = Sym(M1) Sym(M2), and Sym(T) is upper
   triangular when T is: the image of x^a only involves monomials x^c with c >= a
   componentwise, hence lexicographically. */
CMatrix
sym_power(const Basis& b, const CMatrix& M, int j)
{
  const int nv = b.nvars;
  std::vector<Eigen::VectorXcd> img(b.offsets[j + 1]);
  img[0] = Eigen::VectorXcd::Ones(1);
  for (int p = 1; p < b.offsets[j + 1]; p++)
    {
      const int q = b.parent[p], v = b.var[p], d = b.degree[p];
      Eigen::VectorXcd out = Eigen::VectorXcd::Zero(b.offsets[d + 1] - b.offsets[d]);
      for (int r = b.offsets[d - 1]; r < b.offsets[d]; r++)
        {
          const std::complex<double> c = img[q](r - b.offsets[d - 1]);
          if (c == 0.0)
            continue;
          for (int l = 0; l < nv; l++)
            if (M(v, l) != 0.0)
              out(b.lin[static_cast<size_t>(r) * nv + l] - b.offsets[d]) += c * M(v, l);
        }
      img[p] = std::move(out);
    }
  const int Nj = b.offsets[j + 1] - b.offsets[j];
  CMatrix S(Nj, Nj);
  for (int p = 0; p < Nj; p++)
    S.row(p) = img[b.offsets[j] + p].transpose();
  return S;
}

} // namespace

Matrix
solve_gen_sylvester(const Matrix& A, const Matrix& B, const Matrix& S, const Matrix& R)
{
  Eigen::ComplexSchur<CMatrix> schur(S.cast<std::complex<double>>());
  const CMatrix& Q = schur.matrixU();
  const CMatrix D = solve_triangular_sylvester(A, B, schur.matrixT(),
                                               R.cast<std::complex<double>>() * Q);
  return (D * Q.adjoint()).real();
}

SymmetricSylvester::SymmetricSylvester(const Matrix& G) : nx {static_cast<int>(G.rows())}
{
  Eigen::ComplexSchur<CMatrix> schur(G.cast<std::complex<double>>());
  Q = schur.matrixU();
  T = schur.matrixT();
}

const SymmetricSylvester::Block&
SymmetricSylvester::block(int j)
{
  auto it = blocks.find(j);
  if (it == blocks.end())
    {
      auto basis = std::make_unique<Basis>(nx, j, false);
      Block blk {sym_power(*basis, Q, j), sym_power(*basis, T, j),
                 sym_power(*basis, Q.adjoint(), j), std::move(basis)};
      it = blocks.emplace(j, std::move(blk)).first;
    }
  return it->second;
}

Matrix
SymmetricSylvester::solve(const Matrix& A, const Matrix& B, const Matrix& R, int j)
{
  // S = Sym(Q) Sym(T) Sym(Q^H): with D~ = D Sym(Q), A D~ + B D~ Sym(T) = -R Sym(Q)
  const Block& blk = block(j);
  const CMatrix D = solve_triangular_sylvester(A, B, blk.Ts, R.cast<std::complex<double>>() * blk.P);
  return (D * blk.Pinv).real();
}

namespace
{

/* Solve A D + B D S = -R exploiting the block-triangular structure of S (see the
   Python implementation, dyno.experimental.tpsa.solve_order_block). Columns are
   degree-k monomials x^a u^b sigma^c, jdeg = |a| + |b|; the expectation maps them into
   the block jdeg = |a|, so blocks are solved from jdeg = k down to 0. Within a block,
   only the x-only monomials (X) feed back into themselves, through Sym^j(g_xx): this is
   the Sylvester equation, solved by `solve_x`. The others (Y) need a plain solve. */
Matrix
solve_order_block(
    const Matrix& A, const Matrix& B, const Matrix& S, const Matrix& R,
    const std::vector<int>& jdeg, const std::vector<char>& has_u, const std::vector<char>& skip,
    const std::function<Matrix(int j, const std::vector<int>& X, const Matrix& RX)>& solve_x)
{
  const Eigen::Index n = R.rows(), N = R.cols();
  Matrix D = Matrix::Zero(n, N);
  Eigen::PartialPivLU<Matrix> lu(A);
  const int jmax = jdeg.empty() ? -1 : *std::max_element(jdeg.begin(), jdeg.end());
  auto select = [&](const Matrix& Mat, const std::vector<int>& idx) {
    Matrix out(Mat.rows(), idx.size());
    for (size_t q = 0; q < idx.size(); q++)
      out.col(q) = Mat.col(idx[q]);
    return out;
  };
  for (int j = jmax; j >= 0; j--)
    {
      std::vector<int> X, Y;
      for (Eigen::Index q = 0; q < N; q++)
        if (jdeg[q] == j && !skip[q])
          (has_u[q] ? Y : X).push_back(static_cast<int>(q));
      if (!X.empty())
        {
          const Matrix RX = select(R, X) + B * (D * select(S, X));
          const Matrix DX = solve_x(j, X, RX);
          for (size_t q = 0; q < X.size(); q++)
            D.col(X[q]) = DX.col(q);
        }
      if (!Y.empty())
        {
          const Matrix RY = select(R, Y) + B * (D * select(S, Y));
          const Matrix DY = lu.solve(-RY);
          for (size_t q = 0; q < Y.size(); q++)
            D.col(Y[q]) = DY.col(q);
        }
    }
  return D;
}

std::vector<double>
default_params(const DynareModel& model, std::optional<std::vector<double>>& params)
{
  if (params)
    {
      if (params->size() != model.parameters.size())
        throw ModelSemanticException("Perturbation: wrong number of parameters");
      return *params;
    }
  std::vector<double> p;
  for (const auto& name : model.parameters)
    p.push_back(model.context.contains(name) ? model.context.at(name) : 0.0);
  return p;
}

} // namespace

/*
 * TaylorPerturbation
 */

TaylorPerturbation::TaylorPerturbation(const DynareModel& model, int order,
                                       std::optional<std::vector<double>> ybar_arg,
                                       std::optional<std::vector<double>> params,
                                       std::optional<Matrix> Sigma_arg, bool polish,
                                       bool states_only, const std::string& sylvester) :
    K {order},
    n {static_cast<int>(model.endogenous.size())},
    m {static_cast<int>(model.exogenous.size())},
    tape {model, default_params(model, params)}
{
  if (order < 1)
    throw ModelSemanticException("Perturbation: order must be >= 1");
  if (sylvester != "block" && sylvester != "schur")
    throw ModelSemanticException("Perturbation: unknown sylvester method '" + sylvester + "'");
  block = sylvester == "block";
  if (model.equations.size() != static_cast<size_t>(n))
    throw ModelSemanticException("Perturbation: the number of equations must equal the number of endogenous variables");
  if (m == 0)
    throw ModelSemanticException("Perturbation requires a model with exogenous shocks");

  if (ybar_arg)
    {
      if (ybar_arg->size() != static_cast<size_t>(n))
        throw ModelSemanticException("Perturbation: ybar has the wrong size");
      ybar = *ybar_arg;
    }
  else
    for (const auto& name : model.endogenous)
      ybar.push_back(model.context.contains(name) ? model.context.at(name) : 0.0);

  if (Sigma_arg)
    {
      if (Sigma_arg->rows() != m || Sigma_arg->cols() != m)
        throw ModelSemanticException("Perturbation: Sigma has the wrong shape");
      Sigma = *Sigma_arg;
    }
  else
    {
      Sigma = Matrix::Zero(m, m);
      std::map<std::string, int> pos;
      for (int i = 0; i < m; i++)
        pos[model.exogenous[i]] = i;
      for (const auto& [key, v] : model.covariances)
        if (pos.contains(key.first) && pos.contains(key.second))
          Sigma(pos[key.first], pos[key.second]) = Sigma(pos[key.second], pos[key.first]) = v;
    }

  // Newton steps on the steady state: perturbation accuracy is floored by its residual
  if (polish)
    {
      Basis bs(n, 1);
      Algebra alg(bs, 1);
      for (int it = 0; it < 50; it++)
        {
          set_steady_state();
          std::vector<Value> args(3 * n + m);
          for (int i = 0; i < n; i++)
            {
              std::vector<double> c = alg.variable(i);
              c[0] = ybar[i];
              args[i] = args[n + i] = args[2 * n + i] = {false, 0.0, c};
            }
          for (int i = 0; i < m; i++)
            args[3 * n + i] = {true, 0.0, {}};
          const auto out = tape.evaluate(&alg, args);
          Matrix J = Matrix::Zero(n, n);
          Eigen::VectorXd r(n);
          for (int e = 0; e < n; e++)
            {
              r(e) = out[e].scalar ? out[e].s : out[e].c[0];
              if (!out[e].scalar)
                for (int i = 0; i < n; i++)
                  J(e, i) = out[e].c[1 + i];
            }
          if (!r.allFinite())
            throw EvaluationException(
                "Perturbation: non-finite steady-state residual (pass the steady state as ybar)");
          const Eigen::VectorXd dy = J.partialPivLu().solve(r);
          bool small = true;
          for (int i = 0; i < n; i++)
            {
              ybar[i] -= dy(i);
              small = small && std::abs(dy(i)) <= 1e-15 * std::max(1.0, std::abs(ybar[i]));
            }
          if (small || r.cwiseAbs().maxCoeff() == 0.0)
            break;
        }
    }

  set_steady_state();
  const Matrix J = jacobian_at(ybar);
  f_p = J.leftCols(n);
  f_0 = J.middleCols(n, n);
  f_m = J.middleCols(2 * n, n);
  f_e = J.rightCols(m);

  for (int i = 0; i < n; i++)
    if (!states_only || f_m.col(i).cwiseAbs().maxCoeff() > 0)
      states.push_back(i);
  nx = static_cast<int>(states.size());

  // z = (xhat[nx], u[m], sigma, u'[m]); w = first nx + m + 1 variables
  Nw = nx + m + 1;
  i_sigma = nx + m;
  for (int i = 0; i < m; i++)
    i_up.push_back(nx + m + 1 + i);
  const auto t0 = std::chrono::steady_clock::now();
  basis = std::make_unique<Basis>(Nw + m, K);
  timings[0]["basis"] = seconds_since(t0);
  is_w.assign(basis->M, 1);
  for (int j = 0; j < basis->M; j++)
    for (int i : i_up)
      if (basis->exponents(j)[i] > 0)
        is_w[j] = 0;
  E = expectation_operator(*basis, i_sigma, i_up, Sigma);
  C = RowMatrix::Zero(n, basis->M);
}

void
TaylorPerturbation::set_steady_state()
{
  tape.steady_state.assign(3 * n + m, 0.0);
  for (int i = 0; i < 3 * n; i++)
    tape.steady_state[i] = ybar[i % n];
}

Matrix
TaylorPerturbation::jacobian_at(const std::vector<double>& y) const
{
  const int nJ = 3 * n + m;
  Basis bJ(nJ, 1);
  Algebra alg(bJ, 1);
  std::vector<Value> args(nJ);
  for (int i = 0; i < nJ; i++)
    {
      std::vector<double> c = alg.variable(i);
      c[0] = i < 3 * n ? y[i % n] : 0.0;
      args[i] = {false, 0.0, c};
    }
  const auto out = tape.evaluate(&alg, args);
  Matrix J = Matrix::Zero(n, nJ);
  for (int e = 0; e < n; e++)
    if (!out[e].scalar)
      for (int i = 0; i < nJ; i++)
        J(e, i) = out[e].c[1 + i];
  return J;
}

std::vector<int>
TaylorPerturbation::w_monomials(int k) const
{
  std::vector<int> rows;
  for (int j = 0; j < basis->size(k); j++)
    if (is_w[j])
      rows.push_back(j);
  return rows;
}

RowMatrix
TaylorPerturbation::monomial_table(int k, const std::vector<Inner>& inner,
                                   const std::vector<int>& rows, bool linear) const
{
  const Basis& b = *basis;
  const Algebra alg(b, k);
  const int Mk = alg.M;
  std::vector<int> pos(Mk, -1);
  std::vector<std::vector<double>> full(inner.size());
  if (!linear)
    for (size_t v = 0; v < inner.size(); v++)
      if (inner[v].single < 0)
        full[v].assign(inner[v].c, inner[v].c + Mk);
  RowMatrix T(rows.size(), Mk);
  std::vector<double> prev(Mk);
  for (size_t r = 0; r < rows.size(); r++)
    {
      const int j = rows[r];
      pos[j] = static_cast<int>(r);
      if (j == 0)
        {
          T.row(r).setZero();
          T(r, 0) = 1.0;
          continue;
        }
      const int p = b.parent[j], v = b.var[j], d = b.degree[j];
      Eigen::Map<Eigen::RowVectorXd>(prev.data(), Mk) = T.row(pos[p]);
      std::vector<double> out;
      if (inner[v].single >= 0)
        out = alg.shift(prev, inner[v].single);
      else if (linear)
        out = alg.mul_linear(prev, inner[v].c, d - 1);
      else
        out = alg.mul(prev, full[v]);
      T.row(r) = Eigen::Map<const Eigen::RowVectorXd>(out.data(), Mk);
    }
  return T;
}

RowMatrix
TaylorPerturbation::residual(int k) const
{
  const Basis& b = *basis;
  const Algebra alg(b, k);
  const int Mk = alg.M;
  const RowMatrix Ck = C.leftCols(Mk);

  // g(g_x(w), u', sigma): compose the w-monomials with the inner arguments
  std::vector<Inner> inner;
  for (int s : states)
    inner.push_back({-1, Ck.row(s).data()});
  for (int i : i_up)
    inner.push_back({i, nullptr});
  inner.push_back({i_sigma, nullptr});
  const std::vector<int> rows = w_monomials(k);
  const RowMatrix T = monomial_table(k, inner, rows, false);
  RowMatrix Crows(n, rows.size());
  for (size_t r = 0; r < rows.size(); r++)
    Crows.col(r) = Ck.col(rows[r]);
  const RowMatrix Y1 = Crows * T;

  std::vector<Value> args(3 * n + m);
  for (int i = 0; i < n; i++)
    {
      std::vector<double> c1(Y1.row(i).data(), Y1.row(i).data() + Mk);
      c1[0] += ybar[i];
      args[i] = {false, 0.0, std::move(c1)};
      std::vector<double> c0(Ck.row(i).data(), Ck.row(i).data() + Mk);
      c0[0] += ybar[i];
      args[n + i] = {false, 0.0, std::move(c0)};
      args[2 * n + i] = {true, ybar[i], {}};
    }
  for (int p = 0; p < nx; p++)
    {
      std::vector<double> c = alg.variable(p);
      c[0] = ybar[states[p]];
      args[2 * n + states[p]] = {false, 0.0, std::move(c)};
    }
  for (int i = 0; i < m; i++)
    args[3 * n + i] = {false, 0.0, alg.variable(nx + i)};

  const auto out = tape.evaluate(&alg, args);
  RowMatrix F(n, Mk);
  for (int e = 0; e < n; e++)
    {
      if (out[e].scalar)
        {
          F.row(e).setZero();
          F(e, 0) = out[e].s;
        }
      else
        F.row(e) = Eigen::Map<const Eigen::RowVectorXd>(out[e].c.data(), Mk);
    }
  return E.apply(F);
}

std::pair<Matrix, std::vector<int>>
TaylorPerturbation::linear_operator(int k) const
{
  const Basis& b = *basis;
  const int Mk = b.size(k);
  std::vector<Inner> inner;
  for (int s : states)
    inner.push_back({-1, C.row(s).data() + 1}); // linear part
  for (int i : i_up)
    inner.push_back({i, nullptr});
  inner.push_back({i_sigma, nullptr});
  const std::vector<int> rows = w_monomials(k);
  std::vector<int> cols;
  std::vector<int> colpos(Mk, -1);
  for (int j : rows)
    if (b.degree[j] == k)
      {
        colpos[j] = static_cast<int>(cols.size());
        cols.push_back(j);
      }
  const RowMatrix T = monomial_table(k, inner, rows, true);
  std::vector<int> rowpos(Mk, -1);
  for (size_t r = 0; r < rows.size(); r++)
    rowpos[rows[r]] = static_cast<int>(r);

  // S = E(T)[cols, cols]
  const auto N = static_cast<Eigen::Index>(cols.size());
  Matrix S = Matrix::Zero(N, N);
  for (Eigen::Index p = 0; p < N; p++)
    {
      const auto r = rowpos[cols[p]];
      for (size_t q = 0; q < E.src.size() && E.src[q] < Mk; q++)
        if (colpos[E.tgt[q]] >= 0)
          S(p, colpos[E.tgt[q]]) += E.val[q] * T(r, E.src[q]);
    }
  return {S, cols};
}

void
TaylorPerturbation::solve(bool check)
{
  const Basis& b = *basis;
  auto t0 = std::chrono::steady_clock::now();
  Matrix X = solve_quadratic(f_m, f_0, f_p);
  // Columns of non-state variables are zero (f_m has zero columns there)
  std::vector<char> is_state(n, 0);
  for (int s : states)
    is_state[s] = 1;
  const double xmax = X.cwiseAbs().maxCoeff();
  for (int i = 0; i < n; i++)
    if (!is_state[i])
      {
        if (X.col(i).cwiseAbs().maxCoeff() > 1e-10 * std::max(xmax, 1.0))
          throw EvaluationException("First-order solution depends on non-state variables");
        X.col(i).setZero();
      }
  C.setZero();
  for (int p = 0; p < nx; p++)
    for (int i = 0; i < n; i++)
      C(i, 1 + p) = X(i, states[p]);
  const Matrix A = f_0 + f_p * X;
  const Matrix& B = f_p;
  // Transition of the states, x_t = g_xx x_{t-1} + ...: S_k = Sym^j(g_xx) on x-only monomials
  Matrix g_xx(nx, nx);
  for (int p = 0; p < nx; p++)
    for (int q = 0; q < nx; q++)
      g_xx(p, q) = X(states[p], states[q]);
  SymmetricSylvester sym(g_xx);
  const double basis_time = timings[0]["basis"];
  timings.clear();
  timings[0] = {{"basis", basis_time}, {"first_order", seconds_since(t0)}};
  unknowns.clear();
  max_residual = 0.0;

  for (int k = 1; k <= K; k++)
    {
      t0 = std::chrono::steady_clock::now();
      const RowMatrix R = residual(k);
      timings[k]["residual"] = seconds_since(t0);

      t0 = std::chrono::steady_clock::now();
      const auto [S, cols] = linear_operator(k);
      timings[k]["operator"] = seconds_since(t0);

      t0 = std::chrono::steady_clock::now();
      const auto N = static_cast<Eigen::Index>(cols.size());
      std::vector<char> known(N, 0), has_u(N, 0);
      std::vector<int> jdeg(N);
      Matrix Rc(n, N);
      for (Eigen::Index q = 0; q < N; q++)
        {
          const int* e = b.exponents(cols[q]);
          known[q] = k == 1 && cols[q] < 1 + nx; // g_x at order 1
          jdeg[q] = k - e[i_sigma];
          for (int i = 0; i < m; i++)
            has_u[q] = has_u[q] || e[nx + i] > 0;
          Rc.col(q) = R.col(cols[q]);
        }
      Matrix D;
      if (block)
        {
          auto solve_x = [&](int j, const std::vector<int>& X, const Matrix& RX) {
            // Columns X are the x-only monomials of degree j, permuted to basis order
            const Basis& bj = *sym.block(j).basis;
            if (static_cast<int>(X.size()) != bj.offsets[j + 1] - bj.offsets[j])
              throw EvaluationException("Perturbation: unexpected Sylvester block");
            std::vector<int> perm(X.size());
            Matrix Rb(n, X.size());
            for (size_t p = 0; p < X.size(); p++)
              {
                const int* e = b.exponents(cols[X[p]]);
                perm[p] = bj.index(std::vector<int>(e, e + nx)) - bj.offsets[j];
                Rb.col(perm[p]) = RX.col(p);
              }
            const Matrix Db = sym.solve(A, B, Rb, j);
            Matrix DX(n, X.size());
            for (size_t p = 0; p < X.size(); p++)
              DX.col(p) = Db.col(perm[p]);
            return DX;
          };
          D = solve_order_block(A, B, S, Rc, jdeg, has_u, known, solve_x);
        }
      else
        {
          std::vector<int> u;
          for (Eigen::Index q = 0; q < N; q++)
            if (!known[q])
              u.push_back(static_cast<int>(q));
          Matrix Su(u.size(), u.size()), Ru(n, u.size());
          for (size_t p = 0; p < u.size(); p++)
            {
              Ru.col(p) = Rc.col(u[p]);
              for (size_t q = 0; q < u.size(); q++)
                Su(p, q) = S(u[p], u[q]);
            }
          const Matrix Du = solve_gen_sylvester(A, B, Su, Ru);
          D = Matrix::Zero(n, N);
          for (size_t p = 0; p < u.size(); p++)
            D.col(u[p]) = Du.col(p);
        }
      int nunknown = 0;
      for (Eigen::Index q = 0; q < N; q++)
        if (!known[q])
          {
            C.col(cols[q]) = D.col(q);
            nunknown++;
          }
      timings[k]["sylvester"] = seconds_since(t0);
      unknowns[k] = nunknown;

      if (check)
        {
          const RowMatrix Rk = residual(k);
          const double r
              = Rk.middleCols(b.offsets[k], b.offsets[k + 1] - b.offsets[k]).cwiseAbs().maxCoeff();
          max_residual = std::max(max_residual, r);
        }
    }
  solved = true;
}

RowMatrix
TaylorPerturbation::g(const RowMatrix& w, std::optional<int> order) const
{
  if (w.cols() != Nw)
    throw ModelSemanticException("Perturbation: w must have n_w = nx + m + 1 columns");
  const int k = order.value_or(K);
  if (k < 0 || k > K)
    throw ModelSemanticException("Perturbation: invalid truncation order");
  const std::vector<int> idx = w_monomials(k);
  RowMatrix out(w.rows(), n);
  for (Eigen::Index r = 0; r < w.rows(); r++)
    {
      for (int i = 0; i < n; i++)
        out(r, i) = ybar[i];
      for (int j : idx)
        {
          const int* e = basis->exponents(j);
          double mon = 1.0;
          for (int v = 0; v < Nw; v++)
            for (int p = 0; p < e[v]; p++)
              mon *= w(r, v);
          if (mon != 0.0)
            out.row(r) += mon * C.col(j).transpose();
        }
    }
  return out;
}

} // namespace perturbation

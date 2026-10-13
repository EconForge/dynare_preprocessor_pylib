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

#ifndef PERTURBATION_HH
#define PERTURBATION_HH

/*
 * Arbitrary-order perturbation with truncated power series (TPSA).
 *
 * C++ port of dyno.experimental.tpsa. We solve E_t f(y_{t+1}, y_t, y_{t-1}, e_t) = 0
 * for a decision rule
 *
 *     y_t = ybar + g(w),   w = (xhat_{t-1}, u_t, sigma)
 *
 * where xhat are the deviations of the state variables (endogenous variables that
 * appear with a lag), u_t = e_t and sigma scales future shocks. The residual is
 *
 *     E f(ybar + g(g_x(w), u', sigma), ybar + g(w), ybar + xhat, u),   u' = sigma eps'
 *
 * with eps' ~ N(0, Sigma). All computations happen on polynomials in
 * z = (xhat, u, sigma, u') truncated at total degree k; the expectation maps
 * u'^a -> E[eps^a] sigma^|a|. The model equations are never differentiated: their
 * expression trees are evaluated on truncated power series.
 *
 * Order by order: the degree-k part R_k of the residual (with g truncated at k-1)
 * determines g_k through A D + B D S_k = -R_k, with A = f_{y0} + f_{y+} g_x,
 * B = f_{y+} and S_k the action of the expectation of the composition on degree-k
 * monomials. S_k is block-triangular: only the state-only monomials need a
 * generalized Sylvester solve. Order 1 (g_x) solves the quadratic matrix equation.
 */

#include <map>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include <Eigen/Dense>

class DynareModel;

namespace perturbation
{

using Matrix = Eigen::MatrixXd;
using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

//! Monomials in `nvars` variables up to total degree `order`, sorted by degree, so
//! that the coefficients of a series truncated at order k <= order are a prefix.
class Basis
{
public:
  //! Without `products`, the product table (pair_start, ic) is not built
  Basis(int nvars, int order, bool products = true);

  int nvars, order;
  int M; //!< Number of monomials
  std::vector<int> exps; //!< M x nvars exponents (row-major)
  std::vector<int> degree; //!< M
  std::vector<int> offsets; //!< order + 2: first index of each degree, then M
  //! e_j = e_{parent[j]} + e_{var[j]} (j >= 1), var[j] being the last variable of j
  std::vector<int> parent, var;
  //! lin[j * nvars + i]: index of e_j + e_i, for deg j < order
  std::vector<int> lin;
  //! Product table: ic[pair_start[a] + b] is the index of e_a + e_b, for all b with
  //! deg b <= order - deg a (these b form a prefix of the basis)
  std::vector<long> pair_start;
  std::vector<int> ic;

  [[nodiscard]] int size(int k) const { return offsets[k + 1]; }
  [[nodiscard]] const int* exponents(int j) const { return &exps[static_cast<size_t>(j) * nvars]; }
  //! Index of a monomial; throws if not in the basis
  [[nodiscard]] int index(const std::vector<int>& e) const;

private:
  std::map<std::vector<int>, int> index_map;
};

//! A scalar or a truncated power series
struct Value
{
  bool scalar {true};
  double s {0.0};
  std::vector<double> c;
};

//! Truncated power series arithmetic over a basis, at truncation order k
class Algebra
{
public:
  Algebra(const Basis& basis, int k);

  const Basis& basis;
  const int k, M;

  [[nodiscard]] std::vector<double> constant(double v) const;
  [[nodiscard]] std::vector<double> variable(int i) const;
  [[nodiscard]] std::vector<double> mul(const std::vector<double>& a,
                                        const std::vector<double>& b) const;
  //! Product by the single variable z_i
  [[nodiscard]] std::vector<double> shift(const std::vector<double>& a, int i) const;
  //! Product of a homogeneous polynomial of degree deg by the linear form sum lin[i] z_i
  [[nodiscard]] std::vector<double> mul_linear(const std::vector<double>& a, const double* lin,
                                               int deg) const;
  //! sum_j d[j] (x - x0)^j, j = 0..k
  [[nodiscard]] std::vector<double> series(const std::vector<double>& x,
                                           const std::vector<double>& d) const;

  [[nodiscard]] std::vector<double> exp(const std::vector<double>& x) const;
  [[nodiscard]] std::vector<double> log(const std::vector<double>& x) const;
  [[nodiscard]] std::vector<double> pow(const std::vector<double>& x, double a) const;
  [[nodiscard]] std::vector<double> sin(const std::vector<double>& x) const;
  [[nodiscard]] std::vector<double> cos(const std::vector<double>& x) const;
  [[nodiscard]] std::vector<double> sinh(const std::vector<double>& x) const;
  [[nodiscard]] std::vector<double> cosh(const std::vector<double>& x) const;
};

//! The dynamic model equations f(y+[n], y0[n], y-[n], e[m]) compiled into a flat list
//! of operations (shared subexpressions are evaluated once, parameters are folded).
class Tape
{
public:
  Tape(const DynareModel& model, std::span<const double> params);

  //! Residuals for the 3n + m arguments; each argument is a scalar or a series of
  //! `alg` (pass nullptr when all arguments are scalars).
  [[nodiscard]] std::vector<Value> evaluate(const Algebra* alg, const std::vector<Value>& args) const;

  int n_endo, n_exo;
  //! Arguments (3n + m) at which STEADY_STATE(...) operators are evaluated
  std::vector<double> steady_state;

private:
  enum class Kind
  {
    constant,
    argument,
    steady_state,
    unary,
    binary
  };
  struct Node
  {
    Kind kind;
    double value {0.0};
    int op {0}, a {-1}, b {-1};
  };
  std::vector<Node> nodes;
  std::vector<int> outputs;
  bool has_steady_state {false};

  [[nodiscard]] std::vector<Value> run(const Algebra* alg, const std::vector<Value>& args,
                                       bool at_steady_state) const;
};

//! Expectation operator: E[u'^a w^b] = E[eps^a] sigma^|a| w^b, as (src, tgt, val)
//! triplets sorted by src
struct ExpectationOp
{
  std::vector<int> src, tgt;
  std::vector<double> val;
  //! Apply to the rows of F (n x Mk)
  [[nodiscard]] RowMatrix apply(const RowMatrix& F) const;
};

//! Stable solvent of A0 + A1 X + A2 X^2 = 0 (complex Schur decomposition of the
//! companion form, ordered by Givens swaps, followed by Newton refinement)
Matrix solve_quadratic(const Matrix& A0, const Matrix& A1, const Matrix& A2);

//! Solve A D + B D S = -R for D, by complex Schur decomposition of S
Matrix solve_gen_sylvester(const Matrix& A, const Matrix& B, const Matrix& S, const Matrix& R);

/*! Solver for A D + B D Sym^j(G) = -R, Sym^j(G) being the action of the linear map
    x -> G x on the degree-j monomials of x. From the Schur decomposition G = Q T Q^H,
    Sym^j(G) = Sym^j(Q) Sym^j(T) Sym^j(Q^H) with Sym^j(T) upper triangular, so that no
    Schur decomposition of the (much larger) Sym^j(G) is needed. The columns of D and R
    follow the order of the degree-j monomials of Basis(nx, j). */
class SymmetricSylvester
{
public:
  explicit SymmetricSylvester(const Matrix& G);
  [[nodiscard]] Matrix solve(const Matrix& A, const Matrix& B, const Matrix& R, int j);

  struct Block
  {
    Eigen::MatrixXcd P, Ts, Pinv; //!< Sym^j(Q), Sym^j(T), Sym^j(Q^H)
    std::unique_ptr<Basis> basis; //!< Basis(nx, j)
  };
  const Block& block(int j);

private:
  int nx;
  Eigen::MatrixXcd Q, T;
  std::map<int, Block> blocks;
};

class TaylorPerturbation
{
public:
  TaylorPerturbation(const DynareModel& model, int order,
                     std::optional<std::vector<double>> ybar = std::nullopt,
                     std::optional<std::vector<double>> params = std::nullopt,
                     std::optional<Matrix> Sigma = std::nullopt, bool polish = true,
                     bool states_only = true, const std::string& sylvester = "block");

  //! Solve order by order. With `check`, the degree-k residual is recomputed after
  //! each solve and its maximum stored in `max_residual`.
  void solve(bool check = true);

  //! Coefficients (n x Mk) of E f(...) truncated at order k
  [[nodiscard]] RowMatrix residual(int k) const;
  //! S_k on the degree-k w-monomials, and their indices in the basis
  [[nodiscard]] std::pair<Matrix, std::vector<int>> linear_operator(int k) const;
  //! Decision rule (levels) at the rows of w (N x Nw), truncated at `order`
  [[nodiscard]] RowMatrix g(const RowMatrix& w, std::optional<int> order = std::nullopt) const;
  //! Indices (in the basis) of the w-monomials of degree <= k
  [[nodiscard]] std::vector<int> w_monomials(int k) const;

  int K, n, m, nx, Nw, i_sigma;
  bool block;
  std::vector<int> states, i_up;
  std::vector<double> ybar;
  Matrix Sigma, f_p, f_0, f_m, f_e;
  std::unique_ptr<Basis> basis;
  std::vector<char> is_w; //!< Monomials involving no u'
  RowMatrix C; //!< n x M coefficients of g (deviations from ybar)
  std::map<int, std::map<std::string, double>> timings;
  std::map<int, int> unknowns;
  double max_residual {0.0};
  bool solved {false};

private:
  Tape tape;
  ExpectationOp E;

  struct Inner
  {
    int single {-1}; //!< If >= 0, the inner argument is the variable z_single
    const double* c {nullptr}; //!< Otherwise its coefficients (all of them, or the linear part)
  };
  //! Coefficients of the monomials `rows` (sorted, closed under parent) evaluated at
  //! `inner`; with `linear`, the inner arguments are linear forms
  [[nodiscard]] RowMatrix monomial_table(int k, const std::vector<Inner>& inner,
                                         const std::vector<int>& rows, bool linear) const;
  [[nodiscard]] Matrix jacobian_at(const std::vector<double>& y) const;
  void set_steady_state();
};

} // namespace perturbation

#endif

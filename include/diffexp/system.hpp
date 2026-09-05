#pragma once
#include "diffexp/exact.hpp"
#include "diffexp/kernel/recurrence.hpp"
#include <map>

namespace diffexp {
// A matrix is cleared once per system. Charts specialize the immutable exact
// operator; epsilon/Taylor requests never trigger global denominator algebra.
class RationalSystem {
 public:
  using Matrix=std::vector<std::vector<Exact>>;
  RationalSystem(Matrix matrix,std::size_t variable,std::size_t epsilon)
      : matrix_(std::move(matrix)), x_(variable), eps_(epsilon),
        q_(first().constant(1)) {
    if (x_==eps_ || x_>=first().variable_count() || eps_>=first().variable_count())
      throw std::invalid_argument("invalid system variables");
    for (const auto& row:matrix_) {
      if (row.size()!=dimension()) throw std::invalid_argument("system matrix must be square");
      for (const auto& a:row) q_=q_.polynomial_lcm(a.denominator());
    }
    cleared_=matrix_;
    for (auto& row:cleared_) for (auto& a:row) a=a*q_;
  }
  std::size_t dimension() const { return matrix_.size(); }
  const Exact& denominator() const { return q_; }
  const Matrix& cleared_matrix() const { return cleared_; }

  kernel::RecurrenceProblem<Rational> regular_problem(
      const Rational& center,unsigned order,int low,int high,
      const std::vector<std::vector<Rational>>& boundary) const {
    if (high<low || boundary.size()!=dimension() || order>1000000 ||
        static_cast<std::int64_t>(high)-low+1>1000000)
      throw std::invalid_argument("invalid local series request");
    const auto width=static_cast<unsigned>(static_cast<std::int64_t>(high)-low+1);
    for (const auto& row:boundary)
      if (row.size()!=width) throw std::invalid_argument("boundary epsilon width mismatch");
    auto substitutions=identity_substitution();
    substitutions[x_]=substitutions[x_]+q_.constant(center);
    auto q=q_.substitute(substitutions);
    auto origin=identity_substitution();
    origin[x_]=q_.constant(0); origin[eps_]=q_.constant(0);
    const auto leading=q.substitute(origin);
    if (!leading.is_rational()) throw std::invalid_argument("system has unspecialized parameters");
    if (leading.is_zero())
      throw std::domain_error("chart is singular or not epsilon regular; a Frobenius compiler is required");
    kernel::RecurrenceProblem<Rational> p;
    p.dimension=dimension(); p.nmax=order; p.frame_base=low; p.frame_width=width;
    p.log_max=0; p.a_target=Rational(0); p.b_target=Rational(0); p.a_shift_min=0;
    auto qs=coefficients(q);
    for (const auto& [powers,c]:qs) {
      auto [n,e]=powers;
      if (p.d_lags.size()<=n) p.d_lags.resize(n+1);
      p.d_lags[n].push_back({e,c});
    }
    std::map<unsigned,std::map<int,std::vector<kernel::MatrixEntry<Rational>>>> lags;
    for (unsigned i=0;i<dimension();++i)
      for (unsigned j=0;j<dimension();++j)
        for (const auto& [powers,c]:coefficients(cleared_[i][j].substitute(substitutions))) {
          auto [n,e]=powers;
          lags[n+1][e].push_back({i,j,c}); // theta=t d/dt, so Nhat=t P.
        }
    p.nhat_lags.resize(lags.empty()?1:lags.rbegin()->first+1);
    for (auto& lag:p.nhat_lags) lag.valuations.assign(dimension()*dimension(),kernel::kCompleteInfinity);
    for (const auto& [n,shifts]:lags)
      for (const auto& [e,entries]:shifts) {
        p.nhat_lags[n].polynomial.push_back({e,entries});
        for (const auto& entry:entries) {
          auto& v=p.nhat_lags[n].valuations[entry.row*dimension()+entry.col]; v=std::min(v,e);
        }
      }
    if (p.d_lags[0].size()==1 && p.d_lags[0][0].shift==0)
      p.d0_inverse_scalar=Rational(1)/p.d_lags[0][0].value;
    for (unsigned i=0;i<dimension();++i) p.blocks.push_back({{i}});
    p.schedule.resize(order+1);
    for (unsigned n=0;n<=order;++n) {
      p.a_shifts.emplace_back(n);
      for (unsigned i=0;i<dimension();++i)
        p.schedule[n].push_back({n?kernel::StepCase::Taylor:kernel::StepCase::Resonant,Rational(n),Rational(0)});
    }
    for (const auto& row:boundary) p.initial.insert(p.initial.end(),row.begin(),row.end());
    p.initial_validity.assign(dimension(),high);
    return p;
  }
 private:
  Matrix matrix_,cleared_;
  std::size_t x_,eps_;
  Exact q_;
  const Exact& first() const {
    if (matrix_.empty() || matrix_[0].empty()) throw std::invalid_argument("empty system matrix");
    return matrix_[0][0];
  }
  std::vector<Exact> identity_substitution() const {
    std::vector<Exact> vars;
    for (std::size_t i=0;i<q_.variable_count();++i) vars.push_back(q_.variable(i));
    return vars;
  }
  using Coefficients=std::map<std::pair<unsigned,int>,Rational>;
  Coefficients coefficients(const Exact& p) const {
    if (!p.denominator().is_rational()) throw std::logic_error("operator was not polynomially cleared");
    const auto divisor=p.denominator().rational();
    Coefficients out;
    for (const auto& t:p.numerator_terms()) {
      for (std::size_t i=0;i<t.powers.size();++i)
        if (i!=x_ && i!=eps_ && t.powers[i])
          throw std::invalid_argument("system has unspecialized parameters");
      if (t.powers[x_]>1000000 || t.powers[eps_]>1000000)
        throw std::overflow_error("operator polynomial degree exceeds compilation limit");
      out.emplace(std::pair<unsigned,int>{t.powers[x_],t.powers[eps_]},t.coefficient/divisor);
    }
    return out;
  }
};

struct TaylorSeries {
  unsigned order,dimension;
  int low,high;
  std::vector<Rational> coefficients; // n, component, epsilon
  const Rational& at(unsigned n,unsigned component,int epsilon) const {
    if (n>order || component>=dimension || epsilon<low || epsilon>high)
      throw std::out_of_range("series coefficient outside complete window");
    return coefficients.at((static_cast<std::size_t>(n)*dimension+component)*(high-low+1)+epsilon-low);
  }
  // Explicitly a polynomial evaluation, not an infinite-series enclosure.
  std::vector<std::vector<Rational>> evaluate_polynomial(const Rational& t) const {
    std::vector<std::vector<Rational>> out(dimension,std::vector<Rational>(high-low+1,Rational(0)));
    for (unsigned i=0;i<dimension;++i) for (int e=low;e<=high;++e)
      for (unsigned n=order+1;n-->0;) out[i][e-low]=out[i][e-low]*t+at(n,i,e);
    return out;
  }
};

inline TaylorSeries regular_series(const RationalSystem& system,const Rational& center,
    unsigned order,int low,int high,const std::vector<std::vector<Rational>>& boundary) {
  auto p=system.regular_problem(center,order,low,high,boundary);
  auto r=kernel::RecurrenceSolver<Rational>(p).run();
  if (r.top_valid<high) throw std::runtime_error("recurrence returned incomplete requested window");
  TaylorSeries out{order,p.dimension,low,high,{}};
  const auto logs=p.log_max+2;
  for (unsigned n=0;n<=order;++n) for (unsigned i=0;i<p.dimension;++i)
    for (unsigned e=0;e<p.frame_width;++e)
      out.coefficients.push_back(r.u.at(((static_cast<std::size_t>(n)*logs)*p.dimension+i)*p.frame_width+e));
  return out;
}
}  // namespace diffexp

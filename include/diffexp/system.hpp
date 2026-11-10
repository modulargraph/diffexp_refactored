#pragma once
#include "diffexp/exact.hpp"
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

  // Exact Taylor coefficients about a regular centre. With Q the cleared
  // denominator and P=Q*A, both polynomial in x and epsilon, Q y' = P y gives
  //   q00 (n+1) y[n+1][e] = sum P[a][e1] y[n-a][e-e1]
  //                         - sum_{(a,e1)!=(0,0)} q[a][e1] (n+1-a) y[n+1-a][e-e1],
  // which is causal in both the Taylor order and the epsilon power, so the
  // window low..high is exact for every retained coefficient.
  std::vector<Rational> regular_coefficients(const Rational& center,unsigned order,int low,int high,
      const std::vector<std::vector<Rational>>& boundary) const {
    if (high<low || boundary.size()!=dimension() || order>1000000 ||
        static_cast<std::int64_t>(high)-low+1>1000000)
      throw std::invalid_argument("invalid local series request");
    const auto width=static_cast<unsigned>(static_cast<std::int64_t>(high)-low+1);
    for (const auto& row:boundary)
      if (row.size()!=width) throw std::invalid_argument("boundary epsilon width mismatch");
    auto substitutions=identity_substitution();
    substitutions[x_]=substitutions[x_]+q_.constant(center);
    const auto q=q_.substitute(substitutions);
    auto origin=identity_substitution();
    origin[x_]=q_.constant(0); origin[eps_]=q_.constant(0);
    const auto leading=q.substitute(origin);
    if (!leading.is_rational()) throw std::invalid_argument("system has unspecialized parameters");
    if (leading.is_zero())
      throw std::domain_error("chart is singular or not epsilon regular; a Frobenius compiler is required");
    const auto q00=leading.rational();
    const auto denominator=coefficients(q);
    struct Entry {unsigned row,column;Coefficients terms;};
    std::vector<Entry> entries;
    for (unsigned i=0;i<dimension();++i)
      for (unsigned j=0;j<dimension();++j) {
        auto terms=coefficients(cleared_[i][j].substitute(substitutions));
        if (!terms.empty()) entries.push_back({i,j,std::move(terms)});
      }
    const unsigned d=dimension();
    std::vector<Rational> y(static_cast<std::size_t>(order+1)*d*width,Rational(0));
    const auto at=[&](unsigned n,unsigned i,unsigned e)->Rational& {return y[(static_cast<std::size_t>(n)*d+i)*width+e];};
    for (unsigned i=0;i<d;++i) for (unsigned e=0;e<width;++e) at(0,i,e)=boundary[i][e];
    for (unsigned n=0;n<order;++n)
      for (unsigned e=0;e<width;++e) {
        std::vector<Rational> next(d,Rational(0));
        for (const auto& entry:entries)
          for (const auto& [powers,c]:entry.terms) {
            const auto [a,shift]=powers;
            if (a<=n && static_cast<unsigned>(shift)<=e) next[entry.row]+=c*at(n-a,entry.column,e-shift);
          }
        for (const auto& [powers,c]:denominator) {
          const auto [a,shift]=powers;
          if ((a==0 && shift==0) || a>n+1 || static_cast<unsigned>(shift)>e) continue;
          const auto weight=c*Rational(static_cast<long>(n+1-a));
          for (unsigned i=0;i<d;++i) next[i]-=weight*at(n+1-a,i,e-shift);
        }
        const auto pivot=q00*Rational(static_cast<long>(n+1));
        for (unsigned i=0;i<d;++i) at(n+1,i,e)=next[i]/pivot;
      }
    return y;
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
  return {order,static_cast<unsigned>(system.dimension()),low,high,
    system.regular_coefficients(center,order,low,high,boundary)};
}
}  // namespace diffexp

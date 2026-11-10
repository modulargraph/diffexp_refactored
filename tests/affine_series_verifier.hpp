#pragma once
// Independent exact verifier of a retained affine Frobenius series, used by
// tests as an oracle for AffineFrobeniusSeries::prepare. It does not share the
// recurrence: it substitutes the retained coefficients into the polynomial
// form of the differential equation, and checks the canonical frame, the
// common frontier, the resonant integration constants and the Wronskian.
#include "diffexp/affine_frobenius.hpp"
#include "diffexp/geometry.hpp"
#include <flint/fmpq_poly.h>

namespace diffexp::dense_polynomial_residual {
// Exact Q[epsilon] arithmetic for checking polynomial differential identities.
// This representation is independent of the Frobenius construction recurrence.
class Polynomial {
  fmpq_poly_t value_;
 public:
  Polynomial(){fmpq_poly_init(value_);}
  Polynomial(const Polynomial& p):Polynomial(){fmpq_poly_set(value_,p.value_);}
  Polynomial(Polynomial&& p)noexcept:Polynomial(){fmpq_poly_swap(value_,p.value_);}
  Polynomial& operator=(Polynomial p)noexcept{fmpq_poly_swap(value_,p.value_);return *this;}
  ~Polynomial(){fmpq_poly_clear(value_);}
  Polynomial(const Exact& p,std::size_t variable):Polynomial(){
    fmpz_poly_t n,d;fmpz_poly_init(n);fmpz_poly_init(d);
    try{p.univariate_polynomials(n,d,variable);if(fmpz_poly_degree(d)!=0)throw std::invalid_argument("residual coefficient is not polynomial");
      fmpq_poly_set_fmpz_poly(value_,n);fmpq_poly_scalar_div_fmpz(value_,value_,d->coeffs);
    }catch(...){fmpz_poly_clear(n);fmpz_poly_clear(d);throw;}
    fmpz_poly_clear(n);fmpz_poly_clear(d);
  }
  Exact exact(const Exact& sample,std::size_t variable)const{
    fmpz_poly_t n,d;fmpz_poly_init(n);fmpz_poly_init(d);
    fmpq_poly_get_numerator(n,value_);fmpz_poly_set_fmpz(d,value_->den);
    try{auto result=sample.from_univariate_polynomials(n,d,variable);fmpz_poly_clear(n);fmpz_poly_clear(d);return result;}
    catch(...){fmpz_poly_clear(n);fmpz_poly_clear(d);throw;}
  }
  bool zero()const{return fmpq_poly_is_zero(value_);}
  std::size_t terms()const{std::size_t n=0;for(slong i=0;i<fmpq_poly_length(value_);++i)n+=!fmpz_is_zero(value_->coeffs+i);return n;}
  void set(unsigned exponent,const Rational& coefficient){fmpq_poly_set_coeff_fmpq(value_,exponent,coefficient.raw());}
  Polynomial& operator+=(const Polynomial& p){fmpq_poly_add(value_,value_,p.value_);return *this;}
  Polynomial& operator-=(const Polynomial& p){fmpq_poly_sub(value_,value_,p.value_);return *this;}
  Polynomial scaled(long n)const{Polynomial p;fmpq_poly_scalar_mul_si(p.value_,value_,n);return p;}
  friend Polynomial operator*(const Polynomial& a,const Polynomial& b){Polynomial p;fmpq_poly_mul(p.value_,a.value_,b.value_);return p;}
};
using Bivariate=std::vector<Polynomial>; // coefficients in x, low order first
inline Bivariate split(const Exact& value,std::size_t xi,std::size_t ei,unsigned cutoff){
  if(!value.denominator().is_rational())throw std::invalid_argument("residual expected polynomial");
  const auto denominator=value.denominator().rational();Bivariate out;
  for(const auto& term:value.numerator_terms()){
    for(std::size_t v=0;v<term.powers.size();++v)if(v!=xi&&v!=ei&&term.powers[v])throw std::invalid_argument("unassigned residual polynomial variable");
    if(term.powers[xi]>cutoff)continue;
    if(out.size()<=term.powers[xi])out.resize(term.powers[xi]+1);
    out[term.powers[xi]].set(term.powers[ei],term.coefficient/denominator);
  }
  return out;
}
template<class Check> inline Bivariate multiply(const Bivariate& a,const Bivariate& b,unsigned cutoff,Check&& check){
  if(a.empty()||b.empty())return {};
  Bivariate out(std::min<std::size_t>(cutoff+1,a.size()+b.size()-1));
  for(std::size_t i=0;i<a.size()&&i<out.size();++i)if(!a[i].zero())
    for(std::size_t j=0;j<b.size()&&i+j<out.size();++j)if(!b[j].zero()){
      check(a[i].terms(),b[j].terms());out[i+j]+=a[i]*b[j];
    }
  return out;
}
} // namespace diffexp::dense_polynomial_residual

namespace diffexp::affine_series_verifier {
using Series=AffineFrobeniusSeries;
using Matrix=Series::Matrix;
struct VerificationLimits {
  std::size_t max_polynomial_terms=1000000, max_shared_terms=4000000, max_total_column_terms=4000000,
              max_term_products=2000000000;
  // Residuals over dense Q[eps] coefficients per x order; false keeps them as
  // sparse bivariate polynomials. Same checks, about half the work.
  bool dense_epsilon_polynomials=true;
  // Cumulative conservative monomial-product estimate, not elapsed CPU work.
  std::function<void(unsigned,unsigned,std::size_t)> column_progress;
};
inline void require(bool condition,const char* why) {
  if(!condition)throw std::runtime_error(why);
}
// The series contradicts the exact mathematics. Receiving limits keep the
// plain runtime_error and are never mistaken for a wrong series.
struct Mismatch : std::runtime_error {
  using std::runtime_error::runtime_error;
};
inline void matches(bool condition,const char* why) {
  if(!condition)throw Mismatch(why);
}
inline void input(const Matrix& a,std::size_t xi,std::size_t ei,unsigned n,const Series::Options& options) {
  require(!a.empty() && a.size()<=256 && a.size()<=options.max_dimension &&
      n<=1000 && n<=options.max_x_order && options.max_terms && options.max_epsilon_depth,
      "affine series input exceeds receiving dimensions/order/term budgets");
  require(!a[0].empty(),"affine series empty source row");
  const auto& sample=a[0][0];
  require(xi!=ei && xi<sample.variable_count() && ei<sample.variable_count(),"affine series variable indices");
  for(const auto& row:a) {
    require(row.size()==a.size(),"affine series source matrix must be square");
    for(const auto& c:row) {
      (void)(sample.constant(0)+c); // Also rejects distinct exact field instances.
      for(const auto& terms:{c.numerator_terms(),c.denominator_terms()})
        for(const auto& term:terms)for(unsigned j=0;j<term.powers.size();++j)
          require(!term.powers[j] || j==xi || j==ei,"affine series source has an unsubstituted parameter");
    }
  }
}
struct State {
  Matrix frame;
  std::vector<Series::Exponent> exponents;
  std::vector<int> successors;
  std::vector<Rational> frontiers;
  Series::Expansion expansion;
};
inline State state_of(const Series& series) {
  return {series.residue_frame(),series.exponents(),series.jordan_successors(),series.absolute_x_frontiers(),series.terms()};
}
struct PolynomialBudget {
  VerificationLimits limits;
  std::size_t products=0;
  std::size_t check(const Exact& p) const {
    require(p.denominator().is_rational(),"affine series residual expected a polynomial");
    auto size=p.numerator_terms().size();require(size<=limits.max_polynomial_terms,"affine series polynomial term budget exhausted");return size;
  }
  Exact multiply(const Exact& a,const Exact& b) {
    auto na=check(a),nb=check(b);
    require(!nb || na<=(limits.max_term_products-products)/nb,"affine series residual work budget exhausted");
    products+=na*nb;auto out=a*b;check(out);return out;
  }
};
inline long integer_offset(const Rational& offset,unsigned n) {
  matches(affine_frobenius_detail::is_integer(offset) && offset.sign()>=0 && offset<=Rational(n),"affine series power outside retained integer frontier");
  return affine_frobenius_detail::integer_value(offset);
}
// Independent verification uses polynomial residuals, not the preparation
// recurrence. A common epsilon denominator removes rational-expression gcds
// from the large residual products. No residual beyond the frontier is claimed.
inline void verify(const Matrix& a,std::size_t xi,std::size_t ei,unsigned n,const Series::Options& options,
    const State& state,const VerificationLimits& limits={}) {
  input(a,xi,ei,n,options);
  require(limits.max_polynomial_terms && limits.max_shared_terms && limits.max_total_column_terms && limits.max_term_products,"affine series zero verification budget");
  PolynomialBudget budget{limits};const unsigned d=a.size();const auto z=a[0][0].constant(0),x=z.variable(xi);
  auto connection=fuchsify::detail::zeros(d,d,z),residue=connection;
  for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j){connection[i][j]=x*a[i][j];residue[i][j]=affine_frobenius_detail::taylor(connection[i][j],xi,0)[0];}
  const auto canonical=affine_frobenius_detail::frame(residue,xi,ei);
  matches(state.frame==canonical.transform && state.successors==canonical.successor && state.exponents.size()==d && state.frontiers.size()==d,"affine series canonical Jordan frame mismatch");
  std::vector<unsigned> cutoffs;
  for(unsigned c=0;c<d;++c) {
    matches(state.exponents[c].power==canonical.exponents[c].power && state.exponents[c].slope==canonical.exponents[c].slope,"affine series canonical spectrum mismatch");
    auto minimum=canonical.exponents[c].power;
    for(unsigned j=0;j<d;++j)if(affine_frobenius_detail::is_integer(canonical.exponents[j].power-canonical.exponents[c].power))minimum=std::min(minimum,canonical.exponents[j].power);
    matches(state.frontiers[c]==minimum+Rational(n),"affine series absolute frontier mismatch");
    cutoffs.push_back(integer_offset(state.frontiers[c]-canonical.exponents[c].power,n));
  }
  matches(state.expansion.rows==d && state.expansion.columns==d && state.expansion.coherent_x_frontier && !state.expansion.omitted_tail_certified,"affine series expansion shape/certificate budget");
  require(state.expansion.terms.size()<=options.max_terms,"affine series expansion shape/certificate budget");
  using Coordinate=std::tuple<unsigned,unsigned,unsigned,unsigned>; // column,k,row,log
  std::map<Coordinate,const Series::Term*> coordinates;
  std::vector<std::vector<const Series::Term*>> columns(d);
  for(const auto& term:state.expansion.terms) {
    matches(term.row<d && term.column<d && term.log_degree<=static_cast<std::size_t>(d)*(n+1) && !term.coefficient.is_zero(),"affine series invalid nonzero term coordinate");
    matches(term.slope==canonical.exponents[term.column].slope,"affine series term slope mismatch");
    auto k=integer_offset(term.power-canonical.exponents[term.column].power,cutoffs[term.column]);
    for(const auto& terms:{term.coefficient.numerator_terms(),term.coefficient.denominator_terms()})
      for(const auto& t:terms)for(unsigned j=0;j<t.powers.size();++j)matches(!t.powers[j] || j==ei,"affine series coefficient depends on a non-epsilon variable");
    matches(coordinates.emplace(Coordinate{term.column,k,term.row,term.log_degree},&term).second,"affine series duplicate term coordinate");
    columns[term.column].push_back(&term);
  }
  // Polynomial ODE data are shared across all retained columns and log powers.
  std::vector<Exact> q(d,z.constant(1));auto p=connection;
  std::size_t shared_terms=0;
  auto count_shared=[&](const Exact& value){auto count=budget.check(value);require(count<=limits.max_shared_terms-shared_terms,"affine series shared polynomial storage budget exhausted");shared_terms+=count;};
  for(unsigned i=0;i<d;++i) {
    for(unsigned j=0;j<d;++j){q[i]=q[i].polynomial_lcm(connection[i][j].denominator());budget.check(q[i]);}
    count_shared(q[i]);
    for(unsigned j=0;j<d;++j){p[i][j]=q[i]*connection[i][j];count_shared(p[i][j]);}
  }
  std::vector<dense_polynomial_residual::Bivariate> dense_q;
  std::vector<std::vector<dense_polynomial_residual::Bivariate>> dense_p;
  if(limits.dense_epsilon_polynomials){
    for(unsigned i=0;i<d;++i){dense_q.push_back(dense_polynomial_residual::split(q[i],xi,ei,n));dense_p.emplace_back();for(unsigned j=0;j<d;++j)dense_p.back().push_back(dense_polynomial_residual::split(p[i][j],xi,ei,n));}
  }
  const auto inverse=fuchsify::detail::inverse(canonical.transform);
  std::vector<Exact> powers(n+1,z.constant(1));for(unsigned k=1;k<=n;++k)powers[k]=powers[k-1]*x;
  for(unsigned c=0;c<d;++c) {
    auto common=z.constant(1);unsigned log_high=0;
    // Denominators repeat across physical rows and log powers. Intern their
    // LCM contributions and exact quotients once per column.
    std::map<std::string,Exact> denominators;
    std::vector<std::pair<const Series::Term*,std::string>> indexed;
    for(const auto* t:columns[c]) {
      auto denominator=t->coefficient.denominator();auto name=denominator.str();
      denominators.try_emplace(name,std::move(denominator));indexed.emplace_back(t,std::move(name));
      log_high=std::max(log_high,t->log_degree);
    }
    for(const auto& [name,denominator]:denominators){common=common.polynomial_lcm(denominator);budget.check(common);}
    for(auto& [name,denominator]:denominators)denominator=common/denominator;
    require(d<=limits.max_total_column_terms/(static_cast<std::size_t>(log_high)+1),"affine series column polynomial slot budget exhausted");
    if(limits.dense_epsilon_polynomials) {
      namespace dp=dense_polynomial_residual;
      auto spend=[&](std::size_t na,std::size_t nb){
        require(na<=limits.max_polynomial_terms && nb<=limits.max_polynomial_terms,"affine series polynomial term budget exhausted");
        require(!nb || na<=(limits.max_term_products-budget.products)/nb,"affine series residual work budget exhausted");budget.products+=na*nb;
      };
      auto check=[&](const dp::Bivariate& value){std::size_t count=0;for(const auto& p:value){auto terms=p.terms();require(terms<=limits.max_polynomial_terms-count,"affine series polynomial term budget exhausted");count+=terms;}return count;};
      const auto cutoff=cutoffs[c];
      std::map<std::string,dp::Polynomial> quotients;
      for(const auto& [name,denominator]:denominators)quotients.emplace(name,dp::Polynomial(denominator,ei));
      std::vector<std::vector<dp::Bivariate>> numerator(d,std::vector<dp::Bivariate>(log_high+1));
      for(const auto& [t,name]:indexed) {
        const auto k=integer_offset(t->power-canonical.exponents[c].power,cutoff);
        dp::Polynomial coefficient(t->coefficient.numerator(),ei);const auto& quotient=quotients.at(name);spend(coefficient.terms(),quotient.terms());
        auto& dest=numerator[t->row][t->log_degree];if(dest.size()<=static_cast<std::size_t>(k))dest.resize(k+1);
        dest[k]+=coefficient*quotient;check(dest);
      }
      std::size_t total=0;for(const auto& row:numerator)for(const auto& value:row){auto terms=check(value);require(terms<=limits.max_total_column_terms-total,"affine series total column polynomial budget exhausted");total+=terms;}
      dp::Polynomial lambda;lambda.set(0,canonical.exponents[c].power);lambda.set(1,canonical.exponents[c].slope);
      for(unsigned i=0;i<d;++i)for(unsigned l=0;l<=log_high;++l) {
        auto derivative=numerator[i][l];
        for(std::size_t k=0;k<derivative.size();++k){spend(lambda.terms(),numerator[i][l][k].terms());derivative[k]=numerator[i][l][k].scaled(k);derivative[k]+=lambda*numerator[i][l][k];}
        if(l<log_high){const auto& next=numerator[i][l+1];if(derivative.size()<next.size())derivative.resize(next.size());for(std::size_t k=0;k<next.size();++k)derivative[k]+=next[k].scaled(l+1);}
        check(derivative);auto residual=dp::multiply(dense_q[i],derivative,cutoff,spend);check(residual);
        for(unsigned j=0;j<d;++j)if(!dense_p[i][j].empty()&&!numerator[j][l].empty()){
          auto term=dp::multiply(dense_p[i][j],numerator[j][l],cutoff,spend);check(term);
          if(residual.size()<term.size())residual.resize(term.size());for(std::size_t k=0;k<term.size();++k)residual[k]-=term[k];check(residual);
        }
        for(const auto& coefficient:residual)matches(coefficient.zero(),"affine series defining polynomial residual is nonzero within frontier");
      }
    } else {
    std::vector<std::vector<Exact>> numerator(d,std::vector<Exact>(log_high+1,z));
    for(const auto& [t,name]:indexed) {
      const auto k=integer_offset(t->power-canonical.exponents[c].power,cutoffs[c]);
      auto regular=budget.multiply(t->coefficient.numerator(),denominators.at(name));
      auto term=budget.multiply(regular,powers[k]);
      numerator[t->row][t->log_degree]=numerator[t->row][t->log_degree]+term;
      budget.check(numerator[t->row][t->log_degree]);
    }
    std::size_t total=0;for(const auto& row:numerator)for(const auto& v:row){auto terms=budget.check(v);require(terms<=limits.max_total_column_terms-total,"affine series total column polynomial budget exhausted");total+=terms;}
    const auto lambda=z.constant(canonical.exponents[c].power)+z.constant(canonical.exponents[c].slope)*z.variable(ei);
    for(unsigned i=0;i<d;++i)for(unsigned l=0;l<=log_high;++l) {
      auto derivative=budget.multiply(x,numerator[i][l].derivative(xi))+budget.multiply(lambda,numerator[i][l]);
      if(l<log_high)derivative=derivative+numerator[i][l+1]*z.constant(l+1);
      auto residual=budget.multiply(q[i],derivative);
      for(unsigned j=0;j<d;++j)if(!p[i][j].is_zero() && !numerator[j][l].is_zero()) {
        residual=residual-budget.multiply(p[i][j],numerator[j][l]);budget.check(residual);
      }
      for(const auto& term:residual.numerator_terms())matches(term.powers[xi]>cutoffs[c],"affine series defining polynomial residual is nonzero within frontier");
    }
    }
    // A polynomial ODE alone permits arbitrary resonant homogeneous constants.
    // These checks bind precisely the zero-integration-constant convention used
    // by prepare, including each seeded diagonal constant at k=0.
    for(unsigned i=0;i<d;++i)if(canonical.exponents[i].slope==canonical.exponents[c].slope) {
      const auto offset=canonical.exponents[i].power-canonical.exponents[c].power;
      if(offset.sign()<0 || offset>Rational(cutoffs[c]) || !affine_frobenius_detail::is_integer(offset))continue;
      const auto k=integer_offset(offset,cutoffs[c]);auto value=z;
      for(unsigned j=0;j<d;++j)if(auto found=coordinates.find(Coordinate{c,k,j,0});found!=coordinates.end())value=value+inverse[i][j]*found->second->coefficient;
      matches(value==z.constant(k==0 && i==c?1:0),"affine series resonant integration constant mismatch");
    }
    if(limits.column_progress)limits.column_progress(c+1,d,budget.products);
  }
  auto regular_trace=z;
  for(unsigned i=0;i<d;++i)regular_trace=regular_trace+a[i][i]-(z.constant(canonical.exponents[i].power)+z.constant(canonical.exponents[i].slope)*z.variable(ei))/x;
  std::vector<Exact> origin;for(unsigned j=0;j<z.variable_count();++j)origin.push_back(j==xi || j==ei?z:z.variable(j));
  std::optional<Exact> wronskian;
  if(!regular_trace.denominator().substitute(origin).is_zero())wronskian=affine_frobenius_detail::determinant(canonical.transform);
  matches(state.expansion.wronskian_prefactor==wronskian,"affine series Wronskian metadata mismatch");
}
inline void verify(const Matrix& a,std::size_t xi,std::size_t ei,const Series& series,const Series::Options& options={},
    const VerificationLimits& limits={}) {
  verify(a,xi,ei,series.x_order(),options,state_of(series),limits);
}
} // namespace diffexp::affine_series_verifier

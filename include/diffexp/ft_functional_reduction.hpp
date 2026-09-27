#pragma once
// Exact algebraic certificates at ordinary rational cuts. These do not certify
// singular endpoint prescriptions, branch choices, analytic tails or ball error.
#include "diffexp/epsilon_gauge.hpp"
#include "diffexp/fuchsify.hpp"
#include <chrono>
#include <cmath>
#include <map>
#include <optional>

namespace diffexp::ft_functional_reduction {
using Matrix = ExactEpsilonMatrix;
struct Functional {
  Matrix connection, lower, weight, upper;
  Exact lower_point, upper_point;
  std::size_t x = 0, epsilon = 1;
};
namespace detail {
inline void shape(const Matrix& a, std::size_t rows, std::size_t columns, const Exact& sample) {
  if (!rows || !columns || rows > 256 || columns > 256 || a.size() != rows)
    throw std::invalid_argument("functional matrix shape/budget");
  for (const auto& row : a) {
    if (row.size() != columns) throw std::invalid_argument("functional matrix column count");
    for (const auto& v : row) sample.require_same_field(v);
  }
}
inline const Exact& sample(const Matrix& a) {
  if (a.empty() || a.front().empty()) throw std::invalid_argument("functional empty matrix");
  return a.front().front();
}
inline void indices(const Exact& z, std::size_t x, std::size_t e) {
  if (x == e || x >= z.variable_count() || e >= z.variable_count())
    throw std::invalid_argument("functional variable indices");
}
inline void constant_in_x(const Matrix& a, std::size_t x) {
  for (const auto& row : a) for (const auto& v : row)
    if (!v.derivative(x).is_zero()) throw std::invalid_argument("endpoint map must be evaluated at its cut");
}
inline Matrix at(const Matrix& a, std::size_t x, const Exact& point) {
  const auto& z = sample(a); z.require_same_field(point);
  if (!point.is_rational()) throw std::invalid_argument("functional requires ordinary rational cut coordinates");
  std::vector<Exact> substitutions;
  for (std::size_t i=0;i<z.variable_count();++i) substitutions.push_back(i==x?point:z.variable(i));
  auto out=a;
  for (auto& row:out) for(auto& v:row) v=v.substitute(substitutions);
  return out;
}
inline Matrix add(Matrix a, const Matrix& b) {
  for(std::size_t i=0;i<a.size();++i)for(std::size_t j=0;j<a[i].size();++j)a[i][j]=a[i][j]+b[i][j];
  return a;
}
inline void equal(const Matrix& a, const Matrix& b, const char* message) {
  if(a!=b)throw std::invalid_argument(message);
}
} // namespace detail

inline bool epsilon_regular(const Matrix& a,std::size_t epsilon) {
  for(const auto& row:a)for(const auto& v:row) {
    auto order=exact_epsilon_valuation(v,epsilon);
    if(order && *order<0)return false;
  }
  return true;
}
inline void validate(const Functional& f) {
  const auto& z=detail::sample(f.connection);
  detail::indices(z,f.x,f.epsilon);
  const auto n=f.connection.size(),m=f.weight.size();
  detail::shape(f.connection,n,n,z);
  for(const auto* a:{&f.lower,&f.weight,&f.upper})detail::shape(*a,m,n,z);
  z.require_same_field(f.lower_point);z.require_same_field(f.upper_point);
  if(!f.lower_point.is_rational() || !f.upper_point.is_rational())
    throw std::invalid_argument("functional requires ordinary rational cut coordinates");
  detail::constant_in_x(f.lower,f.x);detail::constant_in_x(f.upper,f.x);
}
inline bool epsilon_regular(const Functional& f) {
  validate(f);
  return epsilon_regular(f.connection,f.epsilon) && epsilon_regular(f.lower,f.epsilon) &&
    epsilon_regular(f.weight,f.epsilon) && epsilon_regular(f.upper,f.epsilon);
}

// Exact total-derivative identity. Equality of prescribed integrals additionally
// requires C to be admissible along the contour and at its boundaries. This
// algebraic certificate does not authorize new spatial poles or regularizations.
inline Functional rewrite_counterterm(const Functional& f,const Matrix& c) {
  validate(f);detail::shape(c,f.weight.size(),f.connection.size(),detail::sample(f.connection));
  auto out=f;
  out.lower=fuchsify::detail::subtract(f.lower,detail::at(c,f.x,f.lower_point));
  out.weight=fuchsify::detail::subtract(f.weight,detail::add(
    fuchsify::detail::derivative(c,f.x),fuchsify::detail::multiply(c,f.connection)));
  out.upper=detail::add(f.upper,detail::at(c,f.x,f.upper_point));
  return out;
}
struct CountertermCertificate { Matrix counterterm; Functional rewritten; };
inline void verify_counterterm(const Functional& original,const CountertermCertificate& certificate) {
  auto expected=rewrite_counterterm(original,certificate.counterterm);
  const auto& actual=certificate.rewritten;validate(actual);
  if(actual.x!=expected.x || actual.epsilon!=expected.epsilon ||
     actual.lower_point!=expected.lower_point || actual.upper_point!=expected.upper_point)
    throw std::invalid_argument("counterterm certificate metadata mismatch");
  detail::equal(actual.connection,expected.connection,"counterterm connection changed");
  detail::equal(actual.lower,expected.lower,"counterterm lower endpoint identity failed");
  detail::equal(actual.weight,expected.weight,"counterterm differential identity failed");
  detail::equal(actual.upper,expected.upper,"counterterm upper endpoint identity failed");
  if(!epsilon_regular(actual))throw std::invalid_argument("counterterm leaves epsilon poles");
}
inline CountertermCertificate certify_counterterm(const Functional& f,const Matrix& c) {
  CountertermCertificate result{c,rewrite_counterterm(f,c)};
  verify_counterterm(f,result);return result;
}

struct SearchOptions {
  unsigned pole_order=2, spatial_degree=2;
  std::size_t max_unknowns=64,max_equations=2048,max_terms=100000;
  double seconds=2;
};
enum class SearchStatus { certified, no_certificate_in_ansatz, budget_exceeded, unsupported };
struct SearchResult {
  SearchStatus status=SearchStatus::no_certificate_in_ansatz;
  std::optional<CountertermCertificate> certificate;
  std::size_t unknowns_per_row=0,equations=0,visited_terms=0;
  std::string reason;
};
namespace detail {
struct BudgetExceeded : std::runtime_error { using std::runtime_error::runtime_error; };
struct SearchBudget {
  SearchOptions options;
  std::chrono::steady_clock::time_point start=std::chrono::steady_clock::now();
  std::size_t terms=0,equations=0;
  void check()const {
    if(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>options.seconds)
      throw BudgetExceeded("counterterm cooperative time budget exceeded");
  }
  void count_terms(std::size_t n) {
    check();if(n>options.max_terms-terms)throw BudgetExceeded("counterterm term budget exceeded");terms+=n;
  }
};
// f[0] + sum u_j*f[j+1] must be epsilon-regular. A common polynomial
// denominator D reduces this EXACTLY to vanishing numerator monomials below
// val_epsilon(D); no sampled point, truncated expansion or modular rank is used.
inline void regular_constraints(const std::vector<Exact>& f,std::size_t epsilon,
                                Matrix& equations,SearchBudget& budget) {
  auto denominator=f[0].constant(1);
  for(const auto& v:f) {
    budget.check();denominator=denominator.polynomial_lcm(v.denominator());
    budget.count_terms(denominator.numerator_terms().size());
  }
  const auto pole=denominator.minimum_exponent(epsilon);
  if(!pole)return;
  std::map<std::vector<unsigned long>,std::vector<Exact>> rows;
  for(std::size_t j=0;j<f.size();++j) {
    budget.check();auto polynomial=f[j]*denominator;
    if(!polynomial.denominator().is_one())throw std::logic_error("counterterm denominator clearing failed");
    auto terms=polynomial.numerator_terms();budget.count_terms(terms.size());
    for(const auto& t:terms)if(t.powers[epsilon]<pole) {
      auto [it,inserted]=rows.try_emplace(t.powers,std::vector<Exact>(f.size(),f[0].constant(0)));
      if(inserted && rows.size()>budget.options.max_equations-budget.equations)
        throw BudgetExceeded("counterterm equation budget exceeded");
      // Last column is the augmented right-hand side.
      const auto column=j?j-1:f.size()-1;
      it->second[column]=it->second[column]+f[0].constant(j?t.coefficient:-t.coefficient);
    }
  }
  for(auto& [monomial,row]:rows) { equations.push_back(std::move(row));++budget.equations; }
}
inline std::optional<std::vector<Exact>> solve_constants(Matrix a,std::size_t unknowns,
                                                        const Exact& z,SearchBudget& budget) {
  std::vector<std::size_t> pivots;std::size_t row=0;
  for(std::size_t c=0;c<unknowns && row<a.size();++c) {
    budget.check();std::size_t p=row;while(p<a.size() && a[p][c].is_zero())++p;
    if(p==a.size())continue;std::swap(a[p],a[row]);auto pivot=a[row][c];
    for(auto& v:a[row])v=v/pivot;
    for(std::size_t i=0;i<a.size();++i)if(i!=row && !a[i][c].is_zero()) {
      budget.check();auto factor=a[i][c];
      for(std::size_t j=c;j<=unknowns;++j)a[i][j]=a[i][j]-factor*a[row][j];
    }
    pivots.push_back(c);++row;
  }
  for(const auto& r:a) {
    bool zero=true;for(std::size_t j=0;j<unknowns;++j)zero&=r[j].is_zero();
    if(zero && !r.back().is_zero())return std::nullopt;
  }
  std::vector<Exact> result(unknowns,z.constant(0));
  for(std::size_t i=0;i<pivots.size();++i)result[pivots[i]]=a[i].back();
  return result; // Free rational coefficients are deterministically set to zero.
}
} // namespace detail

// Search C_ij = sum u_ijpr*x^p/(D(x)*epsilon^r), u rational constants.
// D is an explicitly supplied, epsilon-independent allowed spatial denominator.
// Failure is scoped to this finite ansatz. Budgets are cooperative: one exact
// FLINT operation can outlast a check; callers needing a hard deadline must isolate it.
inline SearchResult search_counterterm(const Functional& f,const Exact& denominator,
                                       const SearchOptions& options={}) {
  validate(f);const auto& z=detail::sample(f.connection);z.require_same_field(denominator);
  if(!options.pole_order || options.pole_order>32 || options.spatial_degree>64 ||
     !options.max_unknowns || !options.max_equations || !options.max_terms ||
     !std::isfinite(options.seconds) || options.seconds<=0 || options.seconds>60)
    throw std::invalid_argument("counterterm search options/budget");
  if(denominator.is_zero() || !denominator.is_univariate(f.x) ||
     !denominator.denominator().is_rational())
    throw std::invalid_argument("counterterm spatial denominator must be a nonzero univariate polynomial");
  SearchResult result;detail::SearchBudget budget{options};
  try {
    if(!epsilon_regular(f.connection,f.epsilon)) {
      result.status=SearchStatus::unsupported;result.reason="counterterm search requires epsilon-regular connection";return result;
    }
    const auto n=f.connection.size();
    result.unknowns_per_row=n*(options.spatial_degree+1)*options.pole_order;
    if(result.unknowns_per_row>options.max_unknowns)throw detail::BudgetExceeded("counterterm unknown budget exceeded");
    std::vector<Matrix> bases;
    for(std::size_t j=0;j<n;++j)for(unsigned r=1;r<=options.pole_order;++r)
      for(unsigned p=0;p<=options.spatial_degree;++p) {
        budget.check();auto c=fuchsify::detail::zeros(1,n,z);
        c[0][j]=z.variable(f.x).pow(p)/(denominator*z.variable(f.epsilon).pow(r));
        // The ordinary-cut ansatz must itself have finite endpoint evaluations.
        (void)detail::at(c,f.x,f.lower_point);(void)detail::at(c,f.x,f.upper_point);
        bases.push_back(std::move(c));
      }
    std::vector<Matrix> bulk,left,right;
    for(const auto& c:bases) {
      bulk.push_back(fuchsify::detail::subtract(fuchsify::detail::zeros(1,n,z),
        detail::add(fuchsify::detail::derivative(c,f.x),fuchsify::detail::multiply(c,f.connection))));
      left.push_back(fuchsify::detail::subtract(fuchsify::detail::zeros(1,n,z),detail::at(c,f.x,f.lower_point)));
      right.push_back(detail::at(c,f.x,f.upper_point));
    }
    auto c=fuchsify::detail::zeros(f.weight.size(),n,z);
    for(std::size_t row=0;row<f.weight.size();++row) {
      Matrix equations;
      for(unsigned part=0;part<3;++part)for(std::size_t j=0;j<n;++j) {
        const auto& original=part==0?f.weight:part==1?f.lower:f.upper;
        const auto& coefficients=part==0?bulk:part==1?left:right;
        std::vector<Exact> expression{original[row][j]};
        for(const auto& b:coefficients)expression.push_back(b[0][j]);
        detail::regular_constraints(expression,f.epsilon,equations,budget);
      }
      auto solution=detail::solve_constants(std::move(equations),bases.size(),z,budget);
      if(!solution) { result.reason="no certificate in specified rational principal-part ansatz";break; }
      for(std::size_t v=0;v<bases.size();++v)for(std::size_t j=0;j<n;++j)
        c[row][j]=c[row][j]+(*solution)[v]*bases[v][0][j];
      if(row+1==f.weight.size()) {
        budget.check();result.certificate=certify_counterterm(f,c);
        result.status=SearchStatus::certified;result.reason="exact bulk and both endpoint principal parts removed";
      }
    }
  } catch(const detail::BudgetExceeded& e) { result.status=SearchStatus::budget_exceeded;result.reason=e.what(); }
    catch(const std::domain_error& e) { result.status=SearchStatus::unsupported;result.reason=e.what(); }
  result.equations=budget.equations;result.visited_terms=budget.terms;
  return result;
}

struct Segment {
  Matrix connection,projection,reduced_connection;
  Exact start,end;
};
struct Interface { Matrix original,reduced; };
struct RectangularRealization {
  std::vector<Segment> segments;
  std::vector<Interface> interfaces;
  Matrix source_initialization,reduced_initialization,output,reduced_output;
  std::size_t x=0,epsilon=1;
};
// Full-state identities only: no unproved annihilator or sampled rank shortcut.
// Point coordinates are ordinary rational cuts. Singular endpoint/branch and
// geometry certificates remain separate obligations of the caller.
inline void verify_regular_realization(const RectangularRealization& r) {
  if(r.segments.empty() || r.segments.size()>256 || r.interfaces.size()+1!=r.segments.size())
    throw std::invalid_argument("realization segment/interface count");
  const auto& z=detail::sample(r.segments.front().connection);detail::indices(z,r.x,r.epsilon);
  for(const auto& s:r.segments) {
    const auto n=s.connection.size(),m=s.projection.size();
    detail::shape(s.connection,n,n,z);detail::shape(s.projection,m,n,z);detail::shape(s.reduced_connection,m,m,z);
    // Also reject undefined evaluation of a singular projection at either cut.
    (void)detail::at(s.projection,r.x,s.start);(void)detail::at(s.projection,r.x,s.end);
    if(!epsilon_regular(s.reduced_connection,r.epsilon))throw std::invalid_argument("realization reduced connection has epsilon pole");
    detail::equal(detail::add(fuchsify::detail::derivative(s.projection,r.x),
        fuchsify::detail::multiply(s.projection,s.connection)),
      fuchsify::detail::multiply(s.reduced_connection,s.projection),"realization differential identity failed");
  }
  for(std::size_t i=0;i<r.interfaces.size();++i) {
    const auto& a=r.segments[i];const auto& b=r.segments[i+1];const auto& j=r.interfaces[i];
    detail::shape(j.original,b.connection.size(),a.connection.size(),z);
    detail::shape(j.reduced,b.projection.size(),a.projection.size(),z);
    detail::constant_in_x(j.original,r.x);detail::constant_in_x(j.reduced,r.x);
    if(!epsilon_regular(j.reduced,r.epsilon))throw std::invalid_argument("realization reduced interface has epsilon pole");
    detail::equal(fuchsify::detail::multiply(detail::at(b.projection,r.x,b.start),j.original),
      fuchsify::detail::multiply(j.reduced,detail::at(a.projection,r.x,a.end)),"realization interface identity failed");
  }
  const auto& first=r.segments.front();const auto& last=r.segments.back();
  (void)detail::sample(r.source_initialization);
  const auto source_columns=r.source_initialization.front().size(),outputs=r.output.size();
  detail::shape(r.source_initialization,first.connection.size(),source_columns,z);
  detail::shape(r.reduced_initialization,first.projection.size(),source_columns,z);
  detail::shape(r.output,outputs,last.connection.size(),z);
  detail::shape(r.reduced_output,outputs,last.projection.size(),z);
  for(const auto* a:{&r.source_initialization,&r.reduced_initialization,&r.output,&r.reduced_output})detail::constant_in_x(*a,r.x);
  if(!epsilon_regular(r.reduced_initialization,r.epsilon) || !epsilon_regular(r.reduced_output,r.epsilon))
    throw std::invalid_argument("realization reduced initialization/output has epsilon pole");
  detail::equal(fuchsify::detail::multiply(detail::at(first.projection,r.x,first.start),r.source_initialization),
    r.reduced_initialization,"realization initialization identity failed");
  detail::equal(r.output,fuchsify::detail::multiply(r.reduced_output,detail::at(last.projection,r.x,last.end)),
    "realization output identity failed");
}
} // namespace diffexp::ft_functional_reduction

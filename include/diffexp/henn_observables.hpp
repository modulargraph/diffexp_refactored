#pragma once
#include "diffexp/ibp.hpp"
#include <array>
#include <map>
#include <set>

namespace diffexp::henn {
using Integral=ibp::Integral;

// Q(eps)[eps5]/(eps5^2+3), with the specified X0 embedding eps5=+i*sqrt(3).
// A pair avoids treating eps5 as an independent transcendental variable.
struct Coefficient {
  Exact rational, parity_odd;
  explicit Coefficient(const Exact& value):rational(value),parity_odd(value.constant(0)) {}
  Coefficient(Exact a,Exact b):rational(std::move(a)),parity_odd(std::move(b)) {}
  bool is_zero() const {return rational.is_zero() && parity_odd.is_zero();}
  friend Coefficient operator+(const Coefficient& a,const Coefficient& b) {
    return {a.rational+b.rational,a.parity_odd+b.parity_odd};
  }
  friend Coefficient operator-(const Coefficient& a) {return {-a.rational,-a.parity_odd};}
  friend Coefficient operator-(const Coefficient& a,const Coefficient& b) {return a+-b;}
  friend Coefficient operator*(const Coefficient& a,const Coefficient& b) {
    return {a.rational*b.rational-a.parity_odd*b.parity_odd*a.rational.constant(3),
            a.rational*b.parity_odd+a.parity_odd*b.rational};
  }
  friend Coefficient operator/(const Coefficient& a,const Coefficient& b) {
    auto norm=b.rational*b.rational+b.parity_odd*b.parity_odd*b.rational.constant(3);
    if(norm.is_zero())throw std::domain_error("zero Henn algebraic coefficient denominator");
    auto product=a*Coefficient(b.rational,-b.parity_odd);
    return {product.rational/norm,product.parity_odd/norm};
  }
  Coefficient pow(long n) const {
    if(n < -32 || n > 32)throw std::invalid_argument("Henn coefficient exponent exceeds limit");
    auto result=Coefficient(rational.constant(1)),base=*this;
    if(n<0) {base=result/base;n=-n;}
    while(n) {if(n&1)result=result*base;n>>=1;if(n)base=base*base;}
    return result;
  }
};
using Observable=std::map<Integral,Coefficient>;

struct CanonicalBasis {
  ExactField field{{"eps"}};
  std::vector<Observable> components;
  // Includes ALL syntactic XB targets, including terms vanishing at X0.
  // FIRE discovery should probe this union before selecting observables.
  std::vector<Integral> scalar_targets;
  std::vector<Integral> nonzero_targets;
  static constexpr unsigned physical_slots=8, total_slots=11;
  static constexpr std::array<int,5> x0_adjacent_invariants={3,-1,1,1,-1};
  static constexpr const char* invariant_order="s12,s23,s34,s45,s15";
  static constexpr std::array<unsigned,3> numerator_slots={9,10,11}; // 1-based
  static constexpr int dimension_constant=4, dimension_epsilon_slope=-2;
  static constexpr int normalization_epsilon_power=4;
  static constexpr int normalization_euler_gamma_exponent=2;
  static constexpr int delta_prescription_sign=-1;
  static constexpr const char* algebraic_embedding="eps5 = +I*Sqrt[3]; eps5^2 = -3";
  static constexpr const char* scalar_convention="D_i = -q_i^2; coefficients include (-1)^Sum(indices)";
  static constexpr const char* final_normalization="eps^4 * Exp[2*eps*EulerGamma]";
};

namespace detail {
struct Linear {
  Coefficient scalar;
  Observable terms;
  bool contains_integral=false; // Syntactic: cancellation must not hide products.
  explicit Linear(const Exact& zero):scalar(zero) {}
};
inline void add_term(Observable& terms,const Integral& integral,const Coefficient& coefficient) {
  auto [it,inserted]=terms.emplace(integral,coefficient);
  if(!inserted)it->second=it->second+coefficient;
}
class Importer {
 public:
  explicit Importer(const ExactField& field):zero_(field,0) {}
  std::set<Integral> targets;
  Observable observable(const data::Expr& expression) {
    auto result=parse(expression,0);
    if(!result.scalar.is_zero())throw std::invalid_argument("Henn observable has a non-integral offset");
    for(auto it=result.terms.begin();it!=result.terms.end();) {
      if(it->second.is_zero())it=result.terms.erase(it);else ++it;
    }
    return std::move(result.terms);
  }
 private:
  Exact zero_;
  std::size_t nodes_=0;
  Linear scalar(const Exact& value) {Linear result(zero_);result.scalar=Coefficient(value);return result;}
  Linear scale(Linear value,const Coefficient& coefficient) {
    value.scalar=value.scalar*coefficient;
    for(auto& [integral,c]:value.terms)c=c*coefficient;
    return value;
  }
  Linear parse(const data::Expr& e,unsigned depth) {
    if(depth>512 || ++nodes_>200000)throw std::invalid_argument("Henn observable expression exceeds resource limit");
    if(e.head=="XB") {
      if(e.args.size()!=11)throw std::invalid_argument("Henn XB requires all 11 slots");
      Integral integral;long sum=0;
      for(std::size_t i=0;i<11;++i) {
        auto n=data::integer(e.args[i]);
        if(n < -1024 || n > 1024 || (i>=8 && n>0))
          throw std::invalid_argument("invalid Henn index or positive irreducible-numerator slot");
        integral.push_back(static_cast<int>(n));sum+=n;
      }
      targets.insert(integral);
      Linear value(zero_);value.contains_integral=true;
      value.terms.emplace(integral,Coefficient(zero_.constant(sum%2 ? -1 : 1)));
      return value;
    }
    if(e.atom()) {
      if(data::number(e)) {
        // Decimal/precision-marked literals are not exact ancillary coefficients.
        if(e.head.find_first_not_of("0123456789")!=std::string::npos)
          throw std::invalid_argument("Henn coefficients require exact integer literals");
        return scalar(zero_.constant(Rational(e.head)));
      }
      if(e.head=="s12")return scalar(zero_.constant(3));
      if(e.head=="s23" || e.head=="s15")return scalar(zero_.constant(-1));
      if(e.head=="s34" || e.head=="s45")return scalar(zero_.constant(1));
      if(e.head=="eps")return scalar(zero_.variable(0));
      if(e.head=="d")return scalar(zero_.constant(4)-zero_.constant(2)*zero_.variable(0));
      if(e.head=="eps5") {Linear value(zero_);value.scalar.parity_odd=zero_.constant(1);return value;}
      throw std::invalid_argument("unsupported Henn coefficient symbol: "+e.head);
    }
    if(e.head=="neg" && e.args.size()==1)
      return scale(parse(e.args[0],depth+1),Coefficient(zero_.constant(-1)));
    if(e.args.size()!=2 || (e.head!="+" && e.head!="-" && e.head!="*" && e.head!="/" && e.head!="^"))
      throw std::invalid_argument("unsupported Henn coefficient function: "+e.head);
    auto a=parse(e.args[0],depth+1);
    if(e.head=="^") {
      const auto power=data::integer(e.args[1]);
      if(a.contains_integral) {
        if(power!=1)throw std::invalid_argument("nonlinear Henn integral power");
        return a;
      }
      a.scalar=a.scalar.pow(power);return a;
    }
    auto b=parse(e.args[1],depth+1);
    if(e.head=="+" || e.head=="-") {
      if(e.head=="-")b=scale(std::move(b),Coefficient(zero_.constant(-1)));
      a.scalar=a.scalar+b.scalar;
      for(const auto& [integral,c]:b.terms)add_term(a.terms,integral,c);
      a.contains_integral=a.contains_integral||b.contains_integral;return a;
    }
    if(e.head=="/") {
      if(b.contains_integral)throw std::invalid_argument("Henn integral in denominator");
      return scale(std::move(a),Coefficient(zero_.constant(1))/b.scalar);
    }
    if(a.contains_integral && b.contains_integral)
      throw std::invalid_argument("nonlinear product of Henn integrals");
    if(a.contains_integral)return scale(std::move(a),b.scalar);
    if(b.contains_integral)return scale(std::move(b),a.scalar);
    a.scalar=a.scalar*b.scalar;return a;
  }
};
} // namespace detail

inline CanonicalBasis import_x0(const data::Expr& basis) {
  if(basis.head!="List" || basis.args.size()!=108)
    throw std::invalid_argument("Henn canonical ancillary basis must contain exactly 108 observables");
  CanonicalBasis result;detail::Importer importer(result.field);
  std::set<Integral> nonzero;
  for(const auto& expression:basis.args) {
    auto observable=importer.observable(expression);
    for(const auto& [integral,c]:observable)nonzero.insert(integral);
    result.components.push_back(std::move(observable));
  }
  result.scalar_targets.assign(importer.targets.begin(),importer.targets.end());
  result.nonzero_targets.assign(nonzero.begin(),nonzero.end());
  return result;
}
inline CanonicalBasis read_x0(const std::string& path) {
  std::ifstream input(path,std::ios::binary);
  if(!input)throw std::runtime_error("cannot read Henn ancillary file: "+path);
  input.seekg(0,std::ios::end);const auto size=input.tellg();input.seekg(0);
  if(size<0 || size>2000000)throw std::invalid_argument("Henn ancillary file exceeds 2 MB limit");
  std::string text(static_cast<std::size_t>(size),'\0');input.read(text.data(),size);
  if(!input)throw std::runtime_error("failed reading Henn ancillary data");
  // Bound nesting before invoking the general data reader. The importer then
  // bounds AST depth/node count and accepts only the arithmetic whitelist.
  long depth=0;
  for(char c:text) {
    if(c=='(' || c=='[' || c=='{') {if(++depth>256)throw std::invalid_argument("Henn ancillary nesting exceeds limit");}
    if(c==')' || c==']' || c=='}')--depth;
  }
  return import_x0(data::Reader(std::move(text)).read());
}
} // namespace diffexp::henn

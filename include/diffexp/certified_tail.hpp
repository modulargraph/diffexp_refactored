#pragma once

#include "diffexp/jet.hpp"
#include "diffexp/kernel/local_solution.hpp"
#include <limits>

namespace diffexp {
using NativeTailMagnitude = kernel::Magnitude;
using NativeTailBoundary = std::vector<std::vector<kernel::ComplexBall>>;

// Coefficients of A(z,epsilon). With dlog=true, coefficient denotes a
// nonvanishing letter W and the matrix entry is multiplier * W'/W.
struct NativeAnalyticEntry {
  unsigned row, column, epsilon;
  data::Expr coefficient;
  bool dlog = false;
  kernel::ComplexBall multiplier = kernel::ComplexBall(1);
};

namespace native_tail_detail {
using B = kernel::ComplexBall;
using M = NativeTailMagnitude;
struct AD { B value, derivative; };
inline AD add(const AD& a,const AD& b) {return {a.value+b.value,a.derivative+b.derivative};}
inline AD neg(const AD& a) {return {-a.value,-a.derivative};}
inline AD mul(const AD& a,const AD& b) {
  return {a.value*b.value,a.derivative*b.value+a.value*b.derivative};
}
inline AD divide(const AD& a,const AD& b) {
  if(b.value.contains_zero())throw std::domain_error("analytic disk denominator not separated from zero");
  return {a.value/b.value,(a.derivative-(a.value/b.value)*b.derivative)/b.value};
}
inline bool avoids_principal_cut(const B& b) {
  return arb_is_positive(acb_realref(b.raw())) || !arb_contains_zero(acb_imagref(b.raw()));
}
inline AD eval(const data::Expr& e,const B& disk,const std::string& variable) {
  if(data::number(e))return {Jet(0,1,B::precision()).decimal(e.head).at(0),B(0)};
  if(e.atom()) {
    if(e.head==variable)return {disk,B(1)};
    if(e.head=="I")return {B::from_strings("0","1"),B(0)};
    if(e.head=="Pi") {B p;acb_const_pi(p.raw(),B::precision());return {p,B(0)};}
    throw std::invalid_argument("unbound analytic disk variable: "+e.head);
  }
  if(e.head=="neg" && e.args.size()==1)return neg(eval(e.args[0],disk,variable));
  if((e.head=="Sqrt" || e.head=="Log") && e.args.size()==1) {
    auto a=eval(e.args[0],disk,variable);
    if(!avoids_principal_cut(a.value))
      throw std::domain_error("analytic disk does not certify the principal algebraic/logarithmic sheet");
    B v;
    if(e.head=="Sqrt") {
      acb_sqrt(v.raw(),a.value.raw(),B::precision());
      if(v.contains_zero())throw std::domain_error("analytic disk square root contains zero");
      return {v,a.derivative/(B(2)*v)};
    }
    acb_log(v.raw(),a.value.raw(),B::precision());return {v,a.derivative/a.value};
  }
  if(e.args.size()!=2)throw std::invalid_argument("unsupported analytic disk expression");
  auto a=eval(e.args[0],disk,variable);
  if(e.head=="^") {
    auto n=data::integer(e.args[1]);
    if(n < -100000 || n > 100000)throw std::invalid_argument("analytic disk exponent limit");
    bool inverse=n<0;if(inverse)n=-n;AD out{B(1),B(0)};
    while(n) {if(n&1)out=mul(out,a);n>>=1;if(n)a=mul(a,a);}
    return inverse?divide({B(1),B(0)},out):out;
  }
  auto b=eval(e.args[1],disk,variable);
  if(e.head=="+")return add(a,b);
  if(e.head=="-")return add(a,neg(b));
  if(e.head=="*")return mul(a,b);
  if(e.head=="/")return divide(a,b);
  throw std::invalid_argument("unsupported analytic disk operation: "+e.head);
}
inline void require_finite(const B& b) {
  if(!b.is_finite())throw std::domain_error("nonfinite analytic disk data");
}
} // namespace native_tail_detail

// An unforgeable whole-disk proof, constructed from expressions, not sampled
// jets. Interval arithmetic covers the enclosing square, hence the entire
// closed disk. Failure is inconclusive (reduce R), never a convergence claim.
// Sqrt/Log refer only to their principal analytic sheets on that square.
class NativeAnalyticDisk {
 public:
  static NativeAnalyticDisk certify(const std::vector<NativeAnalyticEntry>& entries,
      unsigned dimension,unsigned epsilon_order,const kernel::ComplexBall& center,
      const std::string& radius_exact,const std::string& variable="x") {
    using namespace native_tail_detail;
    if(!dimension || variable.empty() || variable=="I" || variable=="Pi")
      throw std::invalid_argument("invalid analytic disk shape/variable");
    kernel::Rational radius(radius_exact);
    if(!(kernel::Rational(0)<radius))throw std::invalid_argument("analytic disk radius must be positive");
    require_finite(center);
    NativeAnalyticDisk result;
    result.center_=center;result.radius_exact_=radius.str();result.dimension_=dimension;
    result.scalar_accumulator_=dimension==2;
    B r=B::from_strings(radius.str());require_finite(r);
    result.radius_upper_=M::upper_abs(r);result.radius_lower_=M::lower_abs(r);
    if(result.radius_lower_.is_zero())throw std::domain_error("analytic disk radius precision insufficient");
    B disk=center;result.radius_upper_.add_error_to(disk);
    std::vector<std::vector<M>> sums(epsilon_order+1,std::vector<M>(dimension));
    for(const auto& e:entries) {
      if(e.row>=dimension || e.column>=dimension)throw std::invalid_argument("analytic matrix entry index");
      // Higher epsilon powers cannot affect any coefficient in this window.
      if(e.epsilon>epsilon_order)continue;
      require_finite(e.multiplier);
      auto a=eval(e.coefficient,disk,variable);
      auto value=e.dlog?divide({a.derivative,B(0)},{a.value,B(0)}).value:a.value;
      value=value*e.multiplier;require_finite(value);
      const auto bound=M::upper_abs(value);
      sums[e.epsilon][e.row]+=bound;
      if(e.column!=0 && !bound.is_zero())result.scalar_accumulator_=false;
    }
    for(const auto& rows:sums) {
      M bound;for(const auto& b:rows)bound=M::maximum(bound,b);
      if(!bound.is_finite())throw std::domain_error("nonfinite analytic matrix norm");
      result.norms_.push_back(bound);
    }
    if(result.scalar_accumulator_)for(const auto& rows:sums) {
      result.source_norms_.push_back(rows[0]);result.forcing_norms_.push_back(rows[1]);
    }
    return result;
  }
  bool scalar_accumulator()const{return scalar_accumulator_;}
  const std::vector<NativeTailMagnitude>& source_norms()const{return source_norms_;}
  const std::vector<NativeTailMagnitude>& forcing_norms()const{return forcing_norms_;}
  // The accumulator is one-way forcing, not an exponential growth rate.
  const NativeTailMagnitude& ordinary_growth_norm()const{return scalar_accumulator_?source_norms_[0]:norms_[0];}
  unsigned dimension()const{return dimension_;}
  unsigned epsilon_order()const{return static_cast<unsigned>(norms_.size()-1);}
  const std::vector<NativeTailMagnitude>& matrix_norms()const{return norms_;}
  const kernel::ComplexBall& center()const{return center_;}
  const std::string& radius_exact()const{return radius_exact_;}
  const NativeTailMagnitude& radius_upper()const{return radius_upper_;}
  const NativeTailMagnitude& radius_lower()const{return radius_lower_;}
 private:
  NativeAnalyticDisk()=default;
  unsigned dimension_=0;
  kernel::ComplexBall center_;
  std::string radius_exact_;
  NativeTailMagnitude radius_upper_,radius_lower_;
  std::vector<NativeTailMagnitude> norms_,source_norms_,forcing_norms_;
  bool scalar_accumulator_=false;
};

struct NativeTaylorTail {
  // These bound only omitted Taylor terms, uniformly over components.
  // Arithmetic/input uncertainty remains in the independently computed
  // polynomial balls and is never relabeled as Taylor error.
  std::vector<NativeTailMagnitude> absolute;
  std::vector<NativeTailMagnitude> circle_upper;
  unsigned retained_order=0;
  std::string radius_exact;
};

// PRECONDITION: initial encloses the true y(center) coefficients, including
// all earlier truncation errors. A previous finite jet without its remainder
// does not satisfy this contract. Correlated inputs remain conservative balls;
// no independence or cancellation is assumed. No negative epsilon powers of
// A are allowed. N denotes the last retained Taylor degree, inclusively.
//
// Proof: radial integral comparison bounds y_k on |z-center|=R by the
// coefficient of exp(R sum_j M_j eps^j) * sum_k ||initial_k|| eps^k.
// The j=0 factor is exp(R M0); positive-epsilon contributions are a finite
// convolution, so a canonical system needs only (R M1)^j/j!, not exp(R M1).
// Cauchy's inequality and geometric summation then give the omitted tail.
inline NativeTaylorTail certify_native_taylor_tail(const NativeAnalyticDisk& disk,
    const NativeTailBoundary& initial,const kernel::ComplexBall& step,unsigned retained_order) {
  using namespace native_tail_detail;
  if(initial.size()!=disk.dimension())throw std::invalid_argument("tail boundary dimension mismatch");
  const unsigned depth=disk.epsilon_order();
  std::vector<M> b(depth+1);
  for(const auto& row:initial) {
    if(row.size()!=depth+1)throw std::invalid_argument("tail boundary epsilon mismatch");
    for(unsigned k=0;k<=depth;++k) {require_finite(row[k]);b[k]=M::maximum(b[k],M::upper_abs(row[k]));}
  }
  require_finite(step);
  auto ratio=M::upper_abs(step)/disk.radius_lower();
  auto gap=M::positive_difference_lower(M::one(),ratio);
  if(gap.is_zero())throw std::domain_error("evaluation is not provably strictly inside analytic disk");
  if(retained_order==std::numeric_limits<unsigned>::max())throw std::invalid_argument("tail order overflow");
  auto geometric=ratio.power_upper(static_cast<ulong>(retained_order)+1)/gap;
  // For y'=a*y, z'=w*y the source is independent of the accumulator.
  // Bound y on the disk using only a, then integrate w*y along each radius:
  // |z_k(z)-z_k(center)| <= R sum_j W_j Y_{k-j}. The constant z(center)
  // has no omitted Taylor coefficients and remains in the polynomial balls.
  // This recognition uses proved whole-disk entry bounds, never small jets.
  if(disk.scalar_accumulator()) {
    const auto& a=disk.source_norms();const auto& w=disk.forcing_norms();
    std::vector<M> exponential(depth+1),source(depth+1);exponential[0]=M::one();
    for(unsigned k=1;k<=depth;++k) {
      for(unsigned j=1;j<=k;++j)exponential[k]+=M::from_ui(j)*disk.radius_upper()*a[j]*exponential[k-j];
      exponential[k]=exponential[k]/M::lower_abs(B::from_strings(std::to_string(k)));
    }
    const auto ordinary=(disk.radius_upper()*a[0]).exponential_upper();
    for(unsigned k=0;k<=depth;++k) {
      for(unsigned j=0;j<=k;++j)source[k]+=M::upper_abs(initial[0][j])*exponential[k-j];
      source[k]=source[k]*ordinary;
    }
    NativeTaylorTail result;result.retained_order=retained_order;result.radius_exact=disk.radius_exact();
    for(unsigned k=0;k<=depth;++k) {
      M variation;for(unsigned j=0;j<=k;++j)variation+=w[j]*source[k-j];variation=variation*disk.radius_upper();
      const auto circle=M::maximum(source[k],M::upper_abs(initial[1][k])+variation);
      if(!circle.is_finite())throw std::domain_error("nonfinite scalar accumulator disk majorant");
      result.circle_upper.push_back(circle);
      result.absolute.push_back(k==0 && disk.matrix_norms()[0].is_zero()?M::zero():M::maximum(source[k],variation)*geometric);
    }
    return result;
  }
  std::vector<M> g(depth+1);g[0]=M::one();
  for(unsigned k=1;k<=depth;++k) {
    for(unsigned j=1;j<=k;++j)
      g[k]+=M::from_ui(j)*disk.radius_upper()*disk.matrix_norms()[j]*g[k-j];
    // k is represented exactly at practical epsilon depths; use a directed
    // lower bound even if a caller supplies a larger integer.
    g[k]=g[k]/M::lower_abs(B::from_strings(std::to_string(k)));
  }
  auto ordinary=(disk.radius_upper()*disk.matrix_norms()[0]).exponential_upper();
  NativeTaylorTail result;result.retained_order=retained_order;result.radius_exact=disk.radius_exact();
  for(unsigned k=0;k<=depth;++k) {
    M circle;for(unsigned j=0;j<=k;++j)circle+=b[j]*g[k-j];circle=circle*ordinary;
    if(!circle.is_finite())throw std::domain_error("nonfinite solution disk majorant");
    result.circle_upper.push_back(circle);
    // With no epsilon-zero connection y_0 is exactly constant.
    result.absolute.push_back(k==0 && disk.matrix_norms()[0].is_zero()?M::zero():circle*geometric);
  }
  return result;
}

// Checks the actual enclosure radius (both Cartesian radii), not midpoint
// agreement or apparent convergence. absolute_tolerance must be positive.
inline bool native_enclosure_meets_tolerance(const NativeTailBoundary& values,
    const std::string& absolute_tolerance) {
  using namespace native_tail_detail;
  kernel::Rational tolerance(absolute_tolerance);
  if(!(kernel::Rational(0)<tolerance))throw std::invalid_argument("nonpositive enclosure tolerance");
  auto limit=M::lower_abs(B::from_strings(tolerance.str()));
  if(values.empty())throw std::invalid_argument("empty enclosure");
  for(const auto& row:values) {
    if(row.empty())throw std::invalid_argument("empty enclosure epsilon window");
    for(const auto& value:row) {
      if(!value.is_finite())return false;
      B error;arf_set_mag(arb_midref(acb_realref(error.raw())),arb_radref(acb_realref(value.raw())));
      arf_set_mag(arb_midref(acb_imagref(error.raw())),arb_radref(acb_imagref(value.raw())));
      if(!(M::upper_abs(error)<=limit))return false;
    }
  }
  return true;
}
} // namespace diffexp

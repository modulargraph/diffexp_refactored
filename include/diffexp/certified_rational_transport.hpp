#pragma once
#include "diffexp/rational_transport.hpp"
#include "diffexp/certified_tail.hpp"
#include <optional>

namespace diffexp {
struct CertifiedRationalOptions {
  slong working_bits=256;
  unsigned requested_digits=20;
  unsigned taylor_order=112;
  unsigned max_charts=512;
  unsigned max_disk_attempts=32;
  std::string initial_radius="1/4";
};
struct CertifiedRationalResult {
  Boundary boundary;
  unsigned charts=0;
  unsigned disk_attempts=0;
  std::string absolute_tolerance;
  // This is an explicit theorem precondition, not an assertion that arbitrary
  // caller-supplied finite jets have acquired a certificate.
  bool conditional_on_initial_enclosure=true;
};

// Transport y'=A(x,epsilon)y from x=0 to x=1. PRECONDITION: initial encloses
// the true initial coefficients, including every earlier omitted-series error.
// A is a finite polynomial in nonnegative epsilon powers; callers expand any
// rational epsilon dependence to the needed order before invoking this API.
// Epsilon indices are relative to the supplied window, which may start below
// zero. Coefficients may be complex rational functions of x and I. A mapped
// complex contour leg is represented by its already-pulled-back connection.
//
// Every accepted disk proves holomorphy on the whole complex disk, every
// retained polynomial is augmented with its rigorous tail, and the final
// actual radii must meet 10^(-requested_digits). Incoming uncertainty is never
// replaced by a midpoint. The absolute guarantee remains conditional on the
// input enclosure precondition above. Interior poles/insufficient work or
// precision budgets fail rather than returning an uncertified partial path.
inline CertifiedRationalResult certified_rational_line(
    const std::vector<RationalLineEntry>& entries,const Boundary& initial,
    const CertifiedRationalOptions& options={}) {
  using B=kernel::ComplexBall;using M=NativeTailMagnitude;
  if(options.working_bits<64 || options.working_bits>1000000 ||
      !options.requested_digits || options.requested_digits>10000 ||
      !options.taylor_order || options.taylor_order>10000 ||
      !options.max_charts || !options.max_disk_attempts)
    throw std::invalid_argument("invalid certified rational precision or finite work budget");
  if(initial.empty() || initial.size()>100000 || initial[0].empty() || initial[0].size()>1001)
    throw std::invalid_argument("invalid certified rational boundary shape");
  const unsigned dimension=static_cast<unsigned>(initial.size());
  const unsigned depth=static_cast<unsigned>(initial[0].size()-1);
  for(const auto& row:initial) {
    if(row.size()!=depth+1)throw std::invalid_argument("certified rational boundary epsilon mismatch");
    for(const auto& value:row)if(!value.is_finite())throw std::invalid_argument("nonfinite certified rational initial enclosure");
  }
  const Rational requested_radius(options.initial_radius);
  if(!(Rational(0)<requested_radius))throw std::invalid_argument("certified rational radius must be positive");
  struct PrecisionScope {
    slong previous;
    explicit PrecisionScope(slong bits):previous(B::precision()){B::set_precision(bits);}
    ~PrecisionScope(){B::set_precision(previous);}
  } precision(options.working_bits);
  // Parse/validate the recurrence only once; separate exact centered copies
  // are used for disk proofs to avoid interval dependency in denominators.
  const auto compiled=compile_rational_entries(entries);
  for(const auto& entry:entries)if(entry.row>=dimension || entry.column>=dimension)
    throw std::invalid_argument("certified rational matrix entry exceeds boundary dimension");
  CertifiedRationalResult result;result.boundary=initial;
  result.absolute_tolerance="1/1"+std::string(options.requested_digits,'0');
  const auto local_limit=M::lower_abs(B::from_strings("1/1"+std::string(options.requested_digits+20,'0')));
  Rational center(0);
  while(center<Rational(1)) {
    if(result.charts>=options.max_charts)throw std::runtime_error("certified rational chart budget exhausted before endpoint");
    std::vector<NativeAnalyticEntry> centered;
    for(const auto& entry:entries) {
      if(entry.epsilon>depth)continue;
      const auto& coefficient=entry.coefficient;
      std::vector<Exact> point;
      for(std::size_t i=0;i<coefficient.variable_count();++i) {
        auto v=coefficient.variable(i);
        point.push_back(coefficient.variables()[i]=="x"?v+coefficient.constant(center):v);
      }
      centered.push_back({entry.row,entry.column,entry.epsilon,
        data::Reader(coefficient.substitute(point).str()).read()});
    }
    Rational radius=requested_radius,advance(0),step_divisor(4);B step;
    std::optional<NativeTaylorTail> tail;
    for(unsigned attempt=0;attempt<options.max_disk_attempts;++attempt) {
      ++result.disk_attempts;
      bool shrink_radius=true;
      try {
        auto disk=NativeAnalyticDisk::certify(centered,dimension,depth,B(0),radius.str());
        const Rational remaining=Rational(1)-center,quarter=radius/step_divisor;
        advance=remaining<quarter?remaining:quarter;
        if(!(Rational(0)<advance))throw std::runtime_error("certified rational chart failed to make positive progress");
        step=B::from_strings(advance.str());
        auto candidate=certify_native_taylor_tail(disk,result.boundary,step,options.taylor_order);
        bool small=true;for(const auto& error:candidate.absolute)if(!(error<=local_limit))small=false;
        if(small){tail.emplace(std::move(candidate));break;}
        // Shrinking only R leaves the Cauchy factor (h/R)^(N+1) fixed.
        // Reduce that ratio independently, so any positive target can be
        // reached within a sufficient finite chart/proof budget.
        if(disk.radius_upper()*disk.ordinary_growth_norm()<=M::from_ui(8)) {
          step_divisor=step_divisor*Rational(2);shrink_radius=false;
        }
      }catch(const std::domain_error&) {}
      if(shrink_radius)radius=radius/Rational(2);
    }
    if(!tail)throw std::runtime_error("certified rational disk/tail-proof budget exhausted; reduce requested accuracy or increase Taylor order/proof budget");
    auto polynomial=rational_chart(compiled,result.boundary,B::from_strings(center.str()),step,options.taylor_order);
    for(auto& row:polynomial)for(unsigned k=0;k<=depth;++k) {
      tail->absolute[k].add_error_to(row[k]);
      if(!row[k].is_finite())throw std::runtime_error("nonfinite certified rational enclosure");
    }
    result.boundary=std::move(polynomial);center+=advance;++result.charts;
  }
  if(!native_enclosure_meets_tolerance(result.boundary,result.absolute_tolerance))
    throw std::runtime_error("certified rational final enclosure does not meet requested absolute digits; incoming uncertainty, arithmetic precision or Taylor order is insufficient");
  return result;
}
} // namespace diffexp

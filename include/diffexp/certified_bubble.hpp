#pragma once
#include "diffexp/bubble.hpp"
#include "diffexp/certified_tail.hpp"
#include <optional>

namespace diffexp {
struct CertifiedBubbleOptions {
  unsigned epsilon_order=4;
  slong working_bits=256;
  unsigned requested_digits=20;
  unsigned taylor_order=112;
  unsigned max_charts=128;
  unsigned max_disk_attempts=16;
  std::string initial_radius="1/4";
};
struct CertifiedBubbleResult {
  std::vector<kernel::ComplexBall> coefficients;
  unsigned charts=0;
  unsigned disk_attempts=0;
  unsigned requested_digits=0;
  std::string absolute_tolerance;
};

// A complete certificate chain for the Euclidean equal-mass bubble at p^2=-1,
// d=2-2 epsilon. The source boundary is a rigorous gamma-series enclosure.
// Every chart attaches an omitted-Taylor-tail bound before becoming the next
// chart's initial enclosure; no approximate matched boundary enters this path.
// The requested digits mean absolute radius <= 10^(-digits), uniformly over
// returned epsilon coefficients, and are checked on the final actual balls.
inline CertifiedBubbleResult certified_feynman_bubble(const CertifiedBubbleOptions& options={}) {
  using B=kernel::ComplexBall;
  if(options.working_bits<64 || options.working_bits>1000000 ||
      !options.requested_digits || options.requested_digits>10000 ||
      !options.taylor_order || options.taylor_order>100000 || options.epsilon_order>1000 ||
      !options.max_charts || !options.max_disk_attempts)
    throw std::invalid_argument("invalid certified bubble precision or finite work budget");
  Rational initial_radius(options.initial_radius);
  if(!(Rational(0)<initial_radius))throw std::invalid_argument("certified bubble radius must be positive");
  // Existing native scalar arithmetic uses a thread-local ambient precision.
  // Restore it on success and exceptions so this public solver owns its choice.
  struct PrecisionScope {
    slong previous;
    explicit PrecisionScope(slong bits):previous(B::precision()){B::set_precision(bits);}
    ~PrecisionScope(){B::set_precision(previous);}
  } precision(options.working_bits);
  ExactField field({"x","I","eps"});Exact x(field,"x");
  auto family=feynman::banana(1,{Rational(1),Rational(1)});
  auto geometry=feynman::symanzik(family,{x,x.constant(1)-x});
  auto request=feynman::merge_request({1,1},0,1);
  auto quadratic=ibp::quadratic_family(family,x);
  ibp::Generator equations(ibp::PropagatorBasis(ibp::merge(quadratic,0,1,x)),Exact(field,"2-2*eps"));
  ibp::ExactReducer reduction(x);
  for(int power=1;power<=4;++power)for(int isp=0;isp>=-2;--isp)
    for(auto& identity:equations.relations({power,isp}))reduction.insert(std::move(identity));
  auto system=ibp::differential_system(equations,reduction,{{2,0}},0,x);
  auto substitutions=std::vector<Exact>{x,x.variable(1),x.constant(0)};
  auto connection0=system.matrix[0][0].substitute(substitutions);
  auto connection1=system.matrix[0][0].derivative(2).substitute(substitutions);
  if(!system.matrix[0][0].derivative(2).derivative(2).is_zero())
    throw std::logic_error("certified bubble connection unexpectedly nonaffine in epsilon");
  auto origin=substitutions;origin[0]=x.constant(0);
  auto seed=feynman::tadpole(geometry.U.substitute(origin).rational(),geometry.F.substitute(origin).rational(),
    request.source_indices[0],family.loops,2,options.epsilon_order);
  if(seed.epsilon_low!=0)throw std::logic_error("certified bubble gamma boundary unexpectedly polar");
  NativeTailBoundary boundary{seed.coefficients,std::vector<B>(options.epsilon_order+1,B(0))};
  auto integrand=x.constant(request.normalization)*x.pow(request.left_power)*(x.constant(1)-x).pow(request.right_power);
  std::vector<RationalLineEntry> entries{{0,0,0,connection0},{0,0,1,connection1},{1,0,0,integrand}};
  // Parse and validate the exact matrix only once, outside the chart loop.
  const auto compiled=compile_rational_entries(entries);
  std::vector<NativeAnalyticEntry> analytic;
  for(const auto& entry:compiled)analytic.push_back({entry.row,entry.column,entry.epsilon,entry.coefficient});
  CertifiedBubbleResult result;result.requested_digits=options.requested_digits;
  result.absolute_tolerance="1/1"+std::string(options.requested_digits,'0');
  Rational center(0);
  while(center<Rational(1)) {
    if(result.charts>=options.max_charts)throw std::runtime_error("certified bubble chart budget exhausted before endpoint");
    const B c=B::from_strings(center.str());
    Rational radius=initial_radius;
    std::optional<NativeAnalyticDisk> disk;
    for(unsigned attempt=0;attempt<options.max_disk_attempts;++attempt) {
      ++result.disk_attempts;
      try {disk.emplace(NativeAnalyticDisk::certify(analytic,2,options.epsilon_order,c,radius.str()));break;}
      catch(const std::domain_error&) {radius=radius/Rational(2);}
    }
    if(!disk)throw std::runtime_error("certified bubble disk-proof budget exhausted");
    const Rational remaining=Rational(1)-center;
    const Rational half=radius/Rational(2);
    const Rational advance=remaining<half?remaining:half;
    if(!(Rational(0)<advance))throw std::runtime_error("certified bubble failed to make positive progress");
    const B step=B::from_strings(advance.str());
    auto tail=certify_native_taylor_tail(*disk,boundary,step,options.taylor_order);
    auto polynomial=rational_chart(compiled,boundary,c,step,options.taylor_order);
    for(auto& row:polynomial)for(unsigned k=0;k<=options.epsilon_order;++k) {
      tail.absolute[k].add_error_to(row[k]);
      if(!row[k].is_finite())throw std::runtime_error("certified bubble produced nonfinite enclosure");
    }
    boundary=std::move(polynomial);center+=advance;++result.charts;
  }
  if(!native_enclosure_meets_tolerance({boundary[1]},result.absolute_tolerance))
    throw std::runtime_error("certified bubble final enclosure does not meet requested absolute digits; increase working precision or Taylor order");
  result.coefficients=std::move(boundary[1]);return result;
}
} // namespace diffexp

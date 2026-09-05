#pragma once
#include "diffexp/frobenius.hpp"
#include "diffexp/rational_transport.hpp"

namespace diffexp {
// Local transfer across u=0 on the upper bank. The same Frobenius constants
// are used on both sides: arg(u) decreases from +pi to zero. This is a local
// singular-sector match, with no ordinary chart spanning the singularity.
struct FrobeniusCrossing {
  FrobeniusSeries series;
  static constexpr bool omitted_tail_certified = false;
  static FrobeniusCrossing prepare(const FrobeniusSeries::ExactMatrix& local,
      std::size_t xi,std::size_t ei,unsigned order,unsigned epsilon_order) {
    if(local.empty() || local.size()>32 || order>256 || epsilon_order>32)
      throw std::invalid_argument("local singular crossing exceeds finite preparation limits");
    return {FrobeniusSeries::prepare(local,xi,ei,order,epsilon_order)};
  }
  Boundary upper_bank(const Boundary& incoming,const Rational& left_radius,
      const Rational& right_radius) const {
    if(left_radius.sign()<=0 || right_radius.sign()<=0)
      throw std::invalid_argument("singular crossing radii must be positive");
    auto left=Jet::Ball::from_strings((-left_radius).str());
    auto right=Jet::Ball::from_strings(right_radius.str());
    return series.solution(right,series.match(left,incoming));
  }
};

namespace singular_transport_detail {
inline double finite_discrepancy(const Jet::Ball& value,const Jet::Ball& reference) {
  if(!value.is_finite() || !reference.is_finite())
    throw std::domain_error("equal banana endpoint or reference coefficient is nonfinite");
  const auto difference=value-reference;
  if(!difference.is_finite())
    throw std::domain_error("equal banana endpoint discrepancy is nonfinite");
  arf_t bound;arf_init(bound);
  acb_get_abs_ubound_arf(bound,difference.raw(),Jet::Ball::precision());
  const auto result=arf_get_d(bound,ARF_RND_CEIL);arf_clear(bound);
  if(!std::isfinite(result))
    throw std::domain_error("equal banana discrepancy exceeds finite reporting range");
  return result;
}
} // namespace singular_transport_detail

struct OriginalEqualBananaRealResult {
  Boundary value;
  Rational endpoint{20};
  // The independent published table is at20, including when continuing to32.
  double maximum_discrepancy=0,seconds=0;
  unsigned singular_matches=0;
  bool omitted_tail_certified=false;
};
inline OriginalEqualBananaRealResult original_banana_equal_real(const std::string& directory,
    const Rational& endpoint=Rational(20)) {
  using B=Jet::Ball;
  if(endpoint<Rational(20))throw std::invalid_argument("equal banana endpoint must include the saved comparison at20");
  struct PrecisionScope {
    slong previous=B::precision();
    ~PrecisionScope(){B::set_precision(previous);}
  } precision_scope;
  B::set_precision(384);
  auto start=std::chrono::steady_clock::now();
  const auto a0=data::read_file(directory+"/Data/Banana/EqualMass/dt_0.m");
  const auto a1=data::read_file(directory+"/Data/Banana/EqualMass/dt_1.m");
  auto current=read_boundary(original_banana_reference(directory,"boundaryAtMinusOne"),4,4,384);
  const auto reference=read_boundary(original_banana_reference(directory,"referenceAt20"),4,4,384);
  ExactField field({"x","eps"});Exact x(field,"x"),eps(field,"eps");
  auto matrix=[&](const Exact& coordinate,const Exact& derivative) {
    FrobeniusSeries::ExactMatrix result(4,std::vector<Exact>(4,x.constant(0)));
    for(unsigned k=0;k<2;++k) {
      const auto& source=k?a1:a0;
      if(source.args.size()!=4)throw std::invalid_argument("equal banana source dimension");
      for(unsigned i=0;i<4;++i) {
        if(source.args[i].args.size()!=4)throw std::invalid_argument("equal banana source row width");
        for(unsigned j=0;j<4;++j)
          result[i][j]=result[i][j]+evaluate_exact(source.args[i].args[j],x,{{"t",coordinate}})*derivative*(k?eps:x.constant(1));
      }
    }
    return result;
  };
  auto ordinary=[&](const Rational& from,const Rational& to) {
    if(from==to)return;
    const auto width=x.constant(to-from);
    auto pulled=matrix(x.constant(from)+width*x,width);
    std::vector<RationalLineEntry> entries;
    for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j) {
      std::vector<Exact> origin{x,x.constant(0)};
      std::vector<Exact> coefficients{pulled[i][j].substitute(origin),pulled[i][j].derivative(1).substitute(origin)};
      for(unsigned k=0;k<2;++k)if(!coefficients[k].is_zero())entries.push_back({i,j,k,coefficients[k]});
    }
    current=rational_line(entries,std::move(current),50);
  };
  Rational position(-1);
  OriginalEqualBananaRealResult result;
  const std::vector<std::pair<Rational,Rational>> patches{
      {Rational(0),Rational("1/2")},{Rational(4),Rational("1/2")},{Rational(16),Rational(1)}};
  for(const auto& [pole,radius]:patches) {
    ordinary(position,pole-radius);
    auto crossing=FrobeniusCrossing::prepare(matrix(x+x.constant(pole),x.constant(1)),0,1,50,4);
    current=crossing.upper_bank(current,radius,radius);
    position=pole+radius;++result.singular_matches;
  }
  ordinary(position,Rational(20));
  for(unsigned i=0;i<4;++i)for(unsigned k=0;k<5;++k) {
    result.maximum_discrepancy=std::max(result.maximum_discrepancy,
        singular_transport_detail::finite_discrepancy(current[i][k],reference[i][k]));
  }
  // Beyond20 the real path has no further pole. This also reproduces the
  // notebook's original march to32 while retaining its saved table at20.
  ordinary(Rational(20),endpoint);
  for(const auto& row:current)for(const auto& value:row)
    if(!value.is_finite())throw std::domain_error("nonfinite equal banana continued endpoint");
  result.value=std::move(current);result.endpoint=endpoint;
  result.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  return result;
}
inline int run_original_banana_equal_real(const std::string& directory,const Rational& endpoint=Rational(20)) {
  auto result=original_banana_equal_real(directory,endpoint);
  std::cout<<"Original equal banana real +i0 route: singular matches="<<result.singular_matches
    <<", endpoint="<<result.endpoint.str()<<", N50/epsilon4/bits384, reference at20 maximum discrepancy="<<result.maximum_discrepancy
    <<", seconds="<<result.seconds<<"\nOmitted Taylor tails are not certified.\n";
  return std::isfinite(result.maximum_discrepancy) && result.maximum_discrepancy<1e-10?0:1;
}
} // namespace diffexp

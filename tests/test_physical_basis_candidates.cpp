#include "diffexp/physical_basis_scan.hpp"
#include <iostream>
using namespace diffexp;
namespace scan=physical_basis_scan;
void check(bool b,const char* message){if(!b)throw std::runtime_error(message);}
int main(){try {
  ExactField field({"x","eps"});
  auto e=[&](const char*s){return Exact(field,s);};
  level::Result level;
  level.ordered_basis={{1,1},{1,0}};
  level.matrix={{e("0"),e("1")},{e("0"),e("0")}};
  level.requested_integrals={{2,1}};
  level.target_rows={{e("eps^-1"),e("0")}};
  scan::Options options; options.seconds=10; options.max_trials=96;
  options.max_candidates=8; options.max_depth=2; options.use_source_identities=false;
  options.coordinate_rows={{e("eps^-1"),e("0")}};
  bool saw_larger=false,saw_lower_mix=false,saw_depth2=false;
  auto evaluate=[&](const scan::Candidate &c)->std::optional<scan::DemandCost>{
    recursion::verify_physical_basis(c.transform,level.matrix,0);
    check(c.metadata.exact_verified,"missing exact certificate");
    auto requested=fuchsify::detail::multiply(level.target_rows,c.transform.from_basis);
    auto v=exact_epsilon_valuation(requested[0][0],1);
    std::uint64_t loss=v&&*v<0 ? -*v : 0;
    if(c.metadata.connection_characters>scan::expression_characters(level.matrix)&&loss==0)
      saw_larger=true;
    if(!c.transform.to_basis[0][1].is_zero()) saw_lower_mix=true;
    if(c.metadata.depth==2) saw_depth2=true;
    // A real consumer supplies complete recursion cost. Here this deliberately
    // small projected-demand score tests scanner ordering, not endpoint claims.
    return scan::DemandCost{loss};
  };
  auto result=scan::scan_candidates(level,2,0,options,evaluate);
  check(saw_larger,"larger expression discarded before epsilon demand evaluation");
  check(saw_lower_mix,"legitimate lower-sector sparse combination missing");
  check(saw_depth2,"beam never composes two row changes");
  check(!result.candidates.empty()&&result.candidates.size()<=options.max_candidates,"frontier bound");
  for(const auto &c:result.candidates) check(c.demand_cost->at(0)==0,"dominated demand retained");
  check(result.stats.trials<=options.max_trials,"trial bound");
  auto again=scan::scan_candidates(level,2,0,options,evaluate);
  check(again.candidates.size()==result.candidates.size(),"deterministic frontier size");
  for(std::size_t i=0;i<result.candidates.size();++i)
    check(result.candidates[i].transform.to_basis==again.candidates[i].transform.to_basis,"deterministic candidate order");
  // Pareto tradeoffs are preserved; equal demand doesn't imply equal endpoint
  // admissibility, so equal-cost representations must remain available.
  check(!scan::dominates({1,9},{2,8})&&scan::dominates({1,7},{2,8}),"Pareto comparison");
  options.max_depth=1; options.max_trials=20; options.allow_lower_sector_mixes=false;
  options.coordinate_rows.clear(); level.requested_integrals.clear();level.target_rows.clear();
  auto no_mix=scan::scan_candidates(level,2,0,options);
  for(const auto &c:no_mix.candidates) {
    check(c.transform.to_basis[0][1].is_zero()&&c.transform.to_basis[1][0].is_zero(),"unpermitted cross-sector mix");
    recursion::verify_physical_basis(c.transform,level.matrix,0);
  }
  // Source recovery reuses only cached exact identities, including numerators.
  ibp::Relation relation{{{1,-1},e("1")},{{1,1},e("-1")},{{1,0},e("-eps")}};
  level.source_identities={relation}; options.use_source_identities=true;
  auto recovered=scan::scan_candidates(level,2,0,options);
  check(recovered.stats.recovered_coordinates>0&&recovered.stats.numerator_candidates>0,"cached numerator recovery");
  // Pole-introducing supplied coordinates are exactly rejected before scoring.
  options.use_source_identities=false;options.coordinate_rows={{e("1/x"),e("0")}};
  auto poles=scan::scan_candidates(level,2,0,options);
  check(poles.stats.new_poles>0,"new parameter pole admitted");
  auto unresolved=scan::scan_candidates(level,2,0,options,[](const auto&)->std::optional<scan::DemandCost>{return std::nullopt;});
  check(unresolved.candidates.empty()&&!unresolved.unresolved.empty(),"unresolved presented as accepted");
  // Parameter-dependent row changes require the derivative term in the gauge
  // identity. Its pole is already present in the original system.
  level.matrix={{e("1/(1+x)"),e("0")},{e("0"),e("0")}};
  options.coordinate_rows={{e("1/(1+x)"),e("0")}};
  bool derivative_candidate=false;
  scan::scan_candidates(level,2,0,options,[&](const scan::Candidate &c)->std::optional<scan::DemandCost>{
    if(c.transform.to_basis[0][0]==e("1/(1+x)") && c.transform.to_basis[1][1]==e("1")) {
      check(c.transform.connection[0][0].is_zero(),"missing transformation derivative");
      derivative_candidate=true;
    }
    return scan::DemandCost{1};
  });
  check(derivative_candidate,"parameter-dependent candidate missing");
  std::cout<<"Demand-first Pareto basis candidates, exact identities, sparse sectors and cache recovery passed\n";
}catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 1;}}

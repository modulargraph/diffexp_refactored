#include "diffexp/recursion_pipeline.hpp"
#include <iostream>
using namespace diffexp;using B=Jet::Ball;
void check(bool value,const char* why){if(!value)throw std::runtime_error(why);}
int main(){try{
 B::set_precision(192);
 auto graph=recursion::prepare(feynman::example_family("bubble"),{{1,1}});
 auto& node=graph.nodes[0];node.scalar_leaf=false;
 const auto one=graph.dimension.constant(1),zero=graph.dimension.constant(0),x=graph.dimension.variable(0);
 ibp::Integral master(node.merged.denominators.size(),0);master[0]=1;auto dotted=master;dotted[0]=2;
 node.closure.ordered_basis={master};node.closure.success=true;
 node.closure.matrix={{one/(one+x)}};
 // Synthetic cached coordinate, not an asserted physical integral identity:
 // U=Y/(1+x) has U'=0 while Y'=Y/(1+x).
 node.closure.requested_integrals={dotted};node.closure.target_rows={{one/(one+x)}};
 node.observable_rows={{one}};node.operations={pullback::Plan{}};
 node.operations[0].operation=feynman::Operation::UpperLimit;
 causal::Prescription prescription;prescription.f_rim=-1;prescription.levels={{1}};
 prescription.provenance="synthetic endpoint-screening identity";
 recursion::NumericalOptions settings;settings.working_bits=192;settings.endpoint_order=6;
 settings.causal_prescription=prescription;settings.automatic_basis_scan=false;settings.exact_functional_reduction=false;
 recursion::Evaluator original(graph,settings);auto old=original.prepare_factored_consumers(0,3);
 settings.automatic_basis_scan=true;settings.basis_scan_seconds=2;
 recursion::Evaluator search(graph,settings);auto selected=search.prepare_factored_consumers(0,3);
 const auto& selection=search.basis_selections().at(0);
 check(selection.accepted,"cheap frontier failed to promote simpler equivalent basis");
 check(selection.candidates_screened>selection.full_candidates_checked&&selection.full_candidates_checked<=8,
       "endpoint rebuild was not deferred to bounded improving frontier");
 check(selected.connection[0][0].is_zero(),"full-validation candidate differs from screened exact coordinate");
 // The upper operator is expressed at physical x=1-h, not reversed x=h.
 // V=1+x, so Enew=Eold*(2-h)=2. Compare complete retained jets, including zeros.
 auto h=graph.dimension.constant(original.endpoint_geometry(0).overlap);
 auto v=exact_laurent_rows({{one+x}},one-h,3);
 auto composed=linear_boundary::compose(old.right,v,3);
 for(int k=0;k<=3;++k){const auto& lhs=selected.right.coefficients[0][0][k-selected.right.low];
   const auto& rhs=composed.coefficients[0][0][k-composed.low];
   check(acb_overlaps(lhs.raw(),rhs.raw()),"upper-cut coordinate pullback disagrees with rebuilt endpoint");
   B exact(k?0:2);check(acb_contains(lhs.raw(),exact.raw()),"requested physical upper functional changed");}
 check(recursion::basis_persistence::revision==std::string("physical-basis-demand-selection-v3"),"screen revision did not invalidate v2 decisions");
 std::cout<<"Cheap endpoint screens, bounded full validation, upper-cut identity and v3 revision passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

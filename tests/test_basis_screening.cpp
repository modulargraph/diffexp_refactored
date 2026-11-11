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
 // The Direct pole is a child epsilon loss that is the same in every
 // internal basis: the candidate must be adopted at an unchanged child window.
 const auto eps=graph.dimension.variable(1);
 node.observable_rows={{one},{one/eps}};node.operations={pullback::Plan{},pullback::Plan{}};
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
 // Pure epsilon rescalings are not offered under component demands (their
 // cost is invariant), so this one-master system screens exactly the cached
 // coordinate it installs. Installation stays bounded and after screening.
 check(selection.candidates_screened>=selection.full_candidates_checked&&selection.full_candidates_checked>=1&&
       selection.full_candidates_checked<=8,"installation was not deferred to the bounded improving frontier");
 check(selected.connection[0][0].is_zero(),"full-validation candidate differs from screened exact coordinate");
 // The upper operator is expressed at physical x=1-h, not reversed x=h.
 // V=1+x, so Enew=Eold*(2-h)=2. Compare complete retained jets, including zeros.
 auto h=graph.dimension.constant(original.endpoint_geometry(0).upper_overlap);
 auto v=exact_laurent_rows({{one+x}},one-h,3);
 auto composed=linear_boundary::compose(old.right,v,3);
 for(int k=0;k<=3;++k){const auto& lhs=selected.right.coefficients[0][0][k-selected.right.low];
   const auto& rhs=composed.coefficients[0][0][k-composed.low];
   check(acb_overlaps(lhs.raw(),rhs.raw()),"upper-cut coordinate pullback disagrees with rebuilt endpoint");
   B exact(k?0:2);check(acb_contains(lhs.raw(),exact.raw()),"requested physical upper functional changed");}
 // The selected basis reuses the original endpoint frames through the exact
 // pullback: selecting it must not build a second pair of endpoint series.
 check(search.statistics().endpoint_series_built==original.statistics().endpoint_series_built &&
       search.statistics().endpoint_series_built==2,"basis selection rebuilt exact endpoint series");
 check(search.statistics().basis_scans==1,"scan accounting");
 {
   // A smaller connection is worth adopting whatever the child window is:
   // without any epsilon loss, and above the closed-form leaf as well.
   auto lossless=graph;lossless.nodes[0].observable_rows={{one}};lossless.nodes[0].operations.resize(1);
   recursion::Evaluator plain(lossless,settings);auto adopted=plain.prepare_factored_consumers(0,0);
   check(plain.statistics().basis_scans==1 && plain.basis_selections().at(0).accepted &&
         plain.basis_selections().at(0).original_child_loss==0 && adopted.connection[0][0].is_zero() &&
         plain.statistics().endpoint_series_built==2,"smaller connection not adopted at zero child loss");
   auto above_leaf=graph;above_leaf.nodes.push_back(graph.nodes[0]);above_leaf.nodes[1].scalar_leaf=true;
   auto two=settings;two.causal_prescription->levels={{1},{1}};
   recursion::Evaluator leafward(above_leaf,two);(void)leafward.prepare_factored_consumers(0,0);
   check(leafward.statistics().basis_scans==1 && leafward.basis_selections().at(0).accepted,
         "smaller connection not adopted above the closed-form leaf");
 }
 // Earlier revisions ranked by demand scores and stored undecided scans as
 // "original retained"; the connection-size policy must not inherit them.
 check(recursion::basis_persistence::revision==std::string("physical-basis-connection-size-v5"),"selection revision did not invalidate earlier decisions");
 std::cout<<"Connection-size discovery, exact pullback installation, upper-cut identity and v5 revision passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

#include "diffexp/ft_spectral.hpp"
#include "diffexp/ft_spectral_checkpoint.hpp"
#include <iostream>
using namespace diffexp;
using B=Jet::Ball;
void require(bool b,const std::string& message){if(!b)throw std::runtime_error(message);}
void near(const B& value,const B& expected,const char* message){require(ft_spectral::sp::le(transport::magnitude(value-expected),B::from_strings("1e-20")),message);}
int main(){try {
 B::set_precision(256);ExactField field({"x","eps","I"});Exact x(field,"x"),eps(field,"eps"),imag(field,"I"),zero(field,0),one(field,1);
 ft_spectral::Options options;options.accuracy_goal=20;options.endpoint_clustering=true;options.seconds_budget=15;
 ft_spectral::Diagnostics stats;
 auto run=[&](const ExactEpsilonMatrix& a,const LaurentRows& initial,const ExactEpsilonMatrix& f,const std::vector<Exact>& path){auto result=ft_spectral::try_transport(a,initial,f,path,options,stats);require(bool(result),stats.reason);return *result;};
 // h=x: the transformed source is 1/x, not x. g'=1+g/x.
 auto logarithm=run({{-one/x}},{0,0,{{{B(0)}}}},{{one}},{one,one.constant(2)});
 B logtwo;acb_log(logtwo.raw(),B(2).raw(),256);
 near(logarithm.coefficients[0][0][0],B(2)*logtwo,"adjoint gauge/source sign");
 require(stats.normalized_diagonals==1,"rational diagonal was not normalized");
 // A non-logarithmic diagonal must retain its scalar collocation operator.
 auto exponential=run({{-one}},{0,0,{{{B(1)}}}},{{zero}},{zero,one});B expone;acb_exp(expone.raw(),B(1).raw(),256);
 near(exponential.coefficients[0][0][0],expone,"ungauged scalar exponential");
 require(stats.normalized_diagonals==0,"constant diagonal incorrectly gauged");
 // Negative Laurent boundary and positive epsilon edge: g1'=eps*g0.
 LaurentRows laurent{-2,1,{{{B(1),B(0),B(0),B(0)},{B(0),B(0),B(0),B(0)}}}};
 auto shifted=run({{zero,-eps},{zero,zero}},laurent,{{zero,zero}},{zero,one});
 require(shifted.low==-2&&shifted.high==1,"Laurent output window");
 for(unsigned k=0;k<4;++k){near(shifted.coefficients[0][0][k],B(k==0),"Laurent initial coefficient");near(shifted.coefficients[0][1][k],B(k==1),"positive epsilon Laurent shift");}
 // Source expands the lower bound independently of the initial window.
 auto sourced=run({{zero}},{0,1,{{{B(0),B(0)}}}},{{one/eps.pow(2)}},{zero,one});
 require(sourced.low==-2&&sourced.high==1,"Laurent forcing lower bound");
 for(unsigned k=0;k<4;++k)near(sourced.coefficients[0][0][k],B(k==0),"negative Laurent forcing");
 // DAG at epsilon zero, reverse feedback at epsilon one. Both SCC order and
 // previous epsilon-layer dependency are needed: cosh(sqrt(eps)*x), sinh()/sqrt(eps).
 LaurentRows feedback{0,3,{{{B(1),B(0),B(0),B(0)},{B(0),B(0),B(0),B(0)}},{{B(0),B(0),B(0),B(0)},{B(1),B(0),B(0),B(0)}}}};
 auto coupled=run({{zero,-one},{-eps,zero}},feedback,{{zero,zero},{zero,zero}},{zero,one});
 long even[]={1,2,24,720},odd[]={1,6,120,5040};
 for(unsigned k=0;k<4;++k){near(coupled.coefficients[0][0][k],B::from_strings("1/"+std::to_string(even[k])),"DAG feedback even coefficient");near(coupled.coefficients[0][1][k],B::from_strings("1/"+std::to_string(odd[k])),"DAG feedback odd coefficient");near(coupled.coefficients[1][1][k],B::from_strings("1/"+std::to_string(even[k])),"shared observable operator");near(coupled.coefficients[1][0][k],k?B::from_strings("1/"+std::to_string(odd[k-1])):B(0),"reverse observable feedback");}
 // Endpoint clustering changes only the parameter, including its Jacobian.
 // Sources exercise that Jacobian even when the diagonal is gauged away.
 auto delta=one/one.constant(50);B delta_ball=B::from_strings("1/50"),big=B(1)+delta_ball;
 B log_ratio;acb_log(log_ratio.raw(),(big/delta_ball).raw(),256);
 B root_product;acb_sqrt(root_product.raw(),(big*delta_ball).raw(),256);
 for(bool right:{false,true}) {
  auto pole=right?one+delta:-delta;
  auto rational_cluster=run({{-one/(x-pole)}},{0,0,{{{B(0)}}}},{{one}},{zero,one});
  near(rational_cluster.coefficients[0][0][0],(right?delta_ball:big)*log_ratio,"clustered rational gauge/source Jacobian");
  require(stats.clustered_legs==1&&stats.normalized_diagonals==1,"near-endpoint rational leg was not clustered and gauged");
  auto half_cluster=run({{-one/(one.constant(2)*(x-pole))}},{0,0,{{{B(0)}}}},{{one}},{zero,one});
  near(half_cluster.coefficients[0][0][0],B(2)*(right?root_product-delta_ball:big-root_product),"clustered half-integer gauge/source Jacobian or branch");
  require(stats.clustered_legs==1&&stats.normalized_diagonals==1,"near-endpoint half-integer leg was not clustered and gauged");
 }
 // The compiled sampler must divide complete epsilon jets, not specialize
 // the denominator at epsilon zero. Also exercise exact imaginary coefficients.
 auto denominator_source=run({{zero}},{0,3,{{{B(0),B(0),B(0),B(0)}}}},{{(one+imag*x)/(one+eps*x)}},{zero,one});
 for(unsigned k=0;k<4;++k) {
  B expected=B::from_strings("1/"+std::to_string(k+1)),imaginary=B::from_strings("1/"+std::to_string(k+2));
  acb_mul_onei(imaginary.raw(),imaginary.raw());expected+=imaginary;if(k%2)expected=-expected;
  near(denominator_source.coefficients[0][0][k],expected,"compiled epsilon-dependent denominator/source coefficients");
 }
 // A genuine two-component epsilon-zero SCC exercises the coupled inverse.
 auto block=run({{zero,-one},{-one,zero}},{0,0,{{{B(1)},{B(0)}}}},{{zero,zero}},{zero,one});
 near(block.coefficients[0][0][0],(expone+B(1)/expone)/B(2),"coupled SCC cosh");
 near(block.coefficients[0][1][0],(expone-B(1)/expone)/B(2),"coupled SCC sinh");
 require(std::find(stats.block_sizes.begin(),stats.block_sizes.end(),2)!=stats.block_sizes.end(),"coupled SCC not detected");
 // A full contour produces square-root monodromy despite resetting normalized
 // gauges at each leg. Straight chords stay away from the singular origin.
 auto monodromy=run({{-one/(x.constant(2)*x)}},{0,0,{{{B(1)}}}},{{zero}},{one,imag,-one,-imag,one});
 near(monodromy.coefficients[0][0][0],B(-1),"half-integer gauge branch continuation");require(stats.legs==4,"contour legs");
 for(const auto& pole:{one/one.constant(3),one}) {
  auto rejected=ft_spectral::try_transport({{-one/(x-pole)}},{0,0,{{{B(1)}}}},{{zero}},{zero,one},options,stats);
  require(!rejected,"canceled integer-gauge singularity accepted");
 }
 // Input uncertainty is propagated, not replaced by agreement of midpoints.
 LaurentRows uncertain{0,0,{{{B(1)}}}};arb_add_error_2exp_si(acb_realref(uncertain.coefficients[0][0][0].raw()),-30);
 auto original=uncertain.coefficients[0][0][0];auto rejected=ft_spectral::try_transport({{zero}},uncertain,{{zero}},{zero,one},options,stats);
 require(!rejected,"large input ball accepted at high accuracy");require(acb_equal(original.raw(),uncertain.coefficients[0][0][0].raw()),"fallback input mutated");
 auto constrained=options;constrained.max_nodes=16;constrained.accuracy_goal=40;
 require(!ft_spectral::try_transport({{-one}},{0,0,{{{B(1)}}}},{{zero}},{zero,one},constrained,stats),"underresolved exponential accepted");
 auto invalid=options;invalid.accuracy_goal=std::numeric_limits<unsigned>::max();
 require(!ft_spectral::try_transport({{zero}},{0,0,{{{B(1)}}}},{{zero}},{zero,one},invalid,stats),"overflowing accuracy goal accepted");
 // Squared source vanishes at every 8/12/16 Chebyshev node but has a
 // strictly positive integral. Resolution agreement alone must not accept it.
 auto u=[&](unsigned n){auto before=one,now=x.constant(2)*(x.constant(2)*x-one);if(!n)return before;for(unsigned k=1;k<n;++k){auto next=x.constant(2)*(x.constant(2)*x-one)*now-before;before=now;now=next;}return now;};
 auto aliased=(x*(x-one)*u(7)*u(11)*u(15)).pow(2);
 auto coarse=options;coarse.max_nodes=16;
 require(!ft_spectral::try_transport({{zero}},{0,0,{{{B(0)}}}},{{aliased}},{zero,one},coarse,stats),"aliased polynomial source accepted below degree");
 require(!ft_spectral::try_transport({{one/eps}},{0,0,{{{B(1)}}}},{{zero}},{zero,one},options,stats),"negative epsilon connection accepted");
 // At a fixed node cap, subdivision must recover the same exponential.
 auto split=options;split.max_nodes=16;split.accuracy_goal=20;split.max_subdivisions=64;split.seconds_budget=20;
 LaurentRows exp_initial{0,0,{{{B(1)}}}};
 auto unsplit=split;unsplit.max_subdivisions=0;
 require(!ft_spectral::try_transport({{-one}},exp_initial,{{zero}},{zero,one},unsplit,stats),"coarse unsplit exponential unexpectedly resolved");
 auto divided=ft_spectral::try_transport({{-one}},exp_initial,{{zero}},{zero,one},split,stats);
 require(bool(divided),"subdivision exponential: "+stats.reason);
 near(divided->coefficients[0][0][0],expone,"subdivided exponential oracle");
 require(stats.subdivisions>0&&stats.accepted_subsegments==stats.subdivisions+1&&stats.legs==1,"subdivision diagnostics");
 auto key=ft_spectral_checkpoint::identity({{-one}},exp_initial,{{zero}},{zero,one},split);
 require(key!=ft_spectral_checkpoint::identity({{-one}},exp_initial,{{zero}},{zero,one},unsplit),"subdivision checkpoint identity collision");
 auto payload=ft_spectral_checkpoint::detail::encode(*divided,key,stats);
 boost::json::object envelope{{"payload",payload},{"sha256",artifacts::detail::sha256(artifacts::detail::canonical(payload))}};
 ft_spectral::Diagnostics restored;
 auto restored_rows=ft_spectral_checkpoint::detail::decode(boost::json::parse(boost::json::serialize(envelope)),key,exp_initial,0,1,split,restored);
 require(restored.subdivisions==stats.subdivisions&&restored.accepted_subsegments==stats.accepted_subsegments,"subdivision checkpoint diagnostics");
 near(restored_rows.coefficients[0][0][0],expone,"subdivision checkpoint result");
 // A near-interior complex pole is not cured by endpoint clustering.
 auto near_pole=one/one.constant(2)+imag/one.constant(64);
 auto nearby=split;nearby.max_nodes=32;nearby.endpoint_clustering=true;
 auto nearby_unsplit=nearby;nearby_unsplit.max_subdivisions=0;
 LaurentRows zero_initial{0,0,{{{B(0)}}}};
 require(!ft_spectral::try_transport({{zero}},zero_initial,{{one/(x-near_pole)}},{zero,one},nearby_unsplit,stats),"near-interior pole unexpectedly resolved without splitting");
 auto near_result=ft_spectral::try_transport({{zero}},zero_initial,{{one/(x-near_pole)}},{zero,one},nearby,stats);
 require(bool(near_result),"near-interior subdivision: "+stats.reason);
 B pball=B::from_strings("1/2","1/64"),logleft,logright;
 acb_log(logleft.raw(),(-pball).raw(),256);acb_log(logright.raw(),(B(1)-pball).raw(),256);
 near(near_result->coefficients[0][0][0],logright-logleft,"near-interior logarithmic integral");
 require(stats.subdivisions>0,"near-interior subdivision missing");
 auto bounded=split;bounded.max_subdivisions=1;bounded.accuracy_goal=40;
 require(!ft_spectral::try_transport({{-one}},exp_initial,{{zero}},{zero,one},bounded,stats)&&stats.subdivisions<=1,"subdivision count bound ignored");
 // Existing uncertainty between local and whole-arm budgets must survive.
 auto inherited_options=options;inherited_options.accuracy_goal=40;inherited_options.max_subdivisions=64;
 LaurentRows inherited{0,0,{{{B(1)}}}};
 arb_add_error_2exp_si(acb_realref(inherited.coefficients[0][0][0].raw()),-146);
 const auto inherited_before=inherited.coefficients[0][0][0];
 auto inherited_result=ft_spectral::try_transport({{zero}},inherited,{{zero}},{zero,one},inherited_options,stats);
 require(bool(inherited_result),"inherited sub-global uncertainty rejected: "+stats.reason);
 require(arb_contains(acb_realref(inherited_result->coefficients[0][0][0].raw()),acb_realref(inherited_before.raw())),"inherited uncertainty stripped");
 require(acb_equal(inherited_before.raw(),inherited.coefficients[0][0][0].raw()),"input radii modified");
 arb_add_error_2exp_si(acb_realref(inherited.coefficients[0][0][0].raw()),-130);
 inherited_options.max_subdivisions=1;
 require(!ft_spectral::try_transport({{zero}},inherited,{{zero}},{zero,one},inherited_options,stats),"above-global inherited uncertainty accepted");
 auto zero_goal=split;zero_goal.accuracy_goal=0;
 require(!ft_spectral::try_transport({{zero}},exp_initial,{{zero}},{zero,one},zero_goal,stats),"subdivision accepted zero accuracy goal");
 auto contour_split=ft_spectral::try_transport({{-one/(x.constant(2)*x)}},exp_initial,{{zero}},{one,imag,-one,-imag,one},split,stats);
 require(bool(contour_split),stats.reason);near(contour_split->coefficients[0][0][0],B(-1),"adaptive original-leg branch continuation");
 require(stats.legs==4&&stats.accepted_subsegments==4+stats.subdivisions,"adaptive original-leg count");
 auto timed=split;timed.seconds_budget=1e-12;
 require(!ft_spectral::try_transport({{-one}},exp_initial,{{zero}},{zero,one},timed,stats),"total subdivision time bound ignored");
 require(!ft_spectral::try_transport({{-one/(x-one/one.constant(3))}},exp_initial,{{zero}},{zero,one},split,stats)&&stats.subdivisions==0,"subdivision hid a true path pole");
 require(!ft_spectral::try_transport({{one/eps}},exp_initial,{{zero}},{zero,one},split,stats)&&stats.subdivisions==0,"subdivision retried structural unsupported input");
 std::cout<<"FT spectral gauges, SCC ordering, Laurent windows, branches and fallback passed\n";
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}

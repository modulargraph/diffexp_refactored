#include "diffexp/ft_spectral.hpp"
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
 std::cout<<"FT spectral gauges, SCC ordering, Laurent windows, branches and fallback passed\n";
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}

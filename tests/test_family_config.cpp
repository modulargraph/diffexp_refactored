#include "diffexp/family_config.hpp"
#include "diffexp/banana_oracle.hpp"
#include <iostream>
using namespace diffexp;
int main(int argc,char** argv) {
 try {
  for(const auto& name:feynman::example_names()) {
    auto example=feynman::example_family(name);auto config=family_config::parse(family_config::describe(example));
    if(config.family.name!=name || config.family.momenta.lines.size()!=example.momenta.lines.size())throw std::runtime_error("family template roundtrip");
  }
  auto fixed_geometry=feynman::example_family("henn_double_pentagon_x0");auto renamed=fixed_geometry;renamed.name="different-label";
  if(!family_config::same_geometry(fixed_geometry,renamed))throw std::runtime_error("geometry comparison depends on a label");
  renamed.momenta.lines[0].mass_squared=Rational(1);
  if(family_config::same_geometry(fixed_geometry,renamed))throw std::runtime_error("fixed basis accepted changed geometry");
  auto document=family_config::describe(feynman::example_family("bubble"));auto& o=document.as_object();
  o["name"]="a-family-never-registered";o["external_gram"]=boost::json::array{boost::json::array{"-3"}};
  for(auto& line:o.at("propagators").as_array())line.as_object()["mass_squared"]="2";
  o["epsilon_order"]=0;o["integrals"]=boost::json::array{boost::json::array{1,1}};
  auto config=family_config::parse(document);
  auto graph=recursion::prepare(config.family,config.integrals,config.preparation);
  auto result=recursion::Evaluator(graph,config.numerical).evaluate(0);
  using B=kernel::ComplexBall;B::set_precision(384);B root,argument,atanh;
  acb_sqrt(root.raw(),B(33).raw(),384);acb_sqrt(argument.raw(),(B(3)/B(11)).raw(),384);acb_atanh(atanh.raw(),argument.raw(),384);
  auto analytic=B(4)*atanh/root;
  if(!acb_overlaps(result.values[0][0].raw(),analytic.raw()) || !result.taylor_tail_certified)throw std::runtime_error("unregistered family analytic enclosure");
  // Even a familiar special-case label cannot select a different ladder.
  auto labelled=config;labelled.family.name="henn_double_pentagon_x0";
  auto labelled_graph=recursion::prepare(labelled.family,labelled.integrals,labelled.preparation);
  auto labelled_value=recursion::Evaluator(labelled_graph,labelled.numerical).evaluate(0);
  if(!acb_overlaps(labelled_value.values[0][0].raw(),analytic.raw()))throw std::runtime_error("family name changed the configured integral");
  auto bad=document;bad.as_object()["numerical"].as_object()["ordnary_order"]=40;
  bool rejected=false;try{family_config::parse(bad);}catch(const std::exception&){rejected=true;}if(!rejected)throw std::runtime_error("unknown field accepted");
  bad=document;bad.as_object()["propagators"].as_array()[0].as_object()["mass_squared"]=0.5;
  rejected=false;try{family_config::parse(bad);}catch(const std::exception&){rejected=true;}if(!rejected)throw std::runtime_error("approximate family coefficient accepted");
  if(argc>1) {
    auto sunrise=family_config::describe(feynman::example_family("sunrise"));auto& definition=sunrise.as_object();
    definition["name"]="configured-three-distinct-masses";definition["external_gram"]=boost::json::array{boost::json::array{"-3"}};definition["epsilon_order"]=0;
    const char* masses[]={"2","3/2","1"};for(unsigned i=0;i<3;++i)definition["propagators"].as_array()[i].as_object()["mass_squared"]=masses[i];
    auto custom=family_config::parse(sunrise);custom.preparation.reduction.provider.executable=argv[1];
    auto begin=std::chrono::steady_clock::now();auto prepared=recursion::prepare(custom.family,custom.integrals,custom.preparation);
    auto prepared_at=std::chrono::steady_clock::now();auto computed=recursion::Evaluator(prepared,custom.numerical).evaluate(0);auto numerical_at=std::chrono::steady_clock::now();
    oracle::BananaOptions limits;limits.target_bits=72;limits.working_bits=192;
    auto reference=oracle::banana_bessel({Rational("2/3"),Rational("1/2"),Rational("1/3")},limits).value/B(3);
    // At D=2 the simultaneous rescaling p²,m² by three gives a factor 1/3.
    if(!acb_overlaps(computed.values[0][-computed.low].raw(),reference.raw()))throw std::runtime_error("configured two-loop family disagrees with independent Bessel integral");
    std::cout<<"Configured two-loop family: preparation_seconds="<<std::chrono::duration<double>(prepared_at-begin).count()
      <<" numerical_seconds="<<std::chrono::duration<double>(numerical_at-prepared_at).count()<<" reference_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-numerical_at).count()<<'\n';
  }
  std::cout<<"Generic family roundtrips and changed-kinematics analytic enclosure passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}

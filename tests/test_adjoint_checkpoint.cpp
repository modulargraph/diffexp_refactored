#include "diffexp/adjoint_checkpoint.hpp"
#include <iostream>
#include <limits>
#include <tuple>
using namespace diffexp;using B=Jet::Ball;namespace cp=adjoint_checkpoint;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* why){try{f();}catch(const std::exception&){return;}throw std::runtime_error(why);}
void equal(const LaurentRows& a,const LaurentRows& b) {
  require(a.low==b.low && a.high==b.high && a.columns()==b.columns() && a.coefficients.size()==b.coefficients.size(),"continuation changed shape");
  for(unsigned i=0;i<a.coefficients.size();++i)for(unsigned j=0;j<a.columns();++j)for(unsigned k=0;k<a.coefficients[i][j].size();++k)
    require(acb_equal(a.coefficients[i][j][k].raw(),b.coefficients[i][j][k].raw()),"continuation changed retained ball bits");
}
int main(){try {
  B::set_precision(384);ExactField field({"x","eps","I"});Exact zero(field,"0"),one(field,"1"),x(field,"x"),e(field,"eps"),two(field,"2");
  const ExactEpsilonMatrix a{{one/(x+two)}},forcing{{one/(e*(x+two))},{-one/(e*(x+two))}};
  B seed=B::from_strings("7/13","2/17");arb_add_error_2exp_si(acb_realref(seed.raw()),-240);
  LaurentRows initial{0,2,{{{seed,B(0),B(0)}},{{-seed,B(0),B(0)}}}};
  const std::vector<Exact> path{zero,Exact(field,"1/2+I/3"),Exact(field,"1+I/3"),two};
  using Chart=std::tuple<unsigned,double,double>;
  std::vector<Chart> sequence;std::vector<AdjointContinuation> snapshots;
  AdjointOptions options;
  options.chart_observer=[&](unsigned leg,double from,double to,const LaurentRows&){sequence.emplace_back(leg,from,to);};
  options.continuation_observer=[&](const AdjointContinuation& state){snapshots.push_back(state);};
  const auto continuous=transport_adjoint_rows(a,initial,forcing,path,options);
  require(snapshots.size()>2 && continuous.low==-1,"forced multi-leg continuation fixture");
  for(unsigned row=0;row<2;++row)for(int k=-1;k<=2;++k) {
    B expected=k==-1?B::from_strings("1/2"):k==0?seed/B(2):B(0);if(row)expected=-expected;
    require(NativeTailMagnitude::upper_abs(continuous.coefficients[row][0][k+1]-expected).approximate_upper()<1e-40,
      "forced Laurent continuation disagrees with independent analytic solution");
  }
  for(const std::size_t stop:{std::size_t(0),snapshots.size()/2,snapshots.size()-1}) {
    std::vector<Chart> resumed_sequence;AdjointConditioningStats work;
    options=AdjointOptions{};options.conditioning_stats=&work;
    options.continuation=std::make_shared<const AdjointContinuation>(snapshots[stop]);
    options.chart_observer=[&](unsigned leg,double from,double to,const LaurentRows&){resumed_sequence.emplace_back(leg,from,to);};
    equal(continuous,transport_adjoint_rows(a,initial,forcing,path,options));
    require(resumed_sequence==std::vector<Chart>(sequence.begin()+stop+1,sequence.end()),"continuation changed original future chart parameters");
    if(stop+1==snapshots.size())require(work.polynomial_charts==0 && work.rational_compilations==0,"completed continuation repeated recurrence work");
  }
  const auto invalid=[&](auto change) {
    auto state=snapshots.front();change(state);AdjointOptions opt;opt.continuation=std::make_shared<const AdjointContinuation>(state);
    rejects([&]{transport_adjoint_rows(a,initial,forcing,path,opt);},"invalid continuation accepted");
  };
  invalid([](auto& s){s.leg=100;});invalid([](auto& s){s.parameter=std::numeric_limits<double>::quiet_NaN();});
  invalid([](auto& s){s.parameter=-0.1;});invalid([](auto& s){s.parameter=1.1;});
  invalid([](auto& s){s.accepted_charts=0;});invalid([](auto& s){s.accepted_charts=20001;});
  invalid([](auto& s){s.rows.low=0;});invalid([](auto& s){s.rows.coefficients.pop_back();});
  invalid([](auto& s){acb_indeterminate(s.rows.coefficients[0][0][0].raw());});
  auto exhausted=snapshots.front();exhausted.parameter=0.1;exhausted.accepted_charts=20000;
  options=AdjointOptions{};options.continuation=std::make_shared<const AdjointContinuation>(exhausted);
  rejects([&]{transport_adjoint_rows(a,initial,forcing,path,options);},"resume reset finite chart budget");

  // Exercise reconstruction of the lazily compiled homogeneous maps too,
  // using inherited uncertainty large enough to enter centered arithmetic.
  B first(1),second(2);arb_add_error_2exp_si(acb_realref(first.raw()),-300);
  arb_add_error_2exp_si(acb_realref(second.raw()),-300);
  LaurentRows wide{0,0,{{{first},{second}},{{second},{first}}}};
  const auto large=one.constant(1000);
  const ExactEpsilonMatrix nilpotent{{-large,-large},{large,large}},unforced{{zero,zero},{zero,zero}};
  for(bool polynomial:{true,false}) {
    AdjointContinuation restart;AdjointConditioningStats work;unsigned count=0;
    options=AdjointOptions{};options.polynomial_recurrence=polynomial;options.max_rows_per_batch=1;
    options.conditioning_stats=&work;
    options.continuation_observer=[&](const AdjointContinuation& state){if(++count==1)restart=state;};
    const auto expected=transport_adjoint_rows(nilpotent,wide,unforced,path,options);
    require(work.centered_charts>0,"checkpoint test did not exercise centered arithmetic");
    options.continuation_observer={};options.continuation=std::make_shared<const AdjointContinuation>(restart);
    equal(expected,transport_adjoint_rows(nilpotent,wide,unforced,path,options));
  }

  const auto directory=std::filesystem::temp_directory_path()/("de3-ordinary-checkpoint-"+std::to_string(::getpid()));
  struct Cleanup {std::filesystem::path path;~Cleanup(){std::error_code ec;std::filesystem::remove_all(path,ec);}} cleanup{directory};
  cp::Statistics stats;cp::Options storage{directory,64*1024*1024,&stats};options=AdjointOptions{};
  const auto key=cp::identity(a,initial,forcing,path,options);const auto file=directory/(key+".json");
  auto different=initial;different.coefficients[0][0][0]+=B(1);
  require(key!=cp::identity(a,different,forcing,path,options),"cache identity omitted original input balls");
  auto alternate=options;alternate.taylor_order=81;
  require(key!=cp::identity(a,initial,forcing,path,alternate),"cache identity omitted order");
  alternate=options;alternate.max_charts_per_leg=19999;
  require(key!=cp::identity(a,initial,forcing,path,alternate),"cache identity omitted finite chart budget");
  require(key!=cp::identity({{two/(x+two)}},initial,forcing,path,options),"cache identity omitted connection");
  require(key!=cp::identity(a,initial,{{zero},{zero}},path,options),"cache identity omitted forcing");
  auto alternate_path=path;alternate_path[1]=Exact(field,"1/2+I/4");
  require(key!=cp::identity(a,initial,forcing,alternate_path,options),"cache identity omitted original contour");
  B::set_precision(320);require(key!=cp::identity(a,initial,forcing,path,options),"cache identity omitted precision");B::set_precision(384);
  options.chart_observer=[](unsigned,double,double,const LaurentRows&){throw std::runtime_error("intentional interruption");};
  rejects([&]{cp::transport(a,initial,forcing,path,options,storage);},"interruption was swallowed");
  require(stats.saved==1 && std::filesystem::exists(file),"accepted chart was not saved before user observer");
  const auto saved_envelope=cp::detail::read(file,storage.max_bytes);
  B::set_precision(64);auto loaded=cp::detail::decode(saved_envelope,key);B::set_precision(384);
  equal(loaded.rows,snapshots.front().rows);
  options=AdjointOptions{};equal(continuous,cp::transport(a,initial,forcing,path,options,storage));
  require(stats.loaded==1,"interrupted checkpoint was not reused");
  const auto saved_count=stats.saved;options.chart_observer=[](unsigned,double,double,const LaurentRows&){throw std::runtime_error("completed arm unexpectedly ran");};
  equal(continuous,cp::transport(a,initial,forcing,path,options,storage));
  require(stats.saved==saved_count && stats.completed_reused==1,"completed arm was recomputed or rewritten");
  auto corrupt=saved_envelope;corrupt.as_object().at("payload").as_object()["leg"]=42;
  cp::detail::publish(file,boost::json::serialize(corrupt));
  rejects([&]{cp::transport(a,initial,forcing,path,options,storage);},"corrupt matching checkpoint silently restarted");
  rejects([&]{cp::detail::decode(saved_envelope,std::string(64,'0'));},"checkpoint for different numerical problem accepted");
  rejects([&]{cp::detail::read(file,1);},"checkpoint file budget ignored");
  std::cout<<"Original-path continuation: bit-identical retained balls/charts, forced analytic solution, finite budgets, atomic interruption/reuse, lossless I/O and identity/corruption checks passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

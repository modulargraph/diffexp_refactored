#include "diffexp/henn_boundary.hpp"
#include <iostream>
using namespace diffexp;using B=kernel::ComplexBall;
void require(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* why){bool failed=false;try{f();}catch(const std::exception&){failed=true;}require(failed,why);}
henn::Integral index(int n){henn::Integral a(11,0);a[0]=n;return a;}
long weight(const henn::Integral& a){long v=19;for(int n:a)v=(v*17+n+1024)%101;return v+1;}
data::Expr assigned(data::Expr e) {
  if(e.head=="XB") {henn::Integral a;long sum=0;for(const auto& q:e.args){a.push_back(data::integer(q));sum+=a.back();}return data::Reader(std::to_string((sum%2?-1:1)*weight(a))).read();}
  for(auto& child:e.args)child=assigned(std::move(child));return e;
}
void small(const B& value,const char* why){mag_t bound;mag_init(bound);acb_get_mag(bound,value.raw());bool ok=value.contains_zero() && mag_get_d(bound)<1e-55;mag_clear(bound);require(ok,why);}
int main(int argc,char** argv){try {
  B::set_precision(256);henn::CanonicalBasis basis;Exact one(basis.field,1),eps(basis.field,"eps");
  basis.components={henn::Observable{{index(1),henn::Coefficient(one)}}};
  require(henn::needed_scalar_high(basis,4)==0,"Henn finite scalar demand");
  B gamma;arb_const_euler(acb_realref(gamma.raw()),256);Jet e(0,5,256);e.set(1,B(1));
  Jet minus(0,5,256);minus.set(1,-B(2)*gamma);auto stripped=minus.exp();
  LaurentBoundary raw{-4,{{}},true};for(unsigned k=0;k<5;++k)raw.values[0].push_back(stripped.at(k));
  auto canonical=henn::project_boundary(basis,{index(1)},raw,4,256);
  require(canonical.low==0 && canonical.high()==4 && canonical.taylor_tail_certified,"normalization window/certificate");
  require(acb_contains(canonical.values[0][0].raw(),B(1).raw()),"global eps^4 normalization");
  for(unsigned k=1;k<5;++k)small(canonical.values[0][k],"global exp(+2 EulerGamma eps) convention");
  raw.taylor_tail_certified=false;require(!henn::project_boundary(basis,{index(1)},raw).taylor_tail_certified,"invented source tail certificate");
  raw.values[0].pop_back();bool demand=false;
  try{henn::project_boundary(basis,{index(1)},raw);}catch(const BoundaryDemand& d){demand=d.required_high==0;}
  require(demand,"missing upper raw coefficient padded with zero");
  raw.values[0].push_back(B(0));
  basis.components[0].begin()->second=henn::Coefficient(Exact(basis.field,"1/(eps*(1-eps))"));
  require(henn::needed_scalar_high(basis,4)==1,"epsilon-polar canonical coefficient lookahead");
  demand=false;try{henn::project_boundary(basis,{index(1)},raw);}catch(const BoundaryDemand& d){demand=d.required_high==1;}
  require(demand,"epsilon pole did not request additional scalar order");
  raw.values[0].push_back(B(3));auto poles=henn::project_boundary(basis,{index(1)},raw);
  require(poles.low==-1,"forbidden canonical negative order discarded");
  auto audit=henn::audit_negative_poles(poles,Rational("1/1000000"));
  require(!audit.pass && audit.failures.size()==1 && audit.failures[0].epsilon_order==-1 && audit.failures[0].component==1,"negative canonical pole audit");
  LaurentBoundary geometric_raw{-4,{std::vector<B>(6)},false};geometric_raw.values[0][0]=B(1);
  auto geometric=henn::project_boundary(basis,{index(1)},geometric_raw,4,256);
  Jet epsilon6(0,6,256);epsilon6.set(1,B(1));Jet exponential6(0,6,256);exponential6.set(1,B(2)*gamma);
  auto direct_geometric=exponential6.exp()/(epsilon6.constant(1)-epsilon6);
  for(unsigned k=0;k<6;++k)small(geometric.values[0][k]-direct_geometric.at(k),"future rational-epsilon coefficient convolution");
  LaurentBoundary uncertain{-1,{{B(0),B(1)}},false};arb_add_error_2exp_si(acb_realref(uncertain.values[0][0].raw()),0);
  require(!henn::audit_negative_poles(uncertain,Rational("1/100")).pass,"zero-containing broad pole ball incorrectly accepted");
  uncertain.values[0][0]=B::from_strings("1/1000000000");require(henn::audit_negative_poles(uncertain,Rational("1/1000000")).pass,"small explicit negative residue rejected");
  basis.components[0].begin()->second=henn::Coefficient(one.constant(0),one);
  LaurentBoundary odd_raw{-4,{{B(1),B(0),B(0),B(0),B(0)}},true};auto odd=henn::project_boundary(basis,{index(1)},odd_raw);
  B sqrt3;arb_set_ui(acb_imagref(sqrt3.raw()),3);arb_sqrt(acb_imagref(sqrt3.raw()),acb_imagref(sqrt3.raw()),256);
  small(odd.values[0][0]-sqrt3,"positive i sqrt3 embedding");
  basis.components[0].begin()->second=henn::Coefficient(Exact(basis.field,"(eps+1)/eps-1/eps"));
  B shared(1);arb_add_error_2exp_si(acb_realref(shared.raw()),-20);odd_raw.values[0][0]=shared;
  auto cancellation=henn::project_boundary(basis,{index(1)},odd_raw,0);
  require(acb_equal(cancellation.values[0][0].raw(),shared.raw()),"exact canonical row cancellation inflated source uncertainty");
  auto structurally_zero=henn::project_boundary(basis,{index(1)},LaurentBoundary{0,{{B(7)}},true},4);
  for(unsigned k=0;k<4;++k)require(structurally_zero.values[0][k].is_zero(),"below-declared-low source values are not structural zeros");
  require(acb_equal_si(structurally_zero.values[0][4].raw(),7),"shifted regular source coefficient");
  rejects([&]{henn::project_boundary(basis,{index(2)},odd_raw);},"missing scalar target accepted");
  rejects([&]{henn::project_boundary(basis,{index(1),index(1)},LaurentBoundary{-4,{odd_raw.values[0],odd_raw.values[0]},true});},"duplicate source targets accepted");
  if(argc>1) {
    auto all=henn::read_x0(argv[1]);auto original=data::read_file(argv[1]);
    require(henn::needed_scalar_high(all,4)==0,"actual all108 require extra positive scalar orders");
    LaurentBoundary inputs{-4,Boundary(all.scalar_targets.size(),std::vector<B>(5)),false};
    for(std::size_t i=0;i<all.scalar_targets.size();++i)inputs.values[i][0]=B(weight(all.scalar_targets[i]));
    auto projected=henn::project_boundary(all,all.scalar_targets,inputs,4,256);
    require(projected.values.size()==108 && projected.low==0 && !projected.taylor_tail_certified,"all108 projected shape");
    Jet context(0,5,256);auto root=evaluate(data::Reader("I*Sqrt[3]").read(),context,{});
    Jet normalization(0,5,256);normalization.set(1,B(2)*gamma);normalization=normalization.exp();
    std::map<std::string,Jet> vars{{"s12",context.constant(3)},{"s23",context.constant(-1)},
      {"s34",context.constant(1)},{"s45",context.constant(1)},{"s15",context.constant(-1)},{"eps5",root}};
    for(unsigned i=0;i<108;++i) {
      auto direct=evaluate(assigned(original.args[i]),context,vars)*normalization;
      for(unsigned k=0;k<5;++k)small(projected.values[i][k]-direct.at(k),"all108 original AST normalization comparison");
    }
    auto reversed_targets=all.scalar_targets;auto reversed_inputs=inputs;
    std::reverse(reversed_targets.begin(),reversed_targets.end());std::reverse(reversed_inputs.values.begin(),reversed_inputs.values.end());
    auto reordered=henn::project_boundary(all,reversed_targets,reversed_inputs,4,256);
    for(unsigned i=0;i<108;++i)for(unsigned k=0;k<5;++k)
      require(acb_equal(reordered.values[i][k].raw(),projected.values[i][k].raw()),"ordered scalar-source remapping failed");
    require(henn::audit_negative_poles(projected,Rational("1/1000000")).pass,"pole-free synthetic all108 failed audit");
    std::cout<<"all108 projected epsilon0..4; raw scalar high0 verified\n";
  }
  std::cout<<"Henn canonical Laurent projection PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}

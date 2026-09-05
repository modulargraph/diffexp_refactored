#pragma once
#include "diffexp/canonical.hpp"
#include "diffexp/geometry.hpp"
#include <array>

namespace diffexp {
inline int run_planar(const std::string& directory,const std::string& family) {
  using B=Jet::Ball; B::set_precision(256);
  if(family!="1loop" && family!="zmz" && family!="mzz" && family!="zzz")
    throw std::invalid_argument("unknown planar family");
  const auto begin=std::chrono::steady_clock::now();
  auto tensor=data::read_file(directory+"/"+family+"/diffEq-"+family+".m");
  auto samples=data::read_file(directory+"/"+family+"/numIntegrals-"+family+".m");
  auto alphabet=data::read_file(directory+"/alphabet.m");
  if(tensor.head!="SparseArray" || tensor.args.size()!=4)throw std::invalid_argument("expected published SparseArray tensor");
  const auto& shape=tensor.args.at(1).args;
  unsigned letters=data::integer(shape.at(0)),d=data::integer(shape.at(1));
  if(data::integer(shape.at(2))!=d || alphabet.args.size()!=letters)throw std::invalid_argument("planar tensor shape mismatch");
  ExactField field({"t","I","r1","r2","r3"});Exact t(field,"t"),imag(field,"I");
  const auto& storage=tensor.args.at(3).args;
  const auto& offsets=storage.at(1).args.at(0).args;
  const auto& coordinates=storage.at(1).args.at(1).args;
  const auto& values=storage.at(2).args;
  if(offsets.size()!=letters+1 || coordinates.size()!=values.size())throw std::invalid_argument("malformed compressed tensor");
  std::vector<CanonicalEntry> entries;std::set<unsigned> active_set;
  for(unsigned l=0;l<letters;++l)
    for(long n=data::integer(offsets[l]);n<data::integer(offsets[l+1]);++n) {
      const auto& ij=coordinates.at(n).args;
      auto c=evaluate_exact(values.at(n),t,{}).rational();if(c.is_zero())continue;
      entries.push_back({static_cast<unsigned>(data::integer(ij.at(0))-1),static_cast<unsigned>(data::integer(ij.at(1))-1),l,c});active_set.insert(l);
    }
  std::vector<unsigned> active(active_set.begin(),active_set.end());
  for(auto& e:entries)e.letter=std::lower_bound(active.begin(),active.end(),e.letter)-active.begin();
  const auto& first=samples.args.at(0).args;
  const auto& last=samples.args.at(5).args;
  const std::array<std::string,6> names{"s12","s15","s23","s34","s45","p1s"};
  auto point=[&](const data::Expr& e) {
    std::map<std::string,Exact> out;
    for(const auto& rule:e.args) {
      if(rule.head!="Rule" || rule.args.size()!=2)throw std::invalid_argument("invalid phase-space point");
      out.emplace(rule.args[0].head,evaluate_exact(rule.args[1],t,{}));
    }
    return out;
  };
  auto initial=point(first.at(0)),final=point(last.at(0));
  auto boundary=[&](const data::Expr& e) {
    if(e.args.size()!=5)throw std::invalid_argument("planar boundary weight count");
    Boundary out(d,std::vector<B>(5,B(0)));Jet scalar(0,1,256);
    for(unsigned k=0;k<5;++k) {
      if(e.args[k].args.size()!=d)throw std::invalid_argument("planar boundary dimension");
      for(unsigned i=0;i<d;++i)out[i][k]=evaluate(e.args[k].args[i],scalar,{}).at(0);
    }
    return out;
  };
  auto current=boundary(first.at(1)),reference=boundary(last.at(1));
  std::map<std::string,Exact> path;
  for(unsigned i=0;i<names.size();++i) {
    auto a=initial.at(names[i]),b=final.at(names[i]);
    auto direction=t.constant(Rational(10+i)/Rational(10));
    path.emplace(names[i],a+(b-a)*t+imag*t.constant(Rational("2/5"))*direction*t*(t.constant(1)-t));
  }
  const std::array<std::string,3> root_names{"tr5","sqrtG3","sqrtG3nc"};
  const std::array<data::Expr,3> radicands{
    data::Reader("(-s12*s15+s12*s23+p1s*s34+s15*s45-s34*s45-s23*s34)^2-4*s23*s34*s45*(p1s-s12-s15+s34)").read(),
    data::Reader("p1s^2+(s23-s45)^2-2*p1s*(s23+s45)").read(),
    data::Reader("(s12+s15)^2-4*p1s*s34").read()};
  std::vector<Exact> root_squares;
  for(unsigned i=0;i<3;++i) {
    root_squares.push_back(reduce_square(evaluate_exact(radicands[i],t,path),1,t.constant(-1)));
    path.emplace(root_names[i],t.variable(i+2));
  }
  std::set<std::string> unique_polynomials;std::vector<B> roots;
  auto add_roots=[&](Exact polynomial) {
    polynomial=reduce_square(polynomial,1,t.constant(-1));
    for(unsigned j=0;j<3;++j) {
      polynomial=polynomial_norm(polynomial,j+2,root_squares[j]);
      polynomial=reduce_square(polynomial,1,t.constant(-1));
    }
    polynomial=polynomial_norm(polynomial,1,t.constant(-1));
    if(polynomial.is_zero())throw std::domain_error("degenerate planar norm polynomial");
    if(polynomial.is_rational())return;
    polynomial=polynomial/polynomial.constant(polynomial.numerator_terms()[0].coefficient);
    if(!unique_polynomials.insert(polynomial.str()).second)return;
    auto found=polynomial_roots(polynomial,0);roots.insert(roots.end(),found.begin(),found.end());
  };
  std::set<unsigned> active_roots;
  for(auto l:active) {
    auto letter=evaluate_exact(alphabet.args[l],t,path);
    for(const auto& terms:{letter.numerator_terms(),letter.denominator_terms()})
      for(const auto& term:terms)for(unsigned j=0;j<3;++j)if(term.powers[j+2])active_roots.insert(j);
    add_roots(letter.numerator());add_roots(letter.denominator());
  }
  for(auto j:active_roots)add_roots(root_squares[j]);
  std::cout<<"Planar "<<family<<": "<<d<<" masters, "<<active.size()<<" active letters, "<<roots.size()<<" isolated singularity candidates"<<std::endl;
  constexpr unsigned order=40;
  unsigned charts=0;double center=0;
  std::vector<B> continued_roots;
  while(center<1) {
    if(++charts>10000)throw std::runtime_error("planar chart budget exhausted");
    double next=clearance_endpoint(center,roots);
    // Both center and endpoint are binary rationals. Use their exact difference
    // as the local evaluation point to avoid gaps from rounded addition.
    B c,end;acb_set_d(c.raw(),center);acb_set_d(end.raw(),next);B step=end-c;
    Jet tau(0,order+4,256);tau.set(0,c);tau.set(1,B(1));
    std::map<std::string,Jet> env;
    for(unsigned i=0;i<names.size();++i) {
      auto a=tau.constant(0),b=tau.constant(0),im=tau.constant(0);
      a.set(0,B::from_strings(initial.at(names[i]).rational().str()));
      b.set(0,B::from_strings(final.at(names[i]).rational().str()));im.set(0,B::from_strings("0","1"));
      auto direction=tau.constant(0);direction.set(0,B::from_strings((Rational(10+i)/Rational(25)).str()));
      env.emplace(names[i],a+(b-a)*tau+im*direction*tau*(tau.constant(1)-tau));
    }
    std::vector<Jet> root_polynomials,root_jets;
    for(unsigned i=0;i<3;++i) {
      auto polynomial=evaluate(radicands[i],tau,env);
      auto root=polynomial.sqrt();
      if(!continued_roots.empty()) {
        auto principal=root.at(0),opposite=-principal;
        bool plus=acb_overlaps(principal.raw(),continued_roots[i].raw());
        bool minus=acb_overlaps(opposite.raw(),continued_roots[i].raw());
        if(plus==minus)throw std::runtime_error("ambiguous square-root sheet at chart transition");
        if(minus)root=-root;
      }
      env.emplace(root_names[i],root);root_polynomials.push_back(std::move(polynomial));root_jets.push_back(std::move(root));
    }
    std::vector<Jet> letter_jets;
    for(auto l:active)letter_jets.push_back(evaluate(alphabet.args[l],tau,env));
    current=canonical_chart(entries,letter_jets,current,B(0),step,order);
    continued_roots.clear();
    for(unsigned i=0;i<3;++i)continued_roots.push_back(continue_polynomial_sqrt(root_polynomials[i],B(0),step,root_jets[i].at(0)));
    center=next;
    if(charts%20==0 || center==1)std::cout<<"Planar "<<family<<": chart "<<charts<<", parameter "<<center<<std::endl;
  }
  double max_error=0;
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<5;++k) {
    max_error=std::max(max_error,finite_reference_error(current[i][k],reference[i][k],256));
  }
  std::cout<<"Planar "<<family<<": maximum reference discrepancy = "<<max_error<<", charts = "<<charts
    <<", seconds = "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<'\n';
  std::cout<<"Root clearance and polynomial root-sheet continuation certified; Taylor remainder is not yet a publication certificate.\n";
  return std::isfinite(max_error) && max_error<1e-8?0:1;
}
} // namespace diffexp

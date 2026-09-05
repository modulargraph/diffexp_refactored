#pragma once
#include "diffexp/rational_transport.hpp"
#include <functional>

namespace diffexp {
// These are the two original contour orderings, not independent boundary seeds.
// Source: Docs/OriginalDiffExpExamples.md, unequal-mass qualification.
enum class OriginalBananaUnequalRoute { momentum_first, mass_first };
inline OriginalBananaUnequalRoute original_banana_unequal_route(const std::string& name) {
  if(name=="momentum-first")return OriginalBananaUnequalRoute::momentum_first;
  if(name=="mass-first")return OriginalBananaUnequalRoute::mass_first;
  throw std::invalid_argument("unequal banana route must be momentum-first or mass-first");
}
struct OriginalBananaUnequalResult {
  Boundary values; // All fifteen original masters, epsilon0..7.
  OriginalBananaUnequalRoute route;
  unsigned taylor_order=64;
  slong working_bits=384;
  bool full_tails_certified=false;
  std::vector<std::string> phases;
};
namespace original_banana_unequal_detail {
inline Exact reduce_i(const Exact& c) {
  return reduce_square(c.numerator(),1,c.constant(-1))/reduce_square(c.denominator(),1,c.constant(-1));
}
inline Boundary equal_leg(const std::string& directory,Boundary boundary,const std::string& from,const std::string& to) {
  ExactField field({"x","I"});Exact x(field,"x"),a(field,from),b(field,to);auto line=a+(b-a)*x;
  std::vector<RationalLineEntry> entries;
  for(unsigned e=0;e<2;++e) {
    auto matrix=data::read_file(directory+"/Data/Banana/EqualMass/dt_"+std::to_string(e)+".m");
    if(matrix.args.size()!=4)throw std::invalid_argument("equal banana matrix dimension");
    for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j) {
      auto c=reduce_i(evaluate_exact(matrix.args.at(i).args.at(j),x,{{"t",line}})*(b-a));
      if(!c.is_zero())entries.push_back({i,j,e,std::move(c)});
    }
  }
  return rational_line(entries,std::move(boundary),64);
}
inline Boundary lift(const Boundary& equal) {
  if(equal.size()!=4)throw std::invalid_argument("equal banana lift requires four masters");
  Boundary out;for(unsigned i:{0,0,0,0,0,0,1,1,1,1,2,3,3,3,3})out.push_back(equal.at(i));return out;
}
inline Boundary mass_leg(const std::string& directory,Boundary boundary,const Rational& psq) {
  ExactField field({"x","I"});Exact x(field,"x"),one(field,1);
  std::map<std::string,Exact> point{{"psq",x.constant(psq)},{"mm1",one+x},{"mm2",one+x/x.constant(2)},{"mm3",one+x/x.constant(3)},{"mm4",one}};
  std::map<std::tuple<unsigned,unsigned,unsigned>,Exact> combined;
  for(unsigned mass=1;mass<=3;++mass)for(unsigned e=0;e<2;++e) {
    auto matrix=data::read_file(directory+"/Data/Banana/UnequalMass/dmm"+std::to_string(mass)+"_"+std::to_string(e)+".m");
    if(matrix.args.size()!=15)throw std::invalid_argument("unequal banana matrix dimension");
    for(unsigned i=0;i<15;++i)for(unsigned j=0;j<15;++j) {
      auto c=evaluate_exact(matrix.args.at(i).args.at(j),x,point)/x.constant(mass);
      if(c.is_zero())continue;
      auto [p,inserted]=combined.try_emplace({i,j,e},x.constant(0));p->second=p->second+c;
    }
  }
  std::vector<RationalLineEntry> entries;
  for(auto& [key,c]:combined)if(!c.is_zero()){auto [i,j,e]=key;entries.push_back({i,j,e,std::move(c)});}
  return rational_line(entries,std::move(boundary),64);
}
inline Boundary momentum_leg(const std::string& directory,Boundary boundary,const std::string& from,const std::string& to) {
  ExactField field({"x","I"});Exact x(field,"x"),a(field,from),b(field,to);auto line=a+(b-a)*x;
  std::map<std::string,Exact> point{{"psq",x},{"mm1",x.constant(2)},{"mm2",x.constant(Rational("3/2"))},{"mm3",x.constant(Rational("4/3"))},{"mm4",x.constant(1)}};
  std::vector<RationalLineEntry> entries;
  for(unsigned e=0;e<2;++e) {
    auto matrix=data::read_file(directory+"/Data/Banana/UnequalMass/dpsq_"+std::to_string(e)+".m");
    if(matrix.args.size()!=15)throw std::invalid_argument("unequal banana momentum matrix dimension");
    for(unsigned i=0;i<15;++i)for(unsigned j=0;j<15;++j) {
      // Specialize fixed masses before the complex affine substitution. This
      // cancels large source polynomials in Q(psq), avoiding repeated I
      // expansions inside the original matrix AST. The substitution is exact.
      auto real_connection=evaluate_exact(matrix.args.at(i).args.at(j),x,point);
      std::array<Exact,2> replacements{line,Exact(field,"I")};
      auto c=reduce_i(real_connection.substitute(replacements)*(b-a));
      if(!c.is_zero())entries.push_back({i,j,e,std::move(c)});
    }
  }
  return rational_line(entries,std::move(boundary),64);
}
}
inline OriginalBananaUnequalResult original_banana_unequal(const std::string& directory,
    OriginalBananaUnequalRoute route=OriginalBananaUnequalRoute::momentum_first,
    const std::function<void(const std::string&)>& progress={}) {
  using B=Jet::Ball;const auto previous=B::precision();struct Restore{slong bits;~Restore(){B::set_precision(bits);}}restore{previous};B::set_precision(384);
  if(route!=OriginalBananaUnequalRoute::momentum_first && route!=OriginalBananaUnequalRoute::mass_first)throw std::invalid_argument("unknown unequal banana route");
  OriginalBananaUnequalResult result;result.route=route;
  auto phase=[&](const std::string& s){result.phases.push_back(s);if(progress)progress(s);};
  auto value=read_boundary(data::read_file(directory+"/Reference/BananaBoundaryAtMinusOneEps7.m"),4,7,384);
  using namespace original_banana_unequal_detail;
  if(route==OriginalBananaUnequalRoute::momentum_first) {
    for(const auto& [a,b]:std::vector<std::pair<std::string,std::string>>{{"-1","-1+10*I"},{"-1+10*I","50+10*I"},{"50+10*I","50"}}){phase("equal momentum "+a+" -> "+b);value=equal_leg(directory,std::move(value),a,b);}
    value=lift(value);phase("mass deformation at p_squared=50");value=mass_leg(directory,std::move(value),Rational(50));
  } else {
    // Explicit connector from the shared external seed; t=0 is avoided on
    // the same upper causal side as the original momentum-first route.
    for(const auto& [a,b]:std::vector<std::pair<std::string,std::string>>{{"-1","-1+I"},{"-1+I","1/2+I"},{"1/2+I","1/2"}}){phase("equal seed connector "+a+" -> "+b);value=equal_leg(directory,std::move(value),a,b);}
    value=lift(value);phase("mass deformation at p_squared=1/2");value=mass_leg(directory,std::move(value),Rational("1/2"));
    for(const auto& [a,b]:std::vector<std::pair<std::string,std::string>>{{"1/2","25+20*I"},{"25+20*I","50"}}){phase("unequal momentum "+a+" -> "+b);value=momentum_leg(directory,std::move(value),a,b);}
  }
  if(value.size()!=15)throw std::runtime_error("unequal endpoint shape");
  for(const auto& row:value){if(row.size()!=8)throw std::runtime_error("unequal endpoint epsilon window");for(const auto& c:row)if(!c.is_finite())throw std::runtime_error("nonfinite unequal endpoint");}
  result.values=std::move(value);return result;
}
inline int run_original_banana_unequal_route(const std::string& directory,
    OriginalBananaUnequalRoute route=OriginalBananaUnequalRoute::momentum_first) {
  using B=Jet::Ball;const auto previous=B::precision();struct Restore{slong p;~Restore(){B::set_precision(p);}}restore{previous};B::set_precision(384);
  std::cout<<"Original unequal banana route="<<(route==OriginalBananaUnequalRoute::momentum_first?"momentum-first":"mass-first")
    <<", Taylor order64, working bits384, evolved epsilon0..7; published comparison0..4.\n";
  auto result=original_banana_unequal(directory,route,[](const std::string& phase){std::cout<<phase<<std::endl;});
  auto reference=read_boundary(data::read_file(directory+"/Reference/BananaUnequalMassAt50.m"),15,4,384);
  double maximum=0;
  for(unsigned i=0;i<15;++i)for(unsigned k=0;k<5;++k) {
    const auto difference=result.values.at(i).at(k)-reference.at(i).at(k);
    if(!difference.is_finite()){std::cout<<"Nonfinite published-reference difference.\n";return 1;}
    arf_t bound;arf_init(bound);acb_get_abs_ubound_arf(bound,difference.raw(),384);
    maximum=std::max(maximum,arf_get_d(bound,ARF_RND_CEIL));arf_clear(bound);
  }
  std::cout<<"Maximum published-reference discrepancy="<<maximum<<"; threshold1e-10.\n"
    <<"Numerical comparison only; final Taylor tails are not certified.\n";
  return std::isfinite(maximum)&&maximum<1e-10?0:1;
}
} // namespace diffexp

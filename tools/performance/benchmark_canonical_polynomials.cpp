// Experimental canonical chart only. Does not select the production solver.
#include "diffexp/transport.hpp"
#include <fstream>
#include <iostream>
using namespace diffexp;
using B=kernel::ComplexBall;
struct Poly {
  acb_poly_t p;Poly(){acb_poly_init(p);}~Poly(){acb_poly_clear(p);}
  Poly(const Poly&)=delete;Poly& operator=(const Poly&)=delete;
};
struct Candidate {Boundary values;std::vector<B> coefficients;};
Candidate polynomial_chart(const transport::Compiled& c,const Boundary& initial,
    const std::vector<B>& roots,const B& center,const B& step,unsigned order) {
  if(!c.canonical)throw std::invalid_argument("canonical input required");
  const auto bits=B::precision();const auto d=c.dimension,w=c.epsilon_order+1;
  auto letters=std::make_unique<Poly[]>(c.letters.size());auto y=std::make_unique<Poly[]>(d*w);
  Jet x(0,order+2,bits);x.set(0,center);x.set(1,B(1));std::map<std::string,Jet> vars{{"x",x}};
  for(unsigned r=0;r<c.squares.size();++r) {
    auto square=evaluate(data::Reader(c.squares[r].str()).read(),x,vars),root=square.sqrt();
    auto factor=root.constant(0);factor.set(0,roots[r]/root.at(0));vars.emplace("r"+std::to_string(r),root*factor);
  }
  for(unsigned l=0;l<c.letters.size();++l) {
    auto value=evaluate(data::Reader(c.letters[l].str()).read(),x,vars);auto jet=value.derivative()/value;
    for(unsigned n=0;n<order;++n){auto a=jet.at(n);acb_poly_set_coeff_acb(letters[l].p,n,a.raw());}
  }
  for(unsigned i=0;i<d;++i)for(unsigned e=0;e<w;++e)acb_poly_set_coeff_acb(y[i*w+e].p,0,initial[i][e].raw());
  std::map<std::pair<unsigned,unsigned>,unsigned> ids;std::vector<unsigned> entry_ids;
  for(const auto& entry:c.canonical_entries){auto [it,added]=ids.try_emplace({entry.letter,entry.column},ids.size());entry_ids.push_back(it->second);}
  auto products=std::make_unique<Poly[]>(ids.size());Poly weighted,rhs;
  for(unsigned e=1;e<w;++e) {
    for(const auto& [key,id]:ids)acb_poly_mullow(products[id].p,letters[key.first].p,y[key.second*w+e-1].p,order,bits);
    // Each epsilon layer depends only on the completed previous layer.
    for(unsigned row=0;row<d;++row) {
      acb_poly_zero(rhs.p);
      for(unsigned p=0;p<c.canonical_entries.size();++p)if(c.canonical_entries[p].row==row) {
        auto weight=B::from_strings(c.canonical_entries[p].coefficient.str());
        acb_poly_scalar_mul(weighted.p,products[entry_ids[p]].p,weight.raw(),bits);
        acb_poly_add(rhs.p,rhs.p,weighted.p,bits);
      }
      acb_poly_integral(y[row*w+e].p,rhs.p,bits);
      acb_poly_set_coeff_acb(y[row*w+e].p,0,initial[row][e].raw());
    }
  }
  Candidate result{Boundary(d,std::vector<B>(w)),std::vector<B>((order+1)*d*w)};
  for(unsigned i=0;i<d;++i)for(unsigned e=0;e<w;++e) {
    acb_poly_evaluate(result.values[i][e].raw(),y[i*w+e].p,step.raw(),bits);
    for(unsigned n=0;n<=order;++n)acb_poly_get_coeff_acb(result.coefficients[(n*d+i)*w+e].raw(),y[i*w+e].p,n);
  }
  return result;
}
int main(int argc,char** argv) {
 try {
  if(argc!=3)throw std::invalid_argument("usage: benchmark_canonical_polynomials REQUEST_JSON EXACT_STEP");
  std::ifstream input(argv[1]);std::string text((std::istreambuf_iterator<char>(input)),{});
  const auto request=boost::json::parse(text).as_object();
  const auto d=transport::integer(request,"dimension",-1,1,1000),k=transport::integer(request,"epsilon_order",0,0,20);
  const auto order=transport::integer(request,"taylor_order",50,2,1000),bits=transport::integer(request,"working_bits",384,64,10000);
  B::set_precision(bits);auto compiled=transport::compile(request,d,k,bits);auto initial=transport::boundary(request,d,k,bits);
  if(!compiled.canonical)throw std::invalid_argument("benchmark requires a canonical request");
  auto entries=transport::numerical_entries(compiled);std::vector<B> roots;
  for(const auto& polynomial:compiled.square_polynomials){B value;acb_sqrt(value.raw(),polynomial.at(0).raw(),bits);roots.push_back(value);}
  const auto step=B::from_strings(Rational(argv[2]).str());
  auto start=std::chrono::steady_clock::now();auto baseline=transport::chart(compiled,entries,initial,roots,B(0),step,order,false,30,true,10000000,true,false);
  auto mid=std::chrono::steady_clock::now();auto candidate=polynomial_chart(compiled,initial,roots,B(0),step,order);
  auto end=std::chrono::steady_clock::now();
  unsigned overlaps=0,finite=0;long accuracy=bits;double max_radius_ratio=0;
  for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e) {
    const auto& a=baseline.values[i][e];const auto& b=candidate.values[i][e];
    finite+=b.is_finite();overlaps+=acb_overlaps(a.raw(),b.raw());accuracy=std::min(accuracy,acb_rel_accuracy_bits(b.raw()));
    const auto ar=transport::arithmetic_error(a),br=transport::arithmetic_error(b);
    const auto denominator=arf_get_d(arb_midref(acb_realref(ar.raw())),ARF_RND_NEAR);
    const auto numerator=arf_get_d(arb_midref(acb_realref(br.raw())),ARF_RND_NEAR);
    if(denominator>0)max_radius_ratio=std::max(max_radius_ratio,numerator/denominator);
  }
  if(finite!=d*(k+1) || overlaps!=finite)throw std::runtime_error("candidate finite/overlap check failed");
  const double old_seconds=std::chrono::duration<double>(mid-start).count(),new_seconds=std::chrono::duration<double>(end-mid).count();
  boost::json::object report{{"status","passed"},{"dimension",d},{"epsilon_order",k},{"taylor_order",order},{"working_bits",bits},{"step",argv[2]},
    {"production_chart_seconds",old_seconds},{"polynomial_chart_seconds",new_seconds},{"observed_chart_speedup",old_seconds/new_seconds},
    {"finite",finite},{"overlapping",overlaps},{"minimum_relative_accuracy_bits_including_zero_balls",accuracy},{"maximum_radius_ratio",max_radius_ratio},
    {"scope","Experimental retained first chart only; overlap is not an independent accuracy or tail certificate. Candidate does not select the production solver. Preparation excluded; polynomial coefficient extraction included in candidate time."}};
  std::cout<<boost::json::serialize(report)<<'\n';
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}

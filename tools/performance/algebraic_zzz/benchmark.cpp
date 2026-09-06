#include "algebra.hpp"
#include <queue>
#include "products.hpp"
int main(int argc,char**argv){try{
 if(argc<4)throw std::runtime_error("usage: benchmark REQUEST DECOMPOSITION ORDER [STEP]");
 auto load=[](const char*p){std::ifstream f(p);return boost::json::parse(std::string(std::istreambuf_iterator<char>(f),{})).as_object();};
 auto req=load(argv[1]),dec=load(argv[2]);unsigned N=std::stoul(argv[3]);if(N<8||N>1000)throw std::runtime_error("Taylor order must be 8..1000");unsigned d=req.at("dimension").as_int64(),k=req.at("epsilon_order").as_int64(),bits=req.at("working_bits").as_int64();B::set_precision(bits);
 auto t=Clock::now();auto c=transport::compile(req,d,k,bits);auto initial=transport::boundary(req,d,k,bits);double original_prep=seconds(t);t=Clock::now();
 std::vector<std::map<unsigned,Exact>> parts;for(auto&l:dec.at("decomposed_letters").as_array()){std::map<unsigned,Exact> p;for(auto&v:l.as_array()){auto&o=v.as_object();p.emplace(unsigned(o.at("mask").as_int64()),Exact(c.field,transport::string(o.at("coefficient"))));}parts.push_back(std::move(p));}
 std::vector<std::vector<std::pair<unsigned,unsigned>>> graph(d);unsigned conflicts=0;
 for(const auto&e:c.canonical_entries)for(auto&[m,a]:parts[e.letter])if(!a.is_zero()){graph[e.row].push_back({e.column,m});graph[e.column].push_back({e.row,m});}
 std::vector<int> gauge(d,-1);unsigned components=0;
 for(unsigned i=0;i<d;++i)if(gauge[i]<0){++components;gauge[i]=0;std::queue<unsigned> todo;todo.push(i);while(!todo.empty()){auto row=todo.front();todo.pop();for(auto[col,m]:graph[row]){int g=gauge[row]^m;if(gauge[col]<0){gauge[col]=g;todo.push(col);}else if(gauge[col]!=g)++conflicts;}}}
 if(conflicts){std::cout<<boost::json::serialize(boost::json::object{{"gauge_conflicts",conflicts}})<<'\n';return 2;}
 std::cerr<<"Root gauge consistent, components="<<components<<"\n";
 std::vector<Exact> products;for(unsigned m=0;m<(1u<<c.squares.size());++m){Exact p(c.field,1);for(unsigned r=0;r<c.squares.size();++r)if(m&(1u<<r))p=p*c.squares[r];products.push_back(p);}
 std::vector<WeightedTerm> rational;
 for(auto&e:c.canonical_entries)for(auto&[m,a]:parts[e.letter]){
  if(unsigned(gauge[e.row]^gauge[e.column])!=m)throw std::runtime_error("gauge parity verification");
  auto coefficient=c.reduced(a*products[unsigned(gauge[e.row])&m]);
  if(!coefficient.is_zero())rational.push_back({e.row,e.column,1,coefficient,e.coefficient});
 }
 for(unsigned i=0;i<d;++i)if(gauge[i])rational.push_back({i,i,0,products[gauge[i]].derivative(0)/(products[gauge[i]]*Exact(c.field,2)),Rational(1)});
 auto compiled=prepare_products(rational);
 unsigned maxdegree=0;for(auto&p:compiled)maxdegree=std::max(maxdegree,unsigned(std::max(p.p.size(),p.q.size())-1));
 double preparation=seconds(t);std::cerr<<"Products compiled "<<preparation<<"s maxdegree="<<maxdegree<<" groups="<<compiled.size()<<"\n";
 std::vector<B> roots;for(auto&p:c.square_polynomials){B r;acb_sqrt(r.raw(),p.at(0).raw(),bits);roots.push_back(r);}
 double next=clearance_endpoint(0,c.singularities);B end;acb_set_d(end.raw(),std::min(1.,next));if(argc>4)end=B::from_strings(Rational(argv[4]).str());
 std::vector<B> start_factors(d,B(1)),end_factors(d,B(1));auto boundary=initial;
 for(unsigned i=0;i<d;++i)for(unsigned r=0;r<roots.size();++r)if(gauge[i]&(1u<<r)){start_factors[i]=start_factors[i]*roots[r];end_factors[i]=end_factors[i]*continue_polynomial_sqrt(c.square_polynomials[r],B(0),end,roots[r]);}
 for(unsigned i=0;i<d;++i)for(auto&b:boundary[i])b=b*start_factors[i];
 auto entries=transport::numerical_entries(c);transport::Chart baseline;Boundary candidate;double base_time=0,candidate_time=0;
 const bool candidate_first=std::getenv("CANDIDATE_FIRST")!=nullptr;
 auto base_run=[&]{t=Clock::now();baseline=transport::chart(c,entries,initial,roots,B(0),end,N,false,40,true,10000000,true,false);base_time=seconds(t);};
 auto candidate_run=[&]{t=Clock::now();candidate=product_chart(compiled,boundary,end,N);candidate_time=seconds(t);};
 if(candidate_first){candidate_run();base_run();}else{base_run();candidate_run();}
 std::cerr<<"baseline "<<base_time<<"s candidate "<<candidate_time<<"s\n";
 unsigned overlap=0;double max_difference=0,max_radius_ratio=0;for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e){candidate[i][e]=candidate[i][e]/end_factors[i];if(!candidate[i][e].is_finite())throw std::runtime_error("nonfinite candidate");overlap+=acb_overlaps(candidate[i][e].raw(),baseline.values[i][e].raw());auto delta=transport::magnitude(candidate[i][e]-baseline.values[i][e]);max_difference=std::max(max_difference,arf_get_d(arb_midref(acb_realref(delta.raw())),ARF_RND_NEAR));auto ar=transport::arithmetic_error(baseline.values[i][e]),br=transport::arithmetic_error(candidate[i][e]);double aa=arf_get_d(arb_midref(acb_realref(ar.raw())),ARF_RND_NEAR),bb=arf_get_d(arb_midref(acb_realref(br.raw())),ARF_RND_NEAR);if(aa>0)max_radius_ratio=std::max(max_radius_ratio,bb/aa);}
 boost::json::array masks;for(int g:gauge)masks.push_back(g);
 std::cout<<boost::json::serialize(boost::json::object{{"status","completed"},{"candidate_first",candidate_first},{"dimension",d},{"epsilon_high",k},{"working_bits",bits},{"taylor_order",N},{"original_prepare_seconds",original_prep},{"rational_prepare_seconds",preparation},{"max_polynomial_degree",maxdegree},{"fallback_rows",0},{"product_groups",compiled.size()},{"rational_entries",rational.size()},{"baseline_seconds",base_time},{"candidate_seconds",candidate_time},{"speedup",base_time/candidate_time},{"overlaps",overlap},{"values",d*(k+1)},{"max_difference_bound_midpoint",max_difference},{"max_radius_ratio",max_radius_ratio},{"step",transport::decimal_mid(acb_realref(end.raw()),30)},{"gauge_masks",masks},{"baseline_values",transport::matrix_json(baseline.values,290)},{"candidate_values",transport::matrix_json(candidate,290)}})<<'\n';
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}

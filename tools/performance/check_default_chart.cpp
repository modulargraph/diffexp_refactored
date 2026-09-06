#include "diffexp/transport.hpp"
int main(int argc,char**argv){try {
 using namespace diffexp;using namespace diffexp::transport;
 if(argc!=2)throw std::invalid_argument("usage: check_default_chart REQUEST");
 std::ifstream file(argv[1]);auto r=json::parse(std::string(std::istreambuf_iterator<char>(file),{})).as_object();
 unsigned d=integer(r,"dimension",0,1,1000),k=integer(r,"epsilon_order",4,0,100),n=integer(r,"taylor_order",50,8,1000);
 B::set_precision(integer(r,"working_bits",384,64,100000));auto c=compile(r,d,k,B::precision());
 auto initial=boundary(r,d,k,B::precision());auto roots=principal_roots_at(c,B(0));auto entries=numerical_entries(c);
 auto p=finite_lag_plan(c);if(!p)throw std::runtime_error("no finite-lag plan");
 B step;acb_set_d(step.raw(),std::min(1.,clearance_endpoint(0,c.singularities)*4./integer(r,"division_order",4,2,100)));
 auto start=std::chrono::steady_clock::now();auto old=chart(c,entries,initial,roots,B(0),step,n,false,40,true,10000000,true,false);
 double old_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 start=std::chrono::steady_clock::now();auto current=chart(c,entries,initial,roots,B(0),step,n,false,40);
 double new_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 if(!current.finite_lag)throw std::runtime_error("default first chart fell back");
 B worst;auto goal=integer(r,"accuracy_goal",0,0,20000);B tolerance=B::from_strings("1e-"+std::to_string(goal+8));
 for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e) {
  B a,b;acb_get_mid(a.raw(),current.values[i][e].raw());acb_get_mid(b.raw(),old.values[i][e].raw());
  auto difference=magnitude(a-b)/(magnitude(b)+B(1));arb_max(acb_realref(worst.raw()),acb_realref(worst.raw()),acb_realref(difference.raw()),B::precision());
  if(goal && !arb_le(acb_realref(current.truncation_errors[i][e].raw()),acb_realref(((magnitude(current.values[i][e])+B(1))*tolerance).raw())))throw std::runtime_error("default first-chart tail check failed");
 }
 std::cout<<json::serialize(json::object{{"finite_lag",true},{"order",n},{"dimension",d},{"series_seconds",old_seconds},{"default_seconds",new_seconds},{"maximum_normalized_midpoint_difference",decimal_mid(acb_realref(worst.raw()),20)}})<<'\n';
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}

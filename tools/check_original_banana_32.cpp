// Historical Banana.nb march to32. The independent saved table is at20;
// the endpoint32 check compares two homotopic paths from the same seed.
#include "diffexp/singular_transport.hpp"
#include <boost/json.hpp>
#include <iostream>
using namespace diffexp;using B=Jet::Ball;namespace json=boost::json;
int main(int argc,char** argv) {
  std::ostream output(std::cout.rdbuf());
  struct Redirect {std::streambuf* old=std::cout.rdbuf(std::cerr.rdbuf());~Redirect(){std::cout.rdbuf(old);}} redirect;
  try {
    if(argc!=2)throw std::invalid_argument("requires original-example directory");
    const auto begin=std::chrono::steady_clock::now();const std::string directory=argv[1];
    auto real=original_banana_equal_real(directory,Rational(32));B::set_precision(384);
    auto contour=read_boundary(original_banana_reference(directory,"boundaryAtMinusOne"),4,4,384);
    const auto a0=data::read_file(directory+"/Data/Banana/EqualMass/dt_0.m"),a1=data::read_file(directory+"/Data/Banana/EqualMass/dt_1.m");
    for(const auto& [from,to]:std::vector<std::pair<std::string,std::string>>{{"-1","-1+5*I"},{"-1+5*I","32+5*I"},{"32+5*I","32"}}) {
      ExactField field({"x","I"});Exact x(field,"x"),a(field,from),b(field,to);auto line=a+(b-a)*x;
      std::vector<RationalLineEntry> entries;
      for(unsigned order=0;order<2;++order) {
        const auto& matrix=order?a1:a0;
        if(matrix.args.size()!=4)throw std::invalid_argument("equal banana matrix dimension");
        for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j) {
          auto coefficient=evaluate_exact(matrix.args.at(i).args.at(j),x,{{"t",line}})*(b-a);
          coefficient=reduce_square(coefficient.numerator(),1,x.constant(-1))/reduce_square(coefficient.denominator(),1,x.constant(-1));
          if(!coefficient.is_zero())entries.push_back({i,j,order,std::move(coefficient)});
        }
      }
      contour=rational_line(entries,std::move(contour),50);
    }
    auto text=[](const B& value) {
      auto part=[](const arb_t b){char* p=arb_get_str(b,45,0);std::string result(p);flint_free(p);return result;};
      return json::object{{"real",part(acb_realref(value.raw()))},{"imaginary",part(acb_imagref(value.raw()))}};
    };
    json::array coefficients;double maximum=0;bool pass=real.maximum_discrepancy<1e-10;
    for(unsigned i=0;i<4;++i)for(unsigned k=0;k<5;++k) {
      const auto error=singular_transport_detail::finite_discrepancy(real.value[i][k],contour[i][k]);maximum=std::max(maximum,error);
      const bool good=std::isfinite(error)&&error<1e-20;pass&=good;
      coefficients.push_back(json::object{{"master",i+1},{"epsilon",k},{"real_path",text(real.value[i][k])},
        {"contour_path",text(contour[i][k])},{"difference",text(real.value[i][k]-contour[i][k])},{"pass",good}});
    }
    const json::object report{{"schema","DiffExp.OriginalEqualBanana32/v1"},{"status",pass?"pass":"fail"},
      {"start","-1"},{"endpoint","32"},{"reference_endpoint","20"},{"masters",4},{"epsilon_low",0},{"epsilon_high",4},
      {"taylor_order",50},{"working_bits",384},{"singular_matches",real.singular_matches},
      {"maximum_reference_discrepancy_at20",real.maximum_discrepancy},{"reference_threshold_at20",1e-10},
      {"maximum_route_discrepancy_at32",maximum},{"route_threshold_at32",1e-20},{"coefficients",coefficients},
      {"real_route","-1 through local +i0 matches at0,4,16; saved comparison at20; continue to32"},
      {"contour_route","-1 -> -1+5i ->32+5i ->32"},{"shared_external_seed","Reference/BananaBoundaryAtMinusOneEps7.m"},
      {"independent_reference_at32",false},{"omitted_tails_certified",false},
      {"seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()}};
    output<<json::serialize(report)<<'\n';return pass?0:1;
  }catch(const std::exception& e){output<<json::serialize(json::object{{"schema","DiffExp.OriginalEqualBanana32/v1"},{"status","error"},{"error",e.what()}})<<'\n';return 1;}
}

#include "algebra.hpp"
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("usage: probe request.json");
 std::ifstream in(argv[1]);auto req=boost::json::parse(std::string(std::istreambuf_iterator<char>(in),{})).as_object();
 unsigned d=req.at("dimension").as_int64(),k=req.at("epsilon_order").as_int64(),bits=req.at("working_bits").as_int64();B::set_precision(bits);
 auto start=Clock::now();auto c=transport::compile(req,d,k,bits);std::cerr<<"compiled "<<seconds(start)<<"s roots="<<c.squares.size()<<" letters="<<c.letters.size()<<"\n";
 std::vector<std::map<unsigned,Exact>> parts;unsigned termcount=0,maxdegree=0;
 for(unsigned l=0;l<c.letters.size();++l){auto q=c.derivative(c.letters[l])/c.letters[l];auto p=decompose(c,c.reduced(q));termcount+=p.size();for(auto&[m,a]:p)for(const auto&ts:{a.numerator_terms(),a.denominator_terms()})for(auto&t:ts)maxdegree=std::max(maxdegree,unsigned(t.powers[0]));parts.push_back(std::move(p));std::cerr<<"letter "<<l<<" parts="<<parts.back().size()<<" time="<<seconds(start)<<"\n";}
 boost::json::array letters;for(auto&p:parts){boost::json::array ts;for(auto&[m,a]:p)ts.push_back(boost::json::object{{"mask",m},{"coefficient",a.str()}});letters.push_back(ts);}
 boost::json::array squares;for(auto&q:c.squares)squares.push_back(boost::json::value(q.str()));
 boost::json::object report{{"roots",c.squares.size()},{"letters",c.letters.size()},{"components",d},{"letter_terms",termcount},{"max_rational_degree",maxdegree},{"preparation_seconds",seconds(start)},{"squares",squares},{"decomposed_letters",letters}};
 std::cout<<boost::json::serialize(report)<<'\n';return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}

#pragma once
#include "ibp/field.hpp"
#include <boost/json.hpp>
#include <flint/fmpq_mat.h>
#include <fstream>
#include <bit>
#include <tuple>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <numeric>
namespace ibp {
namespace json=boost::json;
struct ExactMatrix {
 fmpq_mat_t value;
 ExactMatrix(unsigned r,unsigned c){fmpq_mat_init(value,r,c);}
 ~ExactMatrix(){fmpq_mat_clear(value);}
 ExactMatrix(const ExactMatrix&)=delete;
 ExactMatrix& operator=(const ExactMatrix&)=delete;
};
inline std::string rational_text(const json::value& value) {
 if(value.is_string())return std::string(value.as_string());
 if(value.is_int64())return std::to_string(value.as_int64());
 throw std::invalid_argument("family coefficients must be rational strings or integers");
}
struct Integral {
 std::array<std::int8_t,16> powers{};
 auto operator<=>(const Integral&)const=default;
};
struct IntegralHash {
 std::size_t operator()(const Integral& a)const {std::uint64_t h=1469598103934665603ULL;for(auto v:a.powers){h^=static_cast<unsigned char>(v);h*=1099511628211ULL;}return h;}
};
inline unsigned sector(const Integral& a,unsigned physical) {unsigned mask=0;for(unsigned i=0;i<physical;++i)if(a.powers[i]>0)mask|=1u<<i;return mask;}
inline auto grade(const Integral& a){unsigned sector_size=0,dots=0,numerators=0;for(auto v:a.powers)if(v>0){++sector_size;dots+=v-1;}else numerators-=v;return std::tuple{sector_size,dots,numerators,a.powers};}
using Affine=std::vector<Word>;
struct Geometry {
 unsigned loops,externals,physical,n;std::vector<Affine> denominators;
 std::vector<std::vector<Affine>> contractions;std::vector<bool> trace,zero_sectors;std::vector<std::pair<unsigned,unsigned>> auxiliary_scalar_products;
};
class Family {
 json::object input_;
 public:
 explicit Family(const std::string& file){std::ifstream in(file);if(!in)throw std::invalid_argument("cannot open family");input_=json::parse(std::string(std::istreambuf_iterator<char>(in),{})).as_object();}
 explicit Family(json::object input):input_(std::move(input)){}
 const json::object& json_input()const{return input_;}
 json::object basis_descriptor(const Geometry& g)const {
  json::array slots;const auto count=input_.at("propagators").as_array().size();
  for(unsigned k=0;k<count;++k)slots.push_back(json::object{{"kind","propagator"},{"input_index",k}});
  for(auto [i,j]:g.auxiliary_scalar_products)slots.push_back(json::object{{"kind","scalar_product"},{"vectors",json::array{i,j}}});
  return {{"schema","IBPSolver.Basis/v1"},{"physical_count",g.physical},{"slots",slots},{"convention","D = mass_squared - momentum_squared; vector indices are loops then external vectors; zero-based slots"}};
 }
 Geometry compile(const Field& f)const {
  const auto& lines=input_.at("propagators").as_array();const auto& gram=input_.at("external_gram").as_array();
  auto loops_value=input_.at("loops").as_int64(),physical_value=input_.at("physical_count").as_int64();
  if(loops_value<1||loops_value>4||physical_value<1||physical_value>12||gram.size()>16)throw std::invalid_argument("family size budget");
  unsigned l=loops_value,e=gram.size(),physical=physical_value;
  if(!l || l>4 || physical>12 || !physical || physical>lines.size())throw std::invalid_argument("family size budget");
  unsigned n=l*(l+1)/2+l*e;if(n>16 || physical>n || lines.size()>n)throw std::invalid_argument("at most16 independent scalar products; partial fractions needed for overcomplete denominators");
  for(const auto& row:gram)if(row.as_array().size()!=e)throw std::invalid_argument("Gram shape");
  ExactMatrix gram_exact(e,e);
  for(unsigned i=0;i<e;++i)for(unsigned j=0;j<e;++j){auto text=rational_text(gram[i].as_array()[j]);auto q=fmpq_mat_entry(gram_exact.value,i,j);if(fmpq_set_str(q,text.c_str(),10)||fmpz_is_zero(fmpq_denref(q)))throw std::invalid_argument("invalid Gram rational");fmpq_canonicalise(q);}
  for(unsigned i=0;i<e;++i)for(unsigned j=0;j<i;++j)if(!fmpq_equal(fmpq_mat_entry(gram_exact.value,i,j),fmpq_mat_entry(gram_exact.value,j,i)))throw std::invalid_argument("Gram matrix must be symmetric over Q");
  Geometry out{l,e,physical,n};std::vector<std::pair<unsigned,unsigned>> pairs;
  for(unsigned i=0;i<l;++i)for(unsigned j=i;j<l+e;++j)pairs.emplace_back(i,j);
  auto dot=[&](unsigned a,unsigned b){Affine q(n+1);if(a>b)std::swap(a,b);if(a>=l)q[n]=f.rational(rational_text(gram[a-l].as_array()[b-l]));else q[std::find(pairs.begin(),pairs.end(),std::pair{a,b})-pairs.begin()]=1;return q;};
  auto axpy=[&](Affine& a,const Affine& b,Word c){for(unsigned i=0;i<=n;++i)a[i]=f.add(a[i],f.mul(c,b[i]));};
  for(const auto& item:lines){const auto& line=item.as_object();auto lc=line.at("loop_coefficients").as_array(),ec=line.at("external_coefficients").as_array();if(lc.size()!=l || ec.size()!=e)throw std::invalid_argument("momentum shape");
   std::vector<Word> c;for(const auto& v:lc)c.push_back(f.rational(rational_text(v)));for(const auto& v:ec)c.push_back(f.rational(rational_text(v)));
   Affine q(n+1);q[n]=f.rational(rational_text(line.at("mass_squared")));for(unsigned i=0;i<l+e;++i)for(unsigned j=0;j<l+e;++j)axpy(q,dot(i,j),f.sub(0,f.mul(c[i],c[j])));out.denominators.push_back(std::move(q));
  }
  std::map<unsigned,Affine> echelon;
  auto independent=[&](Affine a){a[n]=0;for(const auto& [p,row]:echelon)axpy(a,row,f.sub(0,a[p]));unsigned p=0;while(p<n&&!a[p])++p;if(p==n)return false;auto inverse=f.inv(a[p]);for(auto& v:a)v=f.mul(v,inverse);echelon.emplace(p,std::move(a));return true;};
  for(const auto& q:out.denominators)if(!independent(q))throw std::domain_error("dependent propagators or exceptional prime");
  // Choose the auxiliary completion over Q. A modular minor may vanish at
  // an exceptional prime; it must not silently change an integral's meaning.
  ExactMatrix exact(n,n),reduced(n,n),routing(1,l+e);
  for(unsigned k=0;k<lines.size();++k){const auto& line=lines[k].as_object();const auto& lc=line.at("loop_coefficients").as_array();const auto& ec=line.at("external_coefficients").as_array();
   for(unsigned j=0;j<l+e;++j){auto text=rational_text(j<l?lc[j]:ec[j-l]);auto q=fmpq_mat_entry(routing.value,0,j);if(fmpq_set_str(q,text.c_str(),10)||fmpz_is_zero(fmpq_denref(q)))throw std::invalid_argument("invalid routing rational");fmpq_canonicalise(q);}
   for(unsigned z=0;z<n;++z){auto [i,j]=pairs[z];auto q=fmpq_mat_entry(exact.value,k,z);fmpq_mul(q,fmpq_mat_entry(routing.value,0,i),fmpq_mat_entry(routing.value,0,j));if(i!=j)fmpq_mul_2exp(q,q,1);fmpq_neg(q,q);}
  }
  unsigned exact_rank=lines.size();if(fmpq_mat_rref(reduced.value,exact.value)!=exact_rank)throw std::invalid_argument("dependent propagators over Q");
  for(unsigned i=0;i<n&&exact_rank<n;++i){fmpq_one(fmpq_mat_entry(exact.value,exact_rank,i));auto rank=fmpq_mat_rref(reduced.value,exact.value);if(rank>exact_rank){Affine q(n+1);q[i]=1;if(!independent(q))throw std::domain_error("exceptional prime for rational auxiliary completion");out.denominators.push_back(std::move(q));out.auxiliary_scalar_products.push_back(pairs[i]);++exact_rank;}else fmpq_zero(fmpq_mat_entry(exact.value,exact_rank,i));}
  std::vector<std::vector<Word>> matrix(n,std::vector<Word>(2*n));for(unsigned i=0;i<n;++i){for(unsigned j=0;j<n;++j)matrix[i][j]=out.denominators[i][j];matrix[i][n+i]=1;}
  for(unsigned i=0;i<n;++i){unsigned pivot=i;while(pivot<n&&!matrix[pivot][i])++pivot;if(pivot==n)throw std::domain_error("exceptional propagator basis");std::swap(matrix[i],matrix[pivot]);auto inverse=f.inv(matrix[i][i]);for(auto& v:matrix[i])v=f.mul(v,inverse);for(unsigned r=0;r<n;++r)if(r!=i){auto c=matrix[r][i];if(c)for(unsigned j=0;j<2*n;++j)matrix[r][j]=f.sub(matrix[r][j],f.mul(c,matrix[i][j]));}}
  std::vector<Affine> sp(n,Affine(n+1));for(unsigned i=0;i<n;++i)for(unsigned j=0;j<n;++j){sp[i][j]=matrix[i][n+j];sp[i][n]=f.sub(sp[i][n],f.mul(sp[i][j],out.denominators[j][n]));}
  for(unsigned i=0;i<l;++i)for(unsigned v=0;v<l+e;++v){out.trace.push_back(i==v);out.contractions.emplace_back();for(const auto& q:out.denominators){Affine derivative(n+1),rewritten(n+1);for(unsigned z=0;z<n;++z){auto [a,b]=pairs[z];if(a==i)axpy(derivative,dot(v,b),q[z]);if(b==i)axpy(derivative,dot(v,a),q[z]);}rewritten[n]=derivative[n];for(unsigned z=0;z<n;++z)axpy(rewritten,sp[z],derivative[z]);out.contractions.back().push_back(std::move(rewritten));}}
  // Prove free-loop sectors over Q, never infer zero sectors from modular rank.
  out.zero_sectors.resize(1u<<physical);
  for(unsigned mask=0;mask<(1u<<physical);++mask){unsigned count=std::popcount(mask);if(!count){out.zero_sectors[mask]=true;continue;}fmpq_mat_t exact,reduced;fmpq_mat_init(exact,count,l);fmpq_mat_init(reduced,count,l);unsigned r=0;
   for(unsigned k=0;k<physical;++k)if(mask&(1u<<k)){const auto& c=lines[k].as_object().at("loop_coefficients").as_array();for(unsigned j=0;j<l;++j){auto text=rational_text(c[j]);if(fmpq_set_str(fmpq_mat_entry(exact,r,j),text.c_str(),10)){fmpq_mat_clear(exact);fmpq_mat_clear(reduced);throw std::invalid_argument("invalid loop coefficient");}fmpq_canonicalise(fmpq_mat_entry(exact,r,j));}++r;}
   out.zero_sectors[mask]=fmpq_mat_rref(reduced,exact)<l;fmpq_mat_clear(exact);fmpq_mat_clear(reduced);
  }
  return out;
 }
};
}

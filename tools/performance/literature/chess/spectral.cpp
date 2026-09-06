// Standalone experimental sequential-epsilon Chebyshev collocation.
// Method credit: CHESS method by Yuanche Liu and Yang Zhang; release 0.1.0, pinned commit
// 4a701fc1332f29f6237d14427336e60615b966e3. This is an independent C++
// implementation of ordinary Lobatto differentiation/collocation, not copied
// CHESS source. It is not a production backend or an adaptive error estimator.
#include "diffexp/transport.hpp"
#include <fstream>
#include <sstream>
using namespace diffexp;
using B=kernel::ComplexBall;
namespace j=boost::json;
struct Mat {
 acb_mat_t v;
 Mat(long r,long c){acb_mat_init(v,r,c);}
 ~Mat(){acb_mat_clear(v);}
 Mat(const Mat&)=delete;
};
int main(int argc,char**argv){try{
 if(argc!=3)throw std::runtime_error("usage: spectral request.json nodes");
 const auto start=std::chrono::steady_clock::now();
 const unsigned n=std::stoul(argv[2]),m=n+1;
 if(n<2||n>128)throw std::runtime_error("nodes outside 2..128");
 std::ifstream in(argv[1]);std::stringstream buffer;buffer<<in.rdbuf();
 auto request=j::parse(buffer.str()).as_object();
 unsigned d=transport::integer(request,"dimension",0,1,1000),k=transport::integer(request,"epsilon_order",0,0,20);
 // CHESS's release precision=60 and WorkingPrecisionA=80 are rounded upward
 // to whole binary bits. Node construction and the recurrence use 200 bits.
 constexpr long bits=200,bitsA=266;
 B::set_precision(bitsA);
 auto c=transport::compile(request,d,k,bitsA);
 if(!c.canonical)throw std::runtime_error("experiment requires epsilon-linear scalar forms or dlog letters");
 if(request.if_contains("asymptotic")||request.if_contains("basis_prefactors"))throw std::runtime_error("experiment supports ordinary fixed-sheet boundaries only");
 for(const auto&s:c.singularities)if(arb_contains_zero(acb_imagref(s.raw()))&&arb_is_nonnegative(acb_realref(s.raw()))&&arb_le(acb_realref(s.raw()),acb_realref(B(1).raw())))throw std::runtime_error("singular point on the real path");
 auto bc=transport::boundary(request,d,k,bits);
 std::vector<B> nodes(m),weights(m);
 B::set_precision(bits);
 for(unsigned i=0;i<m;++i){B arg=B::from_strings(std::to_string(i)+"/"+std::to_string(n)),cosine;acb_cos_pi(cosine.raw(),arg.raw(),bits);nodes[i]=(B(1)-cosine)/B(2);weights[i]=B((i%2)?-1:1);if(i==0||i==n)weights[i]=weights[i]/B(2);}
 // Replace the first differentiation equation with the exact boundary row.
 // Its inverse maps sampled derivatives + the boundary to nodal values.
 Mat D(m,m),Q(m,m),previous(m,d),rhs(m,d),next(m,d);
 acb_one(acb_mat_entry(D.v,0,0));
 for(unsigned i=1;i<m;++i){B diagonal(0);for(unsigned a=0;a<m;++a)if(a!=i){B entry=weights[a]/weights[i]/(nodes[i]-nodes[a]);acb_set(acb_mat_entry(D.v,i,a),entry.raw());diagonal=diagonal-entry;}acb_set(acb_mat_entry(D.v,i,i),diagonal.raw());}
 if(!acb_mat_inv(Q.v,D.v,bits))throw std::runtime_error("scalar collocation inverse failed");
 const auto prepared=std::chrono::steady_clock::now();
 // Continue every root from its principal start through ordered nodes.
 B::set_precision(bitsA);
 std::vector<data::Expr> forms;for(unsigned l=0;l<c.letters.size();++l)forms.push_back(data::Reader(c.letter_derivative(l).str()).read());
 std::vector<std::vector<B>> sampled(m,std::vector<B>(forms.size()));
 auto root_values=transport::principal_roots_at(c,B(0));
 for(unsigned i=1;i<m;++i){Jet x(0,1,bitsA);x.set(0,nodes[i]);std::map<std::string,Jet> vars{{"x",x}};
  for(unsigned r=0;r<c.squares.size();++r){root_values[r]=continue_polynomial_sqrt(c.square_polynomials[r],nodes[i-1],nodes[i],root_values[r]);auto root=x.constant(0);root.set(0,root_values[r]);vars.emplace("r"+std::to_string(r),std::move(root));}
  for(unsigned l=0;l<forms.size();++l)sampled[i][l]=evaluate(forms[l],x,vars).at(0);
 }
 std::vector<B> coefficients;for(const auto&e:c.canonical_entries)coefficients.push_back(B::from_strings(e.coefficient.str()));
 const auto sampled_time=std::chrono::steady_clock::now();
 B::set_precision(bits);
 Boundary final(d,std::vector<B>(k+1,B(0)));
 for(unsigned row=0;row<d;++row){final[row][0]=bc[row][0];for(unsigned node=0;node<m;++node)acb_set(acb_mat_entry(previous.v,node,row),bc[row][0].raw());}
 for(unsigned epsilon=1;epsilon<=k;++epsilon){
  acb_mat_zero(rhs.v);for(unsigned row=0;row<d;++row)acb_set(acb_mat_entry(rhs.v,0,row),bc[row][epsilon].raw());
  B product;
  for(unsigned node=1;node<m;++node)for(unsigned a=0;a<c.canonical_entries.size();++a){const auto&e=c.canonical_entries[a];acb_mul(product.raw(),sampled[node][e.letter].raw(),acb_mat_entry(previous.v,node,e.column),bits);acb_addmul(acb_mat_entry(rhs.v,node,e.row),coefficients[a].raw(),product.raw(),bits);}
  acb_mat_mul(next.v,Q.v,rhs.v,bits);
  for(unsigned row=0;row<d;++row)acb_set(final[row][epsilon].raw(),acb_mat_entry(next.v,n,row));
  acb_mat_swap(previous.v,next.v);
 }
 const auto end=std::chrono::steady_clock::now();
 auto seconds=[](auto a,auto b){return std::chrono::duration<double>(b-a).count();};
 j::object output{{"schema","DiffExp.SpectralExperiment/v1"},{"method","ordinary Chebyshev-Lobatto sequential epsilon collocation"},{"nodes",n},{"working_bits",bits},{"matrix_working_bits",bitsA},{"values",transport::matrix_json(final,65)},{"charts",1},{"recurrence",j::object{{"type","experimental_spectral"}}},{"omitted_tails_certified",false},{"radius_scope","finite collocation arithmetic only; no spectral truncation estimate"},{"timings",j::object{{"preparation_seconds",seconds(start,prepared)},{"sampling_seconds",seconds(prepared,sampled_time)},{"recurrence_seconds",seconds(sampled_time,end)},{"numerical_seconds",seconds(prepared,end)},{"total_seconds",seconds(start,end)}}}};
 std::cout<<j::serialize(output)<<'\n';
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}

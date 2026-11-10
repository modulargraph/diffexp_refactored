#pragma once
#include "diffexp/adjoint_transport.hpp"
#include <set>
namespace diffexp::ft_ultraspherical_detail {
using B=Jet::Ball;using V=std::vector<B>;using M=std::vector<V>;
inline V mulx(const V&v){V r(v.size());for(unsigned k=0;k<v.size();++k){B h;acb_mul_2exp_si(h.raw(),v[k].raw(),-1);if(k+1<v.size())r[k+1]+=h;if(k)r[k-1]+=h;}return r;}
inline V polyapply(const V&p,const V&v){V r(v.size());for(auto it=p.rbegin();it!=p.rend();++it){r=mulx(r);if(!it->is_zero())for(unsigned j=0;j<v.size();++j)if(!v[j].is_zero())r[j]+=*it*v[j];}return r;}
inline V coeff(const Exact&p){if(p.variables()!=std::vector<std::string>{"x","eps","I"})throw std::invalid_argument("ultra unsupported field order");unsigned n=0;for(auto&t:p.numerator_terms())n=std::max(n,unsigned(t.powers[0]));V v(n+1);auto den=p.denominator().rational();for(auto&t:p.numerator_terms()){if(t.powers.size()!=3||t.powers[1]!=0)throw std::invalid_argument("coefficient contains unsupported field symbols");B b=B::from_strings((t.coefficient/den).str());for(unsigned k=0;k<t.powers[2]%4;++k)acb_mul_onei(b.raw(),b.raw());v[t.powers[0]]+=b;}return v;}

inline V power_to_cheb(const V&p){V r(p.size());for(auto it=p.rbegin();it!=p.rend();++it){V next(p.size());for(unsigned j=0;j<r.size();++j){if(r[j].is_zero())continue;if(j==0){if(next.size()>1)next[1]+=r[j];}else{B half;acb_mul_2exp_si(half.raw(),r[j].raw(),-1);next[j-1]+=half;if(j+1<next.size())next[j+1]+=half;}}next[0]+=*it;r=std::move(next);}return r;}
inline V multiplier_column(const V&cheb,unsigned column,unsigned rows){V out(rows);for(unsigned k=0;k<cheb.size();++k){if(cheb[k].is_zero())continue;if(k==0){if(column<rows)out[column]+=cheb[k];continue;}B half;acb_mul_2exp_si(half.raw(),cheb[k].raw(),-1);if(column+k<rows)out[column+k]+=half;int reflected=int(column)-int(k);if(reflected>=0){if(unsigned(reflected)<rows)out[reflected]+=half;}else if(reflected<=-2&&unsigned(-reflected-2)<rows)out[-reflected-2]-=half;}return out;}
inline M build(const std::vector<V>&q,const std::vector<std::vector<V>>&p,unsigned N){unsigned b=q.size();if(!b||b>8||!N||N>256||p.size()!=b)throw std::invalid_argument("ultra fast shape/order");std::vector<V> qc;std::vector<std::vector<V>> pc;for(unsigned i=0;i<b;++i){if(q[i].empty()||q[i].size()>33||p[i].size()!=b)throw std::invalid_argument("ultra fast degree/shape");qc.push_back(power_to_cheb(q[i]));pc.emplace_back();for(auto&v:p[i]){if(v.empty()||v.size()>33)throw std::invalid_argument("ultra fast degree");pc.back().push_back(power_to_cheb(v));}}
 M a(b*(N+1),V(b*(N+1)));for(unsigned n=0;n<=N;++n)for(unsigned j=0;j<b;++j){unsigned col=n?(n-1)*b+j:N*b+j;for(unsigned i=0;i<b;++i){auto v=multiplier_column(pc[i][j],n,N);if(n){if(n>=2){auto prev=multiplier_column(pc[i][j],n-2,N);for(unsigned k=0;k<N;++k)v[k]-=prev[k];}for(auto&c:v)acb_mul_2exp_si(c.raw(),c.raw(),-1);}if(i==j&&n){auto d=multiplier_column(qc[i],n-1,N);for(unsigned k=0;k<N;++k)a[k*b+i][col]=B(n)*d[k]-v[k];}else for(unsigned k=0;k<N;++k)a[k*b+i][col]=-v[k];}a[N*b+j][col]=B(n%2?-1:1);}return a;}

// Sparse LU of the almost-banded operator. Rows are exchanged only where the
// natural pivot is not provably nonzero, so a well-posed system whose
// diagonal entry cancels is still factored; each exchange is replayed on every
// right-hand side in elimination order.
struct Factor {
 M a; std::vector<std::vector<std::pair<unsigned,B>>> elimination;
 std::vector<std::set<unsigned>> cols;
 std::vector<unsigned> exchange;
 explicit Factor(M matrix):a(std::move(matrix)),elimination(a.size()),cols(a.size()),exchange(a.size()){
 unsigned z=a.size();std::vector<std::set<unsigned>> rows(z);
 for(unsigned i=0;i<z;++i)for(unsigned j=0;j<z;++j)if(!a[i][j].is_zero()){cols[i].insert(j);rows[j].insert(i);}
 const auto lower=[](const B&v){arf_t t;arf_init(t);acb_get_abs_lbound_arf(t,v.raw(),B::precision());const double x=arf_get_d(t,ARF_RND_DOWN);arf_clear(t);return x;};
 for(unsigned k=0;k<z;++k){exchange[k]=k;
 if(a[k][k].contains_zero()){unsigned best=k;
 for(auto i:rows[k])if(i>k&&!a[i][k].contains_zero()&&(best==k||lower(a[i][k])>lower(a[best][k])))best=i;
 if(best==k)throw std::runtime_error("ultra zero-containing pivot: whole-arm fallback required");
 for(auto j:cols[k])rows[j].erase(k);for(auto j:cols[best])rows[j].erase(best);
 std::swap(a[k],a[best]);std::swap(cols[k],cols[best]);
 for(auto j:cols[k])rows[j].insert(k);for(auto j:cols[best])rows[j].insert(best);
 exchange[k]=best;}
 std::vector<unsigned> targets;for(auto i:rows[k])if(i>k)targets.push_back(i);
 for(auto i:targets){B f=a[i][k]/a[k][k];elimination[k].push_back({i,f});
 for(auto j:cols[k])if(j>k){acb_submul(a[i][j].raw(),f.raw(),a[k][j].raw(),B::precision());cols[i].insert(j);rows[j].insert(i);}
 a[i][k]=B(0);cols[i].erase(k);}}
 }
 M solve(M rhs)const {unsigned z=a.size(),r=rhs.at(0).size();
 for(unsigned k=0;k<z;++k){if(exchange[k]!=k)std::swap(rhs[k],rhs[exchange[k]]);for(const auto&[i,f]:elimination[k])for(unsigned c=0;c<r;++c)acb_submul(rhs[i][c].raw(),f.raw(),rhs[k][c].raw(),B::precision());}
 M x(z,V(r));for(unsigned ii=z;ii-->0;)for(unsigned c=0;c<r;++c){B v=rhs[ii][c];for(auto j:cols[ii])if(j>ii)acb_submul(v.raw(),a[ii][j].raw(),x[j][c].raw(),B::precision());x[ii][c]=v/a[ii][ii];}return x;
 }
};
using EP=std::vector<V>;
inline EP bivariate(const Exact&p){
 if(p.variables()!=std::vector<std::string>{"x","eps","I"})throw std::invalid_argument("ultra unsupported field order");
 auto den=p.denominator().rational();unsigned dx=0,de=0;for(auto&t:p.numerator_terms()){dx=std::max(dx,unsigned(t.powers[0]));de=std::max(de,unsigned(t.powers[1]));}
 if(dx>32||de>64)throw std::invalid_argument("ultra bivariate degree budget");EP out(de+1,V(dx+1));
 for(auto&t:p.numerator_terms()){B v=B::from_strings((t.coefficient/den).str());for(unsigned i=0;i<t.powers[2]%4;++i)acb_mul_onei(v.raw(),v.raw());out[t.powers[1]][t.powers[0]]+=v;}return out;
}
using SparseRows=std::vector<std::vector<std::pair<unsigned,B>>>;
struct Multipliers {
 std::map<const V*,SparseRows> multiplier_cache;
 std::size_t multiplier_cells=0,multiplier_cap;
 explicit Multipliers(std::size_t cap):multiplier_cap(cap){}
 void add_product(M&rhs,unsigned row,unsigned size,const V&p,const M&source,bool subtract){
 if(std::all_of(p.begin(),p.end(),[](const B&v){return v.is_zero();}))return;
 unsigned r=source[0].size(),n=rhs.size()/size-1,lim=source.size();auto it=multiplier_cache.find(&p);
 if(it==multiplier_cache.end()){SparseRows matrix(n);for(unsigned j=0;j<lim;++j){V unit(lim);unit[j]=B(1);auto col=polyapply(p,unit);for(unsigned k=0;k<n;++k)if(!col[k].is_zero()){if(++multiplier_cells>multiplier_cap)throw std::runtime_error("ultra multiplication storage budget");matrix[k].push_back({j,col[k]});}}it=multiplier_cache.emplace(&p,std::move(matrix)).first;}
 for(unsigned k=0;k<n;++k)for(const auto&[j,v]:it->second[k])for(unsigned a=0;a<r;++a){if(subtract)acb_submul(rhs[k*size+row][a].raw(),v.raw(),source[j][a].raw(),B::precision());else acb_addmul(rhs[k*size+row][a].raw(),v.raw(),source[j][a].raw(),B::precision());}
}

};
}

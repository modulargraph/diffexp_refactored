#include "algebra.hpp"
#include "products.hpp"
void check(bool p,const char*s){if(!p)throw std::runtime_error(s);}
bool contains(const B&v,const Rational&q){fmpq_t r;fmpq_init(r);fmpq_set_str(r,q.str().c_str(),10);bool b=arb_contains_fmpq(acb_realref(v.raw()),r)&&arb_contains_zero(acb_imagref(v.raw()));fmpq_clear(r);return b;}
int main(){try{
 B::set_precision(256);ExactField field({"x","I"});
 std::vector<WeightedTerm> terms{{0,1,1,Exact(field,"1"),Rational(1)},{1,0,1,Exact(field,"1+x"),Rational(1)},{1,1,0,Exact(field,"1/(2*(1+x))"),Rational(1)}};
 auto groups=prepare_products(terms);unsigned N=14,K=3;
 Boundary initial(2,std::vector<B>(K+1,B(0)));initial[0][0]=initial[1][0]=B(1);
 std::vector<Boundary> jets;auto value=product_chart(groups,initial,B::from_strings("1/8"),N,&jets);
 std::vector<Rational> root(N+1,Rational(0));root[0]=Rational(1);for(unsigned n=1;n<=N;++n)root[n]=root[n-1]*(Rational("1/2")-Rational(n-1))/Rational(n);
 std::vector<std::vector<Rational>> y(N+1,std::vector<Rational>(K+1,Rational(0))),z=y;y[0][0]=Rational(1);
 for(unsigned n=0;n<N;++n)for(unsigned k=1;k<=K;++k)for(unsigned m=0;m<=n;++m)y[n+1][k]+=root[m]*y[n-m][k-1]/Rational(n+1);
 for(unsigned n=0;n<=N;++n)for(unsigned k=0;k<=K;++k)for(unsigned m=0;m<=n;++m)z[n][k]+=root[m]*y[n-m][k];
 unsigned checked=0;for(unsigned n=0;n<=N;++n)for(unsigned k=0;k<=K;++k){check(contains(jets[n][0][k],y[n][k]),"original algebraic coefficients");check(contains(jets[n][1][k],z[n][k]),"lifted product coefficients");checked+=2;}
 for(auto&row:initial)arb_add_error_2exp_si(acb_realref(row[0].raw()),-70);
 auto uncertain=product_chart(groups,initial,B::from_strings("1/8"),N);
 for(int sign:{-1,0,1})for(unsigned i=0;i<2;++i)for(unsigned k=0;k<=K;++k){Rational v(0);auto&a=i?z:y;for(unsigned n=N+1;n-->0;)v=v*Rational("1/8")+a[n][k];v=v*(Rational(1)+Rational(sign)*Rational("1/2361183241434822606848"));check(contains(uncertain[i][k],v),"carried uncertainty");}
 std::cout<<"{\"status\":\"passed\",\"exact_coefficients_checked\":"<<checked<<",\"boundary_realizations\":3}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}

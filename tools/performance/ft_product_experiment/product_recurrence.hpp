#pragma once
#include "diffexp/polynomial_transport.hpp"
#include <type_traits>

// Scratch experiment: individual rational-product recurrences, not an SCC
// implementation. Only epsilon-expanded real rational x coefficients at x=0.
// Computes a retained Taylor polynomial; it does NOT certify its omitted tail.
namespace experiment {
using namespace diffexp;
using B = Jet::Ball;
struct Edge { unsigned row, column, epsilon; std::vector<Rational> p,q; };
struct Prepared {
  std::vector<Edge> edges;
  polynomial_transport::Compiled baseline;
  bool supported = true;
};
inline std::vector<Rational> polynomial(const Exact& a) {
  if (!a.denominator().is_rational()) throw std::invalid_argument("not polynomial");
  std::vector<Rational> out(1, Rational(0));
  for (const auto& t : a.numerator_terms()) {
    unsigned degree=0;
    for (unsigned v=0;v<t.powers.size();++v) {
      if (a.variables()[v]=="x") degree=t.powers[v];
      else if (t.powers[v]) throw std::invalid_argument("unexpanded parameter");
    }
    if (degree>256) throw std::invalid_argument("degree budget");
    out.resize(std::max(out.size(),std::size_t(degree+1)),Rational(0));
    out[degree]+=t.coefficient/a.denominator().rational();
  }
  return out;
}
inline Prepared prepare(const std::vector<RationalLineEntry>& input,unsigned d,
                        unsigned order,unsigned high) {
  Prepared out;
  out.baseline=polynomial_transport::compile(input,d,order,high);
  try {
    for(const auto& e:input) {
      if(e.epsilon>high || e.coefficient.is_zero()) continue;
      Edge a{e.row,e.column,e.epsilon,polynomial(e.coefficient.numerator()),
             polynomial(e.coefficient.denominator())};
      if(a.q[0].is_zero()) throw std::invalid_argument("singular center");
      out.edges.push_back(std::move(a));
    }
  } catch(const std::invalid_argument&) {out.supported=false;out.edges.clear();}
  return out;
}
template<class T> T scalar(const Rational& q) {
  if constexpr(std::is_same_v<T,Rational>) return q;
  else return polynomial_transport::detail::ball(q);
}
template<class T> using Jets=std::vector<std::vector<std::vector<T>>>; // n,row,eps
// Independent oracle expands each p/q first and then uses the ordinary dense
// Taylor convolution. Exact rational mode checks all coefficients, not endpoints.
template<class T> Jets<T> coefficients(const Prepared& c,
    const std::vector<std::vector<T>>& boundary,unsigned order,bool dense=false,
    std::size_t* operations=nullptr) {
  const auto d=boundary.size(), width=boundary.at(0).size();
  if(!c.supported || !order || order>1000 || d!=c.baseline.dimension ||
     !width || width>c.baseline.epsilon_high+1)
    throw std::invalid_argument("unsupported prototype dimensions");
  for(const auto& r:boundary) if(r.size()!=width) throw std::invalid_argument("shape");
  // The product jets add E*N*width cells. Explicit finite resource budget.
  polynomial_transport::detail::cells({order+1,d+c.edges.size(),width},2000000);
  Jets<T> y(order+1,std::vector<std::vector<T>>(d,std::vector<T>(width,T(0))));
  y[0]=boundary;
  std::vector<std::vector<std::vector<T>>> f(c.edges.size(),
      std::vector<std::vector<T>>(order,std::vector<T>(width,T(0))));
  std::vector<std::vector<T>> ps,qs,as;
  for(const auto& e:c.edges) {
    std::vector<T> p,q;
    for(const auto& v:e.p)p.push_back(scalar<T>(v));
    for(const auto& v:e.q)q.push_back(scalar<T>(v));
    if constexpr(std::is_same_v<T,B>)
      if(!q[0].is_finite()||q[0].contains_zero())throw std::domain_error("pivot");
    ps.push_back(p);qs.push_back(q);
    if(dense) {
      std::vector<T> a(order,T(0));
      for(unsigned n=0;n<order;++n) {
        if(n<p.size())a[n]=p[n];
        for(unsigned m=1;m<q.size()&&m<=n;++m)a[n]-=q[m]*a[n-m];
        a[n]=a[n]/q[0];
      }
      as.push_back(std::move(a));
    }
  }
  std::size_t work=0;
  auto spend=[&]{if(++work>100000000)throw std::length_error("operation budget");};
  for(unsigned n=0;n<order;++n) {
    for(unsigned j=0;j<c.edges.size();++j) {
      const auto& e=c.edges[j];const auto& p=ps[j];const auto& q=qs[j];
      for(unsigned k=e.epsilon;k<width;++k) {
        auto& v=f[j][n][k-e.epsilon];
        if(dense) {
          for(unsigned m=0;m<=n;++m){spend();v+=as[j][m]*y[n-m][e.column][k-e.epsilon];}
        } else {
          for(unsigned m=0;m<p.size()&&m<=n;++m){spend();v+=p[m]*y[n-m][e.column][k-e.epsilon];}
          for(unsigned m=1;m<q.size()&&m<=n;++m){spend();v-=q[m]*f[j][n-m][k-e.epsilon];}
          v=v/q[0];
        }
        y[n+1][e.row][k]+=v;
      }
    }
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<width;++k)
      y[n+1][i][k]=y[n+1][i][k]/T(n+1);
  }
  if(operations)*operations=work;
  return y;
}
template<class T> std::vector<std::vector<T>> horner(const Jets<T>& y,const T& h) {
  auto out=y[0];for(auto& r:out)for(auto& v:r)v=T(0);
  for(std::size_t n=y.size();n-->0;)
    for(unsigned i=0;i<out.size();++i)for(unsigned k=0;k<out[i].size();++k)
      out[i][k]=out[i][k]*h+y[n][i][k];
  return out;
}
inline Boundary chart(const Prepared& c,const Boundary& b,const B& center,
                      const B& step,unsigned order,bool* fallback=nullptr) {
  for(const auto& r:b)for(const auto& v:r)
    if(!v.is_finite())throw std::invalid_argument("nonfinite boundary");
  if(!step.is_finite()||!center.is_finite())throw std::invalid_argument("nonfinite point");
  if(!c.supported || !center.is_zero() ||
     (std::size_t(order)+1)*(b.size()+c.edges.size())*(c.baseline.epsilon_high+1)>2000000) {
    if(fallback)*fallback=true;
    return polynomial_transport::chart(c.baseline,b,center,step,order);
  }
  if(fallback)*fallback=false;
  return horner(coefficients<B>(c,b,order),step);
}
} // namespace experiment

// Experiment-only copy of adjoint_detail::chart; only Horner temporaries changed.
#pragma once
#include "diffexp/adjoint_transport.hpp"
namespace diffexp::adjoint_detail {
inline Boundary inplace_horner_chart(const std::vector<Entry>& entries,const Boundary& initial,const Jet::Ball& center,
    const Jet::Ball& step,unsigned order) {
  using B=Jet::Ball;const unsigned d=initial.size(),kmax=initial[0].size()-1;
  Jet x(0,order+1,B::precision());x.set(0,center);x.set(1,B(1));
  auto imaginary=x.constant(0);imaginary.set(0,B::from_strings("0","1"));
  std::vector<std::vector<B>> coefficients;
  for(const auto& entry:entries) {
    auto jet=evaluate(entry.coefficient,x,{{"x",x},{"I",imaginary}});coefficients.emplace_back();
    for(unsigned n=0;n<order;++n)coefficients.back().push_back(jet.at(n));
  }
  std::vector<B> values(static_cast<std::size_t>(order+1)*d*(kmax+1),B(0));
  const auto at=[&](unsigned n,unsigned i,unsigned k)->B&{return values[(static_cast<std::size_t>(n)*d+i)*(kmax+1)+k];};
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)at(0,i,k)=initial[i][k];
  for(unsigned n=0;n<order;++n) {
    for(unsigned e=0;e<entries.size();++e)for(const auto& [row,column]:entries[e].positions)
      for(unsigned k=entries[e].epsilon;k<=kmax;++k)for(unsigned m=0;m<=n;++m)
        if(!coefficients[e][m].is_zero() && !at(n-m,column,k-entries[e].epsilon).is_zero())
          acb_addmul(at(n+1,row,k).raw(),coefficients[e][m].raw(),at(n-m,column,k-entries[e].epsilon).raw(),B::precision());
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)acb_div_ui(at(n+1,i,k).raw(),at(n+1,i,k).raw(),n+1,B::precision());
  }
  Boundary out(d,std::vector<B>(kmax+1,B(0)));
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)for(unsigned n=order+1;n-->0;) {
    acb_mul(out[i][k].raw(),out[i][k].raw(),step.raw(),B::precision());
    acb_add(out[i][k].raw(),out[i][k].raw(),at(n,i,k).raw(),B::precision());
  }
  return out;
}
}

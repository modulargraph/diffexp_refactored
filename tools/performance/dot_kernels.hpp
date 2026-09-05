// Experiment kernels copied from tools/benchmark_adjoint_dot.cpp.
#pragma once
// Opt-in experiment only. No production recurrence or conditioning changes.
#include "diffexp/adjoint_transport.hpp"
#include <chrono>
#include <ctime>
#include <iostream>
#include <memory>
using namespace diffexp;
using B=Jet::Ball;
struct AcbArray {
  slong size;acb_ptr values;
  explicit AcbArray(slong n):size(n),values(_acb_vec_init(n)){}
  ~AcbArray(){_acb_vec_clear(values,size);}
  AcbArray(const AcbArray&)=delete;
};
Boundary experimental_chart(const std::vector<adjoint_detail::Entry>& entries,const Boundary& initial,
    const B& center,const B& step,unsigned order,unsigned grouped) {
  const unsigned d=initial.size(),width=initial[0].size();const auto bits=B::precision();
  Jet x(0,order+1,bits);x.set(0,center);x.set(1,B(1));
  auto imaginary=x.constant(0);imaginary.set(0,B::from_strings("0","1"));
  AcbArray coefficients(static_cast<slong>(entries.size())*order);
  for(unsigned e=0;e<entries.size();++e){auto jet=evaluate(entries[e].coefficient,x,{{"x",x},{"I",imaginary}});
    for(unsigned n=0;n<order;++n){auto value=jet.at(n);acb_swap(coefficients.values+e*order+n,value.raw());}}
  const slong stride=static_cast<slong>(d)*width;
  AcbArray values(static_cast<slong>(order+1)*stride);
  const auto at=[&](unsigned n,unsigned i,unsigned k){return values.values+static_cast<slong>(n)*stride+i*width+k;};
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<width;++k)acb_set(at(0,i,k),initial[i][k].raw());
  B sum;
  for(unsigned n=0;n<order;++n){
    for(unsigned e=0;e<entries.size();++e)for(const auto& [row,column]:entries[e].positions)
      for(unsigned k=entries[e].epsilon;k<width;++k){
        auto target=at(n+1,row,k);auto coefficient=coefficients.values+e*order;
        if(grouped){
          // Real Acb arrays permit the documented signed stride. This does
          // not reinterpret C++ wrapper objects as C arrays or copy the lag
          // vectors into temporary gather buffers.
          if(grouped==2)acb_dot_precise(sum.raw(),target,0,coefficient,1,at(n,column,k-entries[e].epsilon),-stride,n+1,bits);
          else acb_dot(sum.raw(),target,0,coefficient,1,at(n,column,k-entries[e].epsilon),-stride,n+1,bits);
          acb_swap(target,sum.raw());
        }else for(unsigned m=0;m<=n;++m)
          if(!acb_is_zero(coefficient+m) && !acb_is_zero(at(n-m,column,k-entries[e].epsilon)))
            acb_addmul(target,coefficient+m,at(n-m,column,k-entries[e].epsilon),bits);
      }
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<width;++k)acb_div_ui(at(n+1,i,k),at(n+1,i,k),n+1,bits);
  }
  Boundary out(d,std::vector<B>(width,B(0)));
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<width;++k){
    acb_zero(out[i][k].raw());
    for(unsigned n=order+1;n-->0;){
      acb_mul(out[i][k].raw(),out[i][k].raw(),step.raw(),bits);
      acb_add(out[i][k].raw(),out[i][k].raw(),at(n,i,k),bits);
    }
  }
  return out;
}
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
struct Timing {double wall,cpu;};
template<class F>Timing time_charts(unsigned repeats,F run){
  auto start=std::chrono::steady_clock::now();auto cpu=std::clock();
  for(unsigned i=0;i<repeats;++i){auto out=run();require(out.back().back().is_finite(),"nonfinite timed result");}
  return {std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/repeats,double(std::clock()-cpu)/CLOCKS_PER_SEC/repeats};
}

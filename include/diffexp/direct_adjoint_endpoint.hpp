#pragma once
#include "diffexp/affine_frobenius.hpp"
#include "diffexp/exact_jet_composition.hpp"

namespace diffexp::direct_adjoint_endpoint {
using Matrix=AffineFrobeniusSeries::Matrix;
struct Options {
  unsigned max_dimension=256,max_rows=5000,max_x_order=256,max_coefficient_matrices=512;
  std::size_t max_entries=2000000;
};
enum class Status { Ready, Unsupported };
struct Result {
  Status status=Status::Unsupported;
  std::string reason;
  // q(x)=sum_{k=0}^{coefficients.size()-1} coefficients[k]*x^(first_power+k).
  // q'+q*A=w through x^(last_power-1), where last_power is the last retained q power.
  long first_power=0;
  std::vector<Matrix> coefficients;
  std::optional<long> epsilon_low;
  bool success()const{return status==Status::Ready;}
};
// If true, every rationally projected endpoint-constant row vanishes under
// the formal DR prescription and has no fixed-sector domain constraints.
inline bool only_epsilon_dependent_exponents(const Matrix& a,std::size_t xi,std::size_t ei) {
  if(a.empty()||a.front().empty()||xi==ei)throw std::invalid_argument("direct adjoint endpoint spectrum dimensions");
  const auto d=a.size();
  if(xi>=a[0][0].variables().size()||ei>=a[0][0].variables().size())throw std::invalid_argument("direct adjoint endpoint variable index");
  auto z=a[0][0].constant(0),x=z.variable(xi);
  auto residue=fuchsify::detail::zeros(d,d,z);
  for(unsigned i=0;i<d;++i){
    if(a[i].size()!=d)throw std::invalid_argument("direct adjoint endpoint spectrum shape");
    for(unsigned j=0;j<d;++j){
      z.require_same_field(a[i][j]);
      if(!a[i][j].is_zero()&&fuchsify::detail::valuation(a[i][j],xi)<-1)return false;
      residue[i][j]=affine_frobenius_detail::taylor(x*a[i][j],xi,0)[0];
    }
  }
  const auto spectrum=affine_frobenius_detail::spectrum(residue,xi,ei);
  return std::all_of(spectrum.begin(),spectrum.end(),[](const auto& e){return !e.slope.is_zero();});
}
// Universal formal-DR primitive in the supplied Fuchsian frame. No physical
// boundary is read and no coefficient of a nonanalytic physical mode is set to
// zero. The guard proves every fixed sector of q*F has strictly positive power.
// Nonzero epsilon slopes use the same meromorphic termwise DR prescription as
// AffineFrobeniusSeries; no common convergence strip or tail bound is claimed.
inline Result prepare(const Matrix& a,const Matrix& w,std::size_t xi,std::size_t ei,
                      unsigned last_power,const Options& options={}) {
  using namespace fuchsify::detail;
  using affine_frobenius_detail::taylor;
  if(a.empty()||a.front().empty()||a.size()>options.max_dimension||w.empty()||w.size()>options.max_rows||
     last_power>options.max_x_order||xi==ei||!options.max_coefficient_matrices||!options.max_entries)
    throw std::invalid_argument("direct adjoint endpoint dimensions or budgets");
  const unsigned d=a.size(),r=w.size();
  if(xi>=a[0][0].variables().size()||ei>=a[0][0].variables().size())throw std::invalid_argument("direct adjoint endpoint variable index");
  auto z=a[0][0].constant(0),x=z.variable(xi);
  for(const auto& row:a)if(row.size()!=d)throw std::invalid_argument("direct adjoint endpoint connection shape");
  for(const auto& row:w)if(row.size()!=d)throw std::invalid_argument("direct adjoint endpoint observable shape");
  for(const auto& row:a)for(const auto& c:row)z.require_same_field(c);
  for(const auto& row:w)for(const auto& c:row)z.require_same_field(c);
  Result out;std::optional<long> valuation_w;
  for(const auto& row:w)for(const auto& c:row)if(!c.is_zero()){
    auto v=valuation(c,xi);valuation_w=valuation_w?std::min(*valuation_w,v):v;
  }
  if(std::size_t(r)*d>options.max_entries){out.reason="direct adjoint endpoint aggregate entry budget";return out;}
  if(!valuation_w){out.status=Status::Ready;out.first_power=0;out.coefficients.push_back(zeros(r,d,z));return out;}
  out.first_power=*valuation_w+1;
  if(out.first_power>long(last_power)){out.reason="requested endpoint order precedes the observable primitive";return out;}
  const long depth=long(last_power)-out.first_power;
  if(depth>=long(options.max_coefficient_matrices)){out.reason="direct adjoint endpoint coefficient budget";return out;}
  const std::size_t entries_per_order=std::size_t(d)*d+2*std::size_t(r)*d;
  const std::size_t scratch_entries=3*std::size_t(d)*d+3*std::size_t(r)*d;
  if(scratch_entries>options.max_entries || std::size_t(depth+1)>(options.max_entries-scratch_entries)/entries_per_order){
    out.reason="direct adjoint endpoint aggregate entry budget";return out;
  }
  auto residue=zeros(d,d,z);
  for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j){
    if(!a[i][j].is_zero()&&valuation(a[i][j],xi)<-1){out.reason="connection is not Fuchsian";return out;}
    residue[i][j]=taylor(x*a[i][j],xi,0)[0];
  }
  // The exact characteristic factors suffice; no eigenvectors or tails are built.
  std::vector<affine_frobenius_detail::Exponent> spectrum;
  try{spectrum=affine_frobenius_detail::spectrum(residue,xi,ei);}catch(const std::domain_error& error){out.reason=error.what();return out;}
  for(const auto& exponent:spectrum)
    if(exponent.slope.is_zero() && exponent.power+Rational(out.first_power)<=Rational(0)){
      out.reason="fixed endpoint sectors require the full physical-domain constraints";return out;
    }
  std::vector<Matrix> g(depth+1,zeros(d,d,z)),forcing(depth+1,zeros(r,d,z));g[0]=residue;
  for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j){
    const auto c=taylor(x*a[i][j],xi,depth);for(long k=1;k<=depth;++k)g[k][i][j]=c[k];
  }
  for(unsigned i=0;i<r;++i)for(unsigned j=0;j<d;++j){
    const auto c=taylor(w[i][j]/power(x,*valuation_w),xi,depth);for(long k=0;k<=depth;++k)forcing[k][i][j]=c[k];
  }
  for(long offset=0;offset<=depth;++offset){
    const long n=out.first_power+offset;auto rhs=forcing[offset];
    for(long k=0;k<offset;++k){const auto p=multiply(out.coefficients[k],g[offset-k]);
      for(unsigned i=0;i<r;++i)for(unsigned j=0;j<d;++j)rhs[i][j]=rhs[i][j]-p[i][j];
    }
    auto op=residue;for(unsigned i=0;i<d;++i)op[i][i]=op[i][i]+z.constant(n);
    Matrix inv;
    try{inv=inverse(op);}catch(const std::domain_error&){
      out.reason="identically resonant adjoint recurrence at x power "+std::to_string(n);
      out.coefficients.clear();out.epsilon_low.reset();return out;
    }
    auto q=multiply(rhs,inv);
    for(const auto& row:q)for(const auto& c:row)if(!c.is_zero()){
      auto low=valuation(c,ei);out.epsilon_low=out.epsilon_low?std::min(*out.epsilon_low,low):low;
    }
    out.coefficients.push_back(std::move(q));
  }
  out.status=Status::Ready;return out;
}
inline Matrix polynomial(const Result& result,std::size_t xi){
  if(!result.success()||result.coefficients.empty())throw std::invalid_argument("direct adjoint endpoint result is not ready");
  using namespace fuchsify::detail;
  auto out=result.coefficients.front();auto z=out[0][0].constant(0),x=z.variable(xi);
  for(auto& row:out)for(auto& c:row)c=z;
  for(unsigned k=0;k<result.coefficients.size();++k)
    for(unsigned i=0;i<out.size();++i)for(unsigned j=0;j<out[i].size();++j)
      out[i][j]=out[i][j]+result.coefficients[k][i][j]*power(x,result.first_power+k);
  return out;
}
// Exact coefficients in a finite spatial series do not make its primitive
// exact. This additional identity is required for full-functional rewriting.
inline bool exact_residual_zero(const Matrix& primitive,const Matrix& connection,
    const Matrix& weight,std::size_t xi,std::size_t ei) {
  const auto product=exact_jet_composition::fuse(primitive,connection,ei).product;
  if(weight.size()!=primitive.size())throw std::invalid_argument("exact primitive weight shape");
  for(std::size_t i=0;i<primitive.size();++i) {
    if(weight[i].size()!=primitive[i].size())throw std::invalid_argument("exact primitive weight width");
    for(std::size_t j=0;j<primitive[i].size();++j)
      if(!(primitive[i][j].derivative(xi)+product[i][j]-weight[i][j]).is_zero())return false;
  }
  return true;
}
} // namespace diffexp::direct_adjoint_endpoint

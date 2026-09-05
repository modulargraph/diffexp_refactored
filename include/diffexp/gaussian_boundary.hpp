#pragma once
#include "diffexp/ibp.hpp"
#include "diffexp/deepest_beta.hpp"

namespace diffexp::gaussian {
struct Budget {std::size_t max_degree=8,max_terms=10000,max_pairings=100000,max_products=100000,max_power=10000;};
struct Boundary {
  Exact U,F,multiplier,dimension;
  unsigned power,loops;
};
// The complete integral is multiplier(d,kinematics) times the scalar tadpole
// (U,F,power,loops). U,F and multiplier stay exact, including a symbolic current
// Feynman parameter. No endpoint substitution is performed here.
inline Boundary reduce(const ibp::PropagatorBasis& basis,const ibp::Integral& target,
    const Exact& dimension,const Budget& budget={},const std::vector<ibp::Affine>& additional_numerators={}) {
  if(basis.physical_count!=1||target.size()!=basis.denominators.size())throw std::invalid_argument("Gaussian boundary needs one physical propagator and complete indices");
  for(std::size_t k=1;k<target.size();++k)if(target[k]>0)throw std::invalid_argument("positive Gaussian auxiliary index");
  if(!budget.max_terms||!budget.max_pairings||!budget.max_products)throw std::invalid_argument("Gaussian expansion requires finite positive budgets");
  const auto zero=dimension.constant(0),one=dimension.constant(1),half=dimension.constant(Rational("1/2"));
  (void)(zero+basis.space.zero().constant);
  const auto L=basis.space.loops,E=basis.space.externals();
  if(target[0]<=0)return {one,one,zero,dimension,0,L}; // Polynomial loop integrals are scaleless.
  if(unsigned(target[0])>budget.max_power)throw std::runtime_error("Gaussian propagator power budget exceeded");
  const auto& denominator=basis.denominators[0];
  std::vector<std::vector<Exact>> Q(L,std::vector<Exact>(L,zero)),R(L,std::vector<Exact>(E,zero)),inverse=Q;
  for(std::size_t k=0;k<basis.space.pairs.size();++k) {
    auto [i,j]=basis.space.pairs[k];auto c=denominator.linear[k];
    if(j<L){Q[i][j]=-c*(i==j?one:half);Q[j][i]=Q[i][j];}
    else R[i][j-L]=-c*half;
  }
  for(unsigned i=0;i<L;++i)inverse[i][i]=one;
  auto det=one;
  for(unsigned i=0;i<L;++i) {
    auto pivot=i;while(pivot<L&&Q[pivot][i].is_zero())++pivot;
    if(pivot==L)throw std::domain_error("singular deepest loop quadratic form");
    if(pivot!=i){std::swap(Q[pivot],Q[i]);std::swap(inverse[pivot],inverse[i]);det=-det;}
    auto c=Q[i][i];det=det*c;
    for(unsigned j=0;j<L;++j){Q[i][j]=Q[i][j]/c;inverse[i][j]=inverse[i][j]/c;}
    for(unsigned k=0;k<L;++k)if(k!=i){auto factor=Q[k][i];for(unsigned j=0;j<L;++j){Q[k][j]=Q[k][j]-factor*Q[i][j];inverse[k][j]=inverse[k][j]-factor*inverse[i][j];}}
  }
  auto h=R;for(unsigned i=0;i<L;++i)for(unsigned e=0;e<E;++e){h[i][e]=zero;for(unsigned j=0;j<L;++j)h[i][e]=h[i][e]+inverse[i][j]*R[j][e];}
  auto M=denominator.constant;
  for(unsigned i=0;i<L;++i)for(unsigned e=0;e<E;++e)for(unsigned f=0;f<E;++f)M=M+R[i][e]*h[i][f]*basis.space.external_gram[e][f];
  if(M.is_zero())return {det,zero,zero,dimension,unsigned(target[0]),L}; // Massless vacuum after translation.
  auto pair_index=[&](unsigned i,unsigned j){if(i>j)std::swap(i,j);auto p=std::find(basis.space.pairs.begin(),basis.space.pairs.end(),std::pair{i,j});if(p==basis.space.pairs.end())throw std::logic_error("Gaussian pair shape");return p-basis.space.pairs.begin();};
  std::vector<ibp::Affine> shifted;
  for(const auto& [i,j]:basis.space.pairs) {
    auto s=basis.space.zero();s.linear[pair_index(i,j)]=one;
    if(j<L) {
      for(unsigned e=0;e<E;++e){s.linear[pair_index(i,L+e)]=s.linear[pair_index(i,L+e)]-h[j][e];s.linear[pair_index(j,L+e)]=s.linear[pair_index(j,L+e)]-h[i][e];}
      for(unsigned e=0;e<E;++e)for(unsigned f=0;f<E;++f)s.constant=s.constant+h[i][e]*h[j][f]*basis.space.external_gram[e][f];
    }else for(unsigned e=0;e<E;++e)s.constant=s.constant-h[i][e]*basis.space.external_gram[e][j-L];
    shifted.push_back(std::move(s));
  }
  std::size_t degree=additional_numerators.size(),products=0,pairings=0;
  if(degree>budget.max_degree)throw std::runtime_error("Gaussian numerator degree budget exceeded");
  for(std::size_t k=1;k<target.size();++k)if(target[k]<0){auto p=std::size_t(-std::int64_t(target[k]));if(p>budget.max_degree-degree)throw std::runtime_error("Gaussian numerator degree budget exceeded");degree+=p;}
  using Monomial=std::vector<unsigned>;std::map<Monomial,Exact> polynomial;
  polynomial.emplace(Monomial(basis.space.size(),0),one);
  std::vector<std::pair<ibp::Affine,std::int64_t>> factors;
  for(std::size_t k=1;k<target.size();++k)if(target[k]<0)factors.emplace_back(basis.denominators[k],-std::int64_t(target[k]));
  for(const auto& a:additional_numerators){if(a.linear.size()!=basis.space.size())throw std::invalid_argument("Gaussian extra numerator shape");factors.emplace_back(a,1);}
  for(const auto& [old,exponent]:factors) {
    auto centered=basis.space.zero();centered.constant=old.constant;
    for(std::size_t j=0;j<shifted.size();++j)centered=ibp::plus(centered,ibp::scaled(shifted[j],old.linear[j]));
    for(std::int64_t power=0;power<exponent;++power) {
      std::map<Monomial,Exact> next;
      auto insert=[&](const Monomial& m,const Exact& c){if(c.is_zero())return;if(++products>budget.max_products)throw std::runtime_error("Gaussian multiplication budget exceeded");ibp::add(next,m,c);if(next.size()>budget.max_terms)throw std::runtime_error("Gaussian monomial budget exceeded");};
      for(const auto& [m,c]:polynomial){insert(m,c*centered.constant);for(std::size_t j=0;j<m.size();++j)if(!centered.linear[j].is_zero()){auto p=m;++p[j];insert(p,c*centered.linear[j]);}}
      polynomial=std::move(next);
    }
  }
  auto multiplier=zero;
  for(const auto& [monomial,coefficient]:polynomial) {
    std::vector<unsigned> types,loop_vertices;
    for(std::size_t j=0;j<monomial.size();++j)for(unsigned k=0;k<monomial[j];++k){auto [a,b]=basis.space.pairs[j];types.push_back(a);types.push_back(b);}
    for(unsigned j=0;j<types.size();++j)if(types[j]<L)loop_vertices.push_back(j);
    if(loop_vertices.size()%2)continue;
    const unsigned pairs=loop_vertices.size()/2;auto moment=zero;
    std::vector<std::pair<unsigned,unsigned>> matching;
    std::function<void(std::vector<unsigned>,Exact)> wick=[&](std::vector<unsigned> remaining,Exact weight) {
      if(remaining.empty()) {
        if(++pairings>budget.max_pairings)throw std::runtime_error("Gaussian Wick pairing budget exceeded");
        std::vector<unsigned> parent(types.size());std::iota(parent.begin(),parent.end(),0);
        auto root=[&](unsigned v){while(parent[v]!=v)v=parent[v];return v;};
        auto join=[&](unsigned a,unsigned b){parent[root(a)]=root(b);};
        for(unsigned k=0;k<types.size();k+=2)join(k,k+1);for(auto [a,b]:matching)join(a,b);
        std::map<unsigned,std::vector<unsigned>> components;
        for(unsigned k=0;k<types.size();++k){auto& external=components[root(k)];if(types[k]>=L)external.push_back(types[k]-L);}
        for(const auto& [r,external]:components){if(external.empty())weight=weight*dimension;else if(external.size()==2)weight=weight*basis.space.external_gram[external[0]][external[1]];else throw std::logic_error("Gaussian Wick contraction graph");}
        moment=moment+weight;return;
      }
      const auto first=remaining.front();
      for(std::size_t j=1;j<remaining.size();++j) {
        auto covariance=-half*inverse[types[first]][types[remaining[j]]];if(covariance.is_zero())continue;
        auto next=remaining;next.erase(next.begin()+j);next.erase(next.begin());matching.emplace_back(first,remaining[j]);wick(std::move(next),weight*covariance);matching.pop_back();
      }
    };
    wick(loop_vertices,one);
    auto radial=one;for(unsigned r=1;r<=pairs;++r)radial=radial*M/(dimension.constant(target[0])-dimension.constant(L)*dimension*half-dimension.constant(r));
    multiplier=multiplier+coefficient*moment*radial;
  }
  return {det,det*M,multiplier,dimension,unsigned(target[0]),L};
}
inline Boundary specialize(const Boundary& b,std::span<const Exact> point) {
  if(!(b.dimension.substitute(point)==b.dimension))throw std::invalid_argument("Gaussian specialization cannot fix integration dimension");
  return {b.U.substitute(point),b.F.substitute(point),b.multiplier.substitute(point),b.dimension,b.power,b.loops};
}
// Geometry must be rationally anchored first. The multiplier's exact epsilon
// valuation is factored before asking the scalar tadpole for lookahead orders.
inline feynman::LaurentValue evaluate(const Boundary& b,std::size_t dimension_variable,int d0,
    unsigned epsilon_top,int rim=1) {
  using B=Jet::Ball;
  if(epsilon_top>100||b.power>10000||b.loops>64||d0< -64||d0>64||dimension_variable>=b.dimension.variable_count()||!(b.dimension==b.dimension.variable(dimension_variable)))throw std::invalid_argument("Gaussian Laurent evaluator needs a dimension variable and bounded order");
  if(b.multiplier.is_zero())return {0,std::vector<B>(epsilon_top+1,B(0)),true};
  auto eps=b.dimension.variable(dimension_variable);std::vector<Exact> point;
  for(std::size_t i=0;i<b.dimension.variable_count();++i)point.push_back(b.dimension.variable(i));
  point[dimension_variable]=b.dimension.constant(d0)-b.dimension.constant(2)*eps;
  auto coefficient=b.multiplier.substitute(point);auto origin=feynman::factor_origin(coefficient,dimension_variable);
  if(origin.power< -100||origin.power>100)throw std::runtime_error("Gaussian epsilon valuation budget exceeded");
  const unsigned lookahead=epsilon_top+unsigned(std::max(0L,-origin.power))+2;
  auto scalar=feynman::tadpole(b.U.rational(),b.F.rational(),b.power,b.loops,d0,lookahead,rim);
  const int low=scalar.epsilon_low+origin.power;
  if(low>int(epsilon_top))return {0,std::vector<B>(epsilon_top+1,B(0)),true};
  Jet context(0,epsilon_top-low+1,B::precision());auto variable=context; if(variable.length()>1)variable.set(1,B(1));
  auto regular=diffexp::evaluate(data::Reader(origin.regular.str()).read(),context,{{b.dimension.variables()[dimension_variable],variable}});
  std::vector<B> coefficients;
  for(int order=low;order<=int(epsilon_top);++order){B sum(0);for(int k=0;k<=order-low;++k)sum+=regular.at(k)*scalar.coefficients.at(order-low-k);coefficients.push_back(sum);}
  return {low,std::move(coefficients),true};
}
} // namespace diffexp::gaussian

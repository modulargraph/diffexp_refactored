#pragma once
#include "diffexp/exact.hpp"
#include <climits>
#include <map>

namespace diffexp::feynman {
struct MomentumLine {
  std::vector<Rational> loop_coefficients,external_coefficients;
  Rational mass_squared;
};
struct Family {
  unsigned loops;
  std::vector<std::vector<Rational>> external_gram;
  std::vector<MomentumLine> lines;
  void validate()const {
    if(!loops || lines.empty())throw std::invalid_argument("empty integral family");
    const auto n=external_gram.size();
    for(const auto& row:external_gram)if(row.size()!=n)throw std::invalid_argument("external Gram shape");
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)
      if(!(external_gram[i][j]==external_gram[j][i]))throw std::invalid_argument("external Gram is not symmetric");
    for(const auto& line:lines)
      if(line.loop_coefficients.size()!=loops || line.external_coefficients.size()!=n)
        throw std::invalid_argument("propagator momentum shape");
  }
};
struct Symanzik {Exact U,F;};

// Complete the square over loop momenta using their bilinear form. External
// scalar products are supplied explicitly; vector multiplication is never
// interpreted as multiplication of ordinary commutative variables.
inline Symanzik symanzik(const Family& family,const std::vector<Exact>& weights) {
  family.validate();
  if(weights.size()!=family.lines.size())throw std::invalid_argument("Symanzik weight count");
  const auto zero=weights[0].constant(0),one=weights[0].constant(1);
  const auto L=family.loops;
  const auto E=family.external_gram.size();
  std::vector<std::vector<Exact>> Q(L,std::vector<Exact>(L,zero)),R(L,std::vector<Exact>(E,zero));
  auto C=zero;
  for(std::size_t n=0;n<weights.size();++n) {
    const auto& line=family.lines[n];const auto& w=weights[n];
    C=C+w*w.constant(line.mass_squared);
    for(unsigned i=0;i<L;++i) {
      for(unsigned j=0;j<L;++j)Q[i][j]=Q[i][j]+w*w.constant(line.loop_coefficients[i]*line.loop_coefficients[j]);
      for(std::size_t j=0;j<E;++j)R[i][j]=R[i][j]+w*w.constant(line.loop_coefficients[i]*line.external_coefficients[j]);
    }
    for(std::size_t i=0;i<E;++i)for(std::size_t j=0;j<E;++j)
      C=C-w*w.constant(line.external_coefficients[i]*family.external_gram[i][j]*line.external_coefficients[j]);
  }
  auto U=one;
  // Symmetric Schur complements preserve the exact quadratic form and require
  // only L pivots. Generic zero leading pivots are handled by row/column swaps.
  for(unsigned i=0;i<L;++i) {
    if(Q[i][i].is_zero()) {
      unsigned j=i+1;while(j<L && Q[j][j].is_zero())++j;
      if(j==L)throw std::domain_error("degenerate quadratic family: no symmetric pivot");
      std::swap(Q[i],Q[j]);for(auto& row:Q)std::swap(row[i],row[j]);std::swap(R[i],R[j]);
    }
    auto pivot=Q[i][i];U=U*pivot;
    for(std::size_t e=0;e<E;++e)for(std::size_t f=0;f<E;++f)
      C=C+R[i][e]*R[i][f]*zero.constant(family.external_gram[e][f])/pivot;
    for(unsigned j=i+1;j<L;++j) {
      for(std::size_t e=0;e<E;++e)R[j][e]=R[j][e]-Q[j][i]*R[i][e]/pivot;
      for(unsigned k=i+1;k<L;++k)Q[j][k]=Q[j][k]-Q[j][i]*Q[i][k]/pivot;
    }
  }
  return {U,U*C};
}

enum class Operation {BetaIntegral,LowerLimit,UpperLimit,Direct};
struct MergeRequest {
  Operation operation;
  std::vector<int> source_indices;
  unsigned left_power=0,right_power=0;
  Rational normalization{1};
};
inline MergeRequest merge_request(std::vector<int> indices,unsigned left,unsigned right) {
  if(left==right || left>=indices.size() || right>=indices.size())throw std::invalid_argument("invalid merge positions");
  int a=indices[left],b=indices[right];
  if(a<0 || b<0)throw std::invalid_argument("reduce physical denominator numerators before merging");
  if(a>INT_MAX-b)throw std::overflow_error("merged propagator power");
  indices[left]=a+b;indices[right]=0;
  if(!a && !b)return {Operation::Direct,std::move(indices)};
  if(!a)return {Operation::LowerLimit,std::move(indices)};
  if(!b)return {Operation::UpperLimit,std::move(indices)};
  Rational factor(1);
  // Gamma(a+b)/(Gamma(a) Gamma(b)), exact for positive integral powers.
  for(int k=1;k<=a+b-1;++k)factor*=Rational(k);
  for(int k=1;k<a;++k)factor=factor/Rational(k);
  for(int k=1;k<b;++k)factor=factor/Rational(k);
  return {Operation::BetaIntegral,std::move(indices),static_cast<unsigned>(a-1),static_cast<unsigned>(b-1),factor};
}

inline std::vector<Exact> merge_weights(const Exact& sample,unsigned lines,
    const std::vector<std::pair<unsigned,unsigned>>& sequence,const std::vector<Exact>& parameters) {
  if(!lines || sequence.size()!=parameters.size())throw std::invalid_argument("merge sequence/parameter mismatch");
  auto zero=sample.constant(0),one=sample.constant(1);
  std::vector<std::vector<Exact>> groups(lines,std::vector<Exact>(lines,zero));
  std::vector<bool> active(lines,true);
  for(unsigned i=0;i<lines;++i)groups[i][i]=one;
  for(std::size_t step=0;step<sequence.size();++step) {
    auto [i,j]=sequence[step];
    if(i==j || i>=lines || j>=lines || !active[i] || !active[j])throw std::invalid_argument("merge sequence reuses an eliminated line");
    const auto& x=parameters[step];
    for(unsigned k=0;k<lines;++k)groups[i][k]=x*groups[i][k]+(one-x)*groups[j][k];
    active[j]=false;
  }
  if(std::count(active.begin(),active.end(),true)!=1)throw std::invalid_argument("deepest boundary requires a complete merge tree");
  return groups[std::find(active.begin(),active.end(),true)-active.begin()];
}

// A row's exact valuation and its primitive's actual epsilon pole depth are
// separate inputs. A log^p / x sector can lose p+1 orders, not always one.
inline int source_epsilon_top(int output_top,int row_valuation,unsigned primitive_pole_depth) {
  auto n=static_cast<std::int64_t>(output_top)-row_valuation+primitive_pole_depth;
  if(n<INT_MIN || n>INT_MAX)throw std::overflow_error("Feynman-trick epsilon demand");
  return static_cast<int>(n);
}

inline Family banana(unsigned loops,const std::vector<Rational>& masses_squared) {
  if(!loops || masses_squared.size()!=loops+1)throw std::invalid_argument("banana mass count");
  Family f{loops,{{Rational(-1)}},{}};
  for(unsigned i=0;i<=loops;++i) {
    MomentumLine line{std::vector<Rational>(loops,Rational(0)),{Rational(0)},masses_squared[i]};
    if(i<loops)line.loop_coefficients[i]=Rational(1);
    else {for(auto& c:line.loop_coefficients)c=Rational(-1);line.external_coefficients[0]=Rational(1);}
    f.lines.push_back(std::move(line));
  }
  return f;
}
}  // namespace diffexp::feynman

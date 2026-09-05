#pragma once
#include "diffexp/ibp.hpp"
#include <climits>

namespace diffexp::pullback {
struct Budget {
  std::size_t max_numerator_power=32,max_terms=10000,max_products=100000;
};
struct Plan {
  feynman::Operation operation=feynman::Operation::Direct;
  ibp::Relation source_row;
  ibp::Integral denominator_indices;
  unsigned left_power=0,right_power=0,merged_slot=0;
  Rational normalization{1};
};
inline bool equal(const ibp::Affine& a,const ibp::Affine& b) {
  return a.constant==b.constant&&a.linear==b.linear;
}
// For BetaIntegral, integrate normalization*x^left_power*(1-x)^right_power
// times source_row. For endpoint operations, take the indicated one-sided limit
// of the entire symbolic row of integrals; its coefficients are not evaluated
// separately. Direct is an exact identity at a regular value of x.
inline Plan plan(const ibp::PropagatorBasis& old_basis,const ibp::PropagatorBasis& new_basis,
    unsigned left,unsigned right,const ibp::Integral& target,const Exact& x,const Budget& budget={}) {
  const auto n=old_basis.denominators.size(),physical=old_basis.physical_count;
  if(!budget.max_terms||!budget.max_products)throw std::invalid_argument("pullback requires positive expansion budgets");
  if(left==right||left>=physical||right>=physical||target.size()!=n)
    throw std::invalid_argument("pullback merge positions or target shape");
  if(new_basis.physical_count+1!=physical||new_basis.denominators.size()!=n ||
     old_basis.space.loops!=new_basis.space.loops||old_basis.space.pairs!=new_basis.space.pairs||
     old_basis.space.external_gram!=new_basis.space.external_gram)
    throw std::invalid_argument("pullback scalar-product spaces or physical slots differ");
  (void)(x+old_basis.space.zero().constant);(void)(x+new_basis.space.zero().constant);
  for(std::size_t k=physical;k<n;++k)if(target[k]>0)throw std::invalid_argument("pullback forbids positive auxiliary indices");
  auto merged=ibp::plus(ibp::scaled(old_basis.denominators[left],x),
                       ibp::scaled(old_basis.denominators[right],x.constant(1)-x));
  auto merged_slot=left-(left>right);
  if(!equal(merged,new_basis.denominators[merged_slot]))throw std::invalid_argument("new basis does not contain the specified merged propagator");
  for(std::size_t k=0;k<physical;++k)if(k!=left&&k!=right)
    if(!equal(old_basis.denominators[k],new_basis.denominators[k-(k>right)]))
      throw std::invalid_argument("nonmerged physical propagator mapping changed");
  Plan out;out.merged_slot=merged_slot;out.denominator_indices.assign(n,0);
  const int a=std::max(0,target[left]),b=std::max(0,target[right]);
  if(a>INT_MAX-b)throw std::overflow_error("pullback merged denominator power");
  out.denominator_indices[merged_slot]=a+b;
  for(std::size_t k=0;k<physical;++k)if(k!=left&&k!=right&&target[k]>0)
    out.denominator_indices[k-(k>right)]=target[k];
  if(a&&b) {
    out.operation=feynman::Operation::BetaIntegral;out.left_power=a-1;out.right_power=b-1;
    // Restrict factorial work as well as polynomial work for hostile indices.
    if(std::size_t(a)+std::size_t(b)>budget.max_products)throw std::runtime_error("pullback beta-normalization budget exceeded");
    for(int k=1;k<=a+b-1;++k)out.normalization*=Rational(k);
    for(int k=1;k<a;++k)out.normalization=out.normalization/Rational(k);
    for(int k=1;k<b;++k)out.normalization=out.normalization/Rational(k);
  } else if(b)out.operation=feynman::Operation::LowerLimit;
  else if(a)out.operation=feynman::Operation::UpperLimit;
  std::size_t degree=0,products=0;
  for(int value:target)if(value<0) {
    auto power=std::size_t(-std::int64_t(value));
    if(power>budget.max_numerator_power-degree)throw std::runtime_error("pullback numerator degree budget exceeded");degree+=power;
  }
  // Carry powers as source-integral indices directly. Multiplication by D_j
  // decrements its source index, preserving cancellations with denominators.
  out.source_row.emplace(out.denominator_indices,x.constant(1));
  for(std::size_t k=0;k<n;++k)if(target[k]<0) {
    const auto rewritten=new_basis.rewrite(old_basis.denominators[k]);
    for(std::int64_t power=0;power< -std::int64_t(target[k]);++power) {
      ibp::Relation next;
      for(const auto& [indices,coefficient]:out.source_row) {
        auto insert=[&](const ibp::Integral& index,const Exact& value) {
          if(value.is_zero())return;
          if(++products>budget.max_products)throw std::runtime_error("pullback multiplication budget exceeded");
          ibp::add(next,index,value);
          if(next.size()>budget.max_terms)throw std::runtime_error("pullback expanded term budget exceeded");
        };
        insert(indices,coefficient*rewritten.constant);
        for(std::size_t j=0;j<n;++j)if(!rewritten.linear[j].is_zero()) {
          auto source=indices;if(source[j]==INT_MIN)throw std::overflow_error("pullback numerator index overflow");--source[j];
          insert(source,coefficient*rewritten.linear[j]);
        }
      }
      out.source_row=std::move(next);
      if(out.source_row.empty())break;
    }
  }
  return out;
}
} // namespace diffexp::pullback

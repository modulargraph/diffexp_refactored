#pragma once
#include "diffexp/feynman.hpp"
#include <functional>
#include <numeric>
#include <tuple>

namespace diffexp::ibp {

// A propagator is affine in independent loop scalar products. Its coefficients
// may themselves be rational functions of kinematics and Feynman parameters.
// The constant is separate: completing the scalar-product basis must not treat
// a mass term as another independent loop variable.
struct Affine {
  Exact constant;
  std::vector<Exact> linear;
};
inline Affine scaled(const Affine& a, const Exact& c) {
  Affine out{a.constant*c,{}};
  for(const auto& x:a.linear) out.linear.push_back(x*c);
  return out;
}
inline Affine plus(const Affine& a,const Affine& b) {
  if(a.linear.size()!=b.linear.size()) throw std::invalid_argument("affine shape mismatch");
  Affine out{a.constant+b.constant,{}};
  for(std::size_t i=0;i<a.linear.size();++i) out.linear.push_back(a.linear[i]+b.linear[i]);
  return out;
}

class ScalarProducts {
 public:
  unsigned loops;
  std::vector<std::vector<Exact>> external_gram;
  std::vector<std::pair<unsigned,unsigned>> pairs;
  ScalarProducts(unsigned l,std::vector<std::vector<Exact>> gram,const Exact& sample)
      : loops(l),external_gram(std::move(gram)),zero_(sample.constant(0)) {
    if(!loops) throw std::invalid_argument("scalar-product space needs loop momenta");
    for(const auto& row:external_gram)
      if(row.size()!=externals()) throw std::invalid_argument("external Gram shape");
    for(unsigned i=0;i<externals();++i) for(unsigned j=0;j<externals();++j)
      if(!(external_gram[i][j]==external_gram[j][i]))
        throw std::invalid_argument("external Gram must be symmetric");
    for(unsigned i=0;i<loops;++i)
      for(unsigned j=i;j<loops+externals();++j) pairs.emplace_back(i,j);
  }
  unsigned externals() const {return external_gram.size();}
  std::size_t size() const {return pairs.size();}
  Affine zero() const {return {zero_,std::vector<Exact>(size(),zero_)};}
  Affine dot(unsigned a,unsigned b) const {
    if(a>=loops+externals() || b>=loops+externals()) throw std::out_of_range("momentum index");
    if(a>b) std::swap(a,b);
    auto out=zero();
    if(a>=loops) out.constant=external_gram[a-loops][b-loops];
    else {
      auto it=std::find(pairs.begin(),pairs.end(),std::pair{a,b});
      out.linear.at(it-pairs.begin())=zero_.constant(1);
    }
    return out;
  }
  Affine contraction(const Affine& expression,unsigned derivative_loop,unsigned vector) const {
    if(expression.linear.size()!=size() || derivative_loop>=loops || vector>=loops+externals())
      throw std::invalid_argument("IBP contraction shape");
    auto result=zero();
    for(std::size_t n=0;n<size();++n) {
      auto [a,b]=pairs[n];
      if(a==derivative_loop) result=plus(result,scaled(dot(vector,b),expression.linear[n]));
      if(b==derivative_loop) result=plus(result,scaled(dot(vector,a),expression.linear[n]));
    }
    return result;
  }
 private:
  Exact zero_;
};

struct QuadraticFamily {
  ScalarProducts scalar_products;
  std::vector<Affine> physical;
  std::vector<Affine> auxiliary;
};
inline QuadraticFamily quadratic_family(const feynman::Family& family,const Exact& sample,
    std::size_t physical_count=0) {
  family.validate();
  if(!physical_count) physical_count=family.lines.size();
  if(physical_count>family.lines.size()) throw std::invalid_argument("physical denominator count exceeds family size");
  std::vector<std::vector<Exact>> gram;
  for(const auto& row:family.external_gram) {
    gram.emplace_back();
    for(const auto& value:row) gram.back().push_back(sample.constant(value));
  }
  ScalarProducts sp(family.loops,std::move(gram),sample);
  QuadraticFamily result{sp,{}, {}};
  for(const auto& line:family.lines) {
    auto d=sp.zero(); d.constant=sample.constant(line.mass_squared);
    auto c=line.loop_coefficients;
    c.insert(c.end(),line.external_coefficients.begin(),line.external_coefficients.end());
    // Convention D=m^2-q^2, matching the native Symanzik/gamma boundary.
    for(unsigned i=0;i<c.size();++i) for(unsigned j=0;j<c.size();++j)
      d=plus(d,scaled(sp.dot(i,j),sample.constant(-c[i]*c[j])));
    if(result.physical.size()<physical_count) result.physical.push_back(std::move(d));
    else result.auxiliary.push_back(std::move(d));
  }
  return result;
}
inline QuadraticFamily merge(const QuadraticFamily& f,unsigned left,unsigned right,const Exact& x) {
  if(left==right || left>=f.physical.size() || right>=f.physical.size())
    throw std::invalid_argument("quadratic merge positions");
  auto out=f;
  out.physical[left]=plus(scaled(f.physical[left],x),scaled(f.physical[right],x.constant(1)-x));
  out.physical.erase(out.physical.begin()+right);
  return out;
}

// Exact rank selection and inversion happen once per family, before numerical
// chart work. ISP slots complete the basis; they are never physical denominators.
class PropagatorBasis {
 public:
  explicit PropagatorBasis(QuadraticFamily family)
      : space(std::move(family.scalar_products)),physical_count(family.physical.size()),
        denominators(std::move(family.physical)),zero_(space.zero().constant) {
    if(!physical_count || physical_count>space.size())
      throw std::invalid_argument("physical propagators require partial fractions before IBP completion");
    std::map<std::size_t,std::vector<Exact>> rank;
    auto insert=[&](std::vector<Exact> row) {
      if(row.size()!=space.size()) throw std::invalid_argument("propagator scalar-product shape");
      for(const auto& [pivot,basis]:rank) {
        auto c=row[pivot];
        for(std::size_t j=pivot;j<row.size();++j) row[j]=row[j]-c*basis[j];
      }
      auto p=std::find_if(row.begin(),row.end(),[](const Exact& x){return !x.is_zero();});
      if(p==row.end()) return false;
      auto pivot=p-row.begin(); auto coefficient=*p;
      for(auto& x:row) x=x/coefficient;
      rank.emplace(pivot,std::move(row)); return true;
    };
    for(const auto& d:denominators)
      if(!insert(d.linear)) throw std::invalid_argument("dependent physical propagators need partial fractions before IBP completion");
    for(auto& d:family.auxiliary) {
      if(!insert(d.linear)) throw std::invalid_argument("explicit numerator slots are linearly dependent");
      denominators.push_back(std::move(d));
    }
    for(std::size_t i=0;i<space.size() && denominators.size()<space.size();++i) {
      auto d=space.zero(); d.linear[i]=zero_.constant(1);
      if(insert(d.linear)) denominators.push_back(std::move(d));
    }
    const auto n=space.size();
    std::vector<std::vector<Exact>> matrix, inverse(n,std::vector<Exact>(n,zero_));
    for(const auto& d:denominators) matrix.push_back(d.linear);
    for(std::size_t i=0;i<n;++i) inverse[i][i]=zero_.constant(1);
    for(std::size_t i=0;i<n;++i) {
      auto pivot=i;
      while(pivot<n && matrix[pivot][i].is_zero()) ++pivot;
      if(pivot==n) throw std::logic_error("completed propagator basis is singular");
      std::swap(matrix[i],matrix[pivot]);std::swap(inverse[i],inverse[pivot]);
      auto c=matrix[i][i];
      for(std::size_t j=0;j<n;++j) {matrix[i][j]=matrix[i][j]/c;inverse[i][j]=inverse[i][j]/c;}
      for(std::size_t k=0;k<n;++k) if(k!=i) {
        auto factor=matrix[k][i];
        if(factor.is_zero()) continue;
        for(std::size_t j=0;j<n;++j) {
          matrix[k][j]=matrix[k][j]-factor*matrix[i][j];
          inverse[k][j]=inverse[k][j]-factor*inverse[i][j];
        }
      }
    }
    for(std::size_t i=0;i<n;++i) {
      Affine s{zero_,inverse[i]};
      for(std::size_t j=0;j<n;++j) s.constant=s.constant-inverse[i][j]*denominators[j].constant;
      scalar_products_in_denominators.push_back(std::move(s));
    }
  }
  Affine rewrite(const Affine& value) const {
    if(value.linear.size()!=space.size()) throw std::invalid_argument("rewrite scalar-product shape");
    Affine out{value.constant,std::vector<Exact>(space.size(),zero_)};
    for(std::size_t i=0;i<space.size();++i)
      if(!value.linear[i].is_zero()) out=plus(out,scaled(scalar_products_in_denominators[i],value.linear[i]));
    return out;
  }
  ScalarProducts space;
  std::size_t physical_count;
  std::vector<Affine> denominators,scalar_products_in_denominators;
 private:
  Exact zero_;
};

using Integral=std::vector<int>;
struct SeedBudget {
  unsigned extra_denominator_powers=1;
  unsigned numerator_degree=1;
  std::size_t maximum_seeds=2000;
};
// Enumerate a finite union of sectors. Zero and negative physical powers are
// included; auxiliary slots are nonpositive. Never drop terms in the resulting
// identities merely because they lie outside this seed range.
template<class Consumer>
std::size_t for_each_seed(std::size_t physical_count,std::size_t total_count,
    const SeedBudget& budget,Consumer consume) {
  if(!physical_count || physical_count>total_count || total_count>256 ||
     budget.extra_denominator_powers>1024 || budget.numerator_degree>1024)
    throw std::invalid_argument("invalid IBP seed range");
  Integral seed(total_count,0);std::size_t count=0;
  auto visit=[&](auto&& self,std::size_t slot,unsigned dots,unsigned numerators)->void {
    if(slot==total_count) {
      if(count>=budget.maximum_seeds) throw std::runtime_error("IBP seed budget exhausted");
      ++count;consume(seed);return;
    }
    seed[slot]=0;self(self,slot+1,dots,numerators);
    for(unsigned n=1;n<=numerators;++n) {
      seed[slot]=-static_cast<int>(n);self(self,slot+1,dots,numerators-n);
    }
    if(slot<physical_count) for(unsigned n=0;n<=dots;++n) {
      seed[slot]=static_cast<int>(n+1);self(self,slot+1,dots-n,numerators);
    }
  };
  visit(visit,0,budget.extra_denominator_powers,budget.numerator_degree);
  return count;
}
struct IntegralOrder {
  bool operator()(const Integral& a,const Integral& b) const {
    auto grade=[](const Integral& x) {
      std::int64_t positive=0,dots=0,numerators=0;
      for(auto n:x) {if(n>0) {++positive;dots+=n-1;}else numerators-=std::int64_t(n);}
      return std::tuple{positive,dots,numerators};
    };
    auto ga=grade(a),gb=grade(b);
    return ga==gb?a<b:ga<gb;
  }
};
using Relation=std::map<Integral,Exact,IntegralOrder>;
template<class Key,class Compare> void add(std::map<Key,Exact,Compare>& to,const Key& key,const Exact& coefficient) {
  if(coefficient.is_zero()) return;
  auto it=to.find(key);
  if(it==to.end()) to.emplace(key,coefficient);
  else {it->second=it->second+coefficient;if(it->second.is_zero()) to.erase(it);}
}
template<class Key,class Compare>
void add_scaled(std::map<Key,Exact,Compare>& to,const std::map<Key,Exact,Compare>& from,const Exact& scale) {
  for(const auto& [key,value]:from) add(to,key,scale*value);
}
inline int shift_index(int value,int delta) {
  const auto result=std::int64_t(value)+delta;
  if(result<INT_MIN || result>INT_MAX) throw std::overflow_error("integral index shift");
  return static_cast<int>(result);
}
inline void insert_derivative(Relation& result,const Integral& indices,std::size_t denominator,
    const Affine& derivative,const Exact& factor) {
  auto raised=indices; raised[denominator]=shift_index(raised[denominator],1);
  add(result,raised,factor*derivative.constant);
  for(std::size_t j=0;j<indices.size();++j) if(!derivative.linear[j].is_zero()) {
    auto lowered=raised;lowered[j]=shift_index(lowered[j],-1);
    add(result,lowered,factor*derivative.linear[j]);
  }
}

class Generator {
 public:
  Generator(PropagatorBasis basis,Exact dimension)
      : basis_(std::move(basis)),dimension_(std::move(dimension)) {
    for(unsigned i=0;i<basis_.space.loops;++i)
      for(unsigned v=0;v<basis_.space.loops+basis_.space.externals();++v) {
        contractions_.emplace_back();
        for(const auto& d:basis_.denominators)
          contractions_.back().push_back(basis_.rewrite(basis_.space.contraction(d,i,v)));
      }
  }
  const PropagatorBasis& basis() const {return basis_;}
  std::vector<Relation> relations(const Integral& indices) const {
    validate(indices);std::vector<Relation> result;std::size_t contraction=0;
    for(unsigned i=0;i<basis_.space.loops;++i)
      for(unsigned v=0;v<basis_.space.loops+basis_.space.externals();++v,++contraction) {
        Relation r;
        if(i==v) add(r,indices,dimension_);
        for(std::size_t k=0;k<indices.size();++k) if(indices[k])
          insert_derivative(r,indices,k,contractions_[contraction][k],dimension_.constant(-std::int64_t(indices[k])));
        result.push_back(std::move(r));
      }
    return result;
  }
  Relation derivative(const Integral& indices,std::size_t parameter) const {
    validate(indices);Relation result;
    if(!dimension_.derivative(parameter).is_zero())
      throw std::invalid_argument("parameter differentiation cannot vary the integration dimension");
    for(const auto& row:basis_.space.external_gram) for(const auto& value:row)
      if(!value.derivative(parameter).is_zero())
        throw std::invalid_argument("varying external Gram requires a kinematic derivative operator");
    for(std::size_t k=0;k<indices.size();++k) if(indices[k]) {
      const auto& d=basis_.denominators[k];
      Affine derivative{d.constant.derivative(parameter),{}};
      for(const auto& c:d.linear) derivative.linear.push_back(c.derivative(parameter));
      insert_derivative(result,indices,k,basis_.rewrite(derivative),dimension_.constant(-std::int64_t(indices[k])));
    }
    return result;
  }
 private:
  void validate(const Integral& indices) const {
    if(indices.size()!=basis_.denominators.size()) throw std::invalid_argument("IBP integral dimension");
    for(std::size_t i=basis_.physical_count;i<indices.size();++i)
      if(indices[i]>0) throw std::invalid_argument("ISP slots cannot become physical denominators");
  }
  PropagatorBasis basis_;
  Exact dimension_;
  std::vector<std::vector<Affine>> contractions_;
};

// A finite exact reducer is useful for small families and as a verification
// provider for external reductions. No truncation of the generated equations is
// allowed. Uneliminated integrals are explicitly unresolved, not silently called
// masters. Each reduction returns its linear combination of input identities.
class ExactReducer {
 public:
  using Witness=std::map<std::size_t,Exact>;
  struct Reduction {Relation remainder;Witness witness;};
  explicit ExactReducer(const Exact& sample,std::size_t equation_limit=20000)
      : one_(sample.constant(1)),equation_limit_(equation_limit) {}
  void insert(Relation equation) {
    if(sources_.size()>=equation_limit_) throw std::runtime_error("exact IBP equation budget exhausted");
    const auto id=sources_.size();sources_.push_back(equation);
    auto reduced=reduce(std::move(equation));
    if(reduced.remainder.empty()) return;
    Witness witness{{id,one_}};
    add_scaled(witness,reduced.witness,-one_);
    auto pivot=std::prev(reduced.remainder.end());const auto key=pivot->first;const auto scale=one_/pivot->second;
    for(auto& [integral,c]:reduced.remainder)c=c*scale;
    for(auto& [source,c]:witness)c=c*scale;
    rows_.emplace(key,Reduction{std::move(reduced.remainder),std::move(witness)});
  }
  Reduction reduce(Relation input) const {
    Witness witness;
    // Every row only introduces smaller terms in the fixed integral ordering.
    for(auto it=rows_.rbegin();it!=rows_.rend();++it) {
      auto found=input.find(it->first);
      if(found==input.end()) continue;
      auto coefficient=found->second;
      add_scaled(input,it->second.remainder,-coefficient);
      add_scaled(witness,it->second.witness,coefficient);
    }
    return {std::move(input),std::move(witness)};
  }
  bool verify(const Relation& input,const Reduction& result) const {
    auto residual=input;add_scaled(residual,result.remainder,-one_);
    for(const auto& [id,coefficient]:result.witness) {
      if(id>=sources_.size()) return false;
      add_scaled(residual,sources_[id],-coefficient);
    }
    return residual.empty();
  }
  std::size_t equation_count() const {return sources_.size();}
  std::size_t rank() const {return rows_.size();}
 private:
  Exact one_;
  std::size_t equation_limit_;
  std::vector<Relation> sources_;
  std::map<Integral,Reduction,IntegralOrder> rows_;
};

struct DifferentialSystem {
  std::vector<Integral> ordered_basis;
  std::vector<std::vector<Exact>> matrix;
};

// Prepare one coordinate map for every derivative and observable consumer.
// This proves the span of the selected basis under the supplied identities;
// it does not establish the global master count.
class BasisReduction {
 public:
  BasisReduction(const ExactReducer& reducer,std::vector<Integral> masters,const Exact& sample)
      : reducer_(reducer),masters_(std::move(masters)),zero_(sample.constant(0)),one_(sample.constant(1)),
        equation_count_(reducer.equation_count()) {
    if(masters_.empty())throw std::invalid_argument("empty selected integral basis");
    for(std::size_t i=0;i<masters_.size();++i) {
      Relation input{{masters_[i],one_}};auto reduced=reducer_.reduce(input);
      if(!reducer_.verify(input,reduced))throw std::logic_error("IBP master reduction witness failed");
      auto residual=coordinates(std::move(reduced.remainder));
      if(residual.form.empty())throw std::invalid_argument("selected basis is zero or linearly dependent under supplied IBP identities");
      auto pivot=std::prev(residual.form.end());const auto key=pivot->first;const auto c=pivot->second;
      for(auto& [integral,value]:residual.form)value=value/c;
      for(auto& value:residual.coefficients)value=-value/c;
      residual.coefficients[i]=residual.coefficients[i]+one_/c;
      rows_.emplace(key,std::move(residual));
    }
  }
  const std::vector<Integral>& ordered_basis() const {return masters_;}
  std::vector<Exact> resolve(const Relation& input) const {
    if(reducer_.equation_count()!=equation_count_)
      throw std::logic_error("IBP provider changed after selected-basis compilation");
    auto reduced=reducer_.reduce(input);
    if(!reducer_.verify(input,reduced))throw std::logic_error("IBP combination reduction witness failed");
    auto residual=coordinates(std::move(reduced.remainder));
    if(!residual.form.empty())throw std::runtime_error("selected-basis closure is unresolved: more IBP identities or basis elements are required");
    auto identity=input;
    for(std::size_t j=0;j<masters_.size();++j)add(identity,masters_[j],-residual.coefficients[j]);
    const auto verified=reducer_.reduce(identity);
    if(!verified.remainder.empty() || !reducer_.verify(identity,verified))
      throw std::logic_error("selected-basis combination failed its exact identity check");
    return std::move(residual.coefficients);
  }
 private:
  struct CoordinateRow {Relation form;std::vector<Exact> coefficients;};
  const ExactReducer& reducer_;
  std::vector<Integral> masters_;
  Exact zero_,one_;
  std::size_t equation_count_;
  std::map<Integral,CoordinateRow,IntegralOrder> rows_;
  CoordinateRow coordinates(Relation form) const {
    std::vector<Exact> coefficients(masters_.size(),zero_);
    for(auto it=rows_.rbegin();it!=rows_.rend();++it) {
      auto found=form.find(it->first);if(found==form.end())continue;
      const auto c=found->second;add_scaled(form,it->second.form,-c);
      for(std::size_t j=0;j<coefficients.size();++j)
        coefficients[j]=coefficients[j]+c*it->second.coefficients[j];
    }
    return {std::move(form),std::move(coefficients)};
  }
};
inline DifferentialSystem differential_system(const Generator& generator,const ExactReducer& reducer,
    std::vector<Integral> masters,std::size_t parameter,const Exact& sample) {
  BasisReduction basis(reducer,std::move(masters),sample);
  DifferentialSystem result{basis.ordered_basis(),{}};
  for(const auto& master:result.ordered_basis)result.matrix.push_back(basis.resolve(generator.derivative(master,parameter)));
  return result;
}
} // namespace diffexp::ibp

#pragma once
// Sectors of a family whose integrals all vanish in dimensional regularization.
//
// With the Symanzik polynomials U and F of the lines of a sector and
// G = U + F, the sector is scaleless if constants k_i exist with
//   sum_i k_i a_i dG/da_i = G
// (R. N. Lee, A. A. Pomeransky, arXiv:1308.6676, section 4): rescaling the
// parameters a_i -> lambda^{k_i} a_i then rescales the parametric integrand
// by a power of lambda, so the integral equals a multiple of itself for every
// lambda and vanishes, with any powers of the lines and any numerator. The
// condition is one linear equation k.e = 1 per monomial a^e of G, so it is
// decided by the exponents of the monomials whose coefficients are not
// identically zero in the field of the family: no finite field, no sample
// point. U and F have different degrees and cannot cancel each other, and
// the test does not depend on the common sign of the quadratic form.
//
// A sector with a free loop momentum has U = 0 and is found as well.
#include "diffexp/ibp.hpp"
#include <map>
#include <vector>

namespace diffexp::scaleless {
namespace detail {
// Polynomials in the line parameters with exact coefficients.
using Exponent=std::vector<unsigned char>;
using Polynomial=std::map<Exponent,Exact>;
inline void add(Polynomial& out,const Exponent& e,const Exact& c) {
  if(c.is_zero())return;
  auto [it,inserted]=out.try_emplace(e,c);
  if(!inserted){it->second=it->second+c;if(it->second.is_zero())out.erase(it);}
}
inline Polynomial sum(Polynomial a,const Polynomial& b,bool negate=false) {
  for(const auto& [e,c]:b)add(a,e,negate?c.constant(0)-c:c);
  return a;
}
inline Polynomial product(const Polynomial& a,const Polynomial& b) {
  Polynomial out;
  for(const auto& [e,c]:a)for(const auto& [f,d]:b){auto g=e;for(std::size_t i=0;i<g.size();++i)g[i]+=f[i];add(out,g,c*d);}
  return out;
}
inline Polynomial determinant(const std::vector<std::vector<Polynomial>>& m,const Polynomial& one) {
  const auto n=m.size();
  if(!n)return one;
  if(n==1)return m[0][0];
  Polynomial out;
  for(std::size_t column=0;column<n;++column) {
    if(m[0][column].empty())continue;
    std::vector<std::vector<Polynomial>> minor;
    for(std::size_t i=1;i<n;++i){minor.emplace_back();for(std::size_t j=0;j<n;++j)if(j!=column)minor.back().push_back(m[i][j]);}
    out=sum(std::move(out),product(m[0][column],determinant(minor,one)),column%2);
  }
  return out;
}
} // namespace detail

// U and F of the physical lines, as polynomials in one parameter per line, up
// to a common sign of each.
struct Symanzik {detail::Polynomial U,F;};
inline Symanzik symanzik(const ibp::PropagatorBasis& b) {
  using namespace detail;
  const auto lines=b.physical_count;const auto L=b.space.loops,E=b.space.externals();
  const auto zero=b.space.zero().constant;const auto half=zero.constant(Rational(1)/Rational(2));
  const Exponent none(lines,0);const Polynomial one{{none,zero.constant(1)}};
  // sum_k a_k D_k = l.A.l + 2 l.B + C in the loop momenta l, B in the span of
  // the external momenta.
  std::vector<std::vector<Polynomial>> A(L,std::vector<Polynomial>(L)),B(L,std::vector<Polynomial>(E));
  Polynomial C;
  for(std::size_t k=0;k<lines;++k) {
    auto e=none;e[k]=1;
    add(C,e,b.denominators[k].constant);
    for(std::size_t q=0;q<b.space.pairs.size();++q) {
      const auto [i,j]=b.space.pairs[q];const auto& value=b.denominators[k].linear[q];
      if(value.is_zero())continue;
      if(j>=L)add(B[i][j-L],e,value*half);
      else if(i==j)add(A[i][i],e,value);
      else {add(A[i][j],e,value*half);add(A[j][i],e,value*half);}
    }
  }
  Symanzik out;out.U=determinant(A,one);
  out.F=product(out.U,C);
  for(unsigned i=0;i<L;++i)for(unsigned j=0;j<L;++j) {
    Polynomial contracted;
    for(unsigned e=0;e<E;++e)for(unsigned f=0;f<E;++f)if(!b.space.external_gram[e][f].is_zero()&&!B[i][e].empty()&&!B[j][f].empty()) {
      auto term=product(B[i][e],B[j][f]);for(auto& [exponent,c]:term)c=c*b.space.external_gram[e][f];
      contracted=sum(std::move(contracted),term);
    }
    if(contracted.empty())continue;
    std::vector<std::vector<Polynomial>> minor;
    for(unsigned r=0;r<L;++r)if(r!=j){minor.emplace_back();for(unsigned c=0;c<L;++c)if(c!=i)minor.back().push_back(A[r][c]);}
    out.F=sum(std::move(out.F),product(determinant(minor,one),contracted),(i+j)%2==0);
  }
  return out;
}

// One flag per sector, indexed by the bit mask of its lines.
inline std::vector<bool> sectors(const ibp::PropagatorBasis& b) {
  if(b.physical_count>=8*sizeof(std::size_t))throw std::invalid_argument("scaleless sector enumeration: more physical denominators than bits of a sector mask");
  const auto lines=b.physical_count;const auto polynomials=symanzik(b);
  std::vector<std::pair<std::size_t,const detail::Exponent*>> monomials;
  for(const auto* polynomial:{&polynomials.U,&polynomials.F})for(const auto& [e,c]:*polynomial) {
    std::size_t mask=0;for(std::size_t k=0;k<lines;++k)if(e[k])mask|=std::size_t(1)<<k;
    monomials.emplace_back(mask,&e);
  }
  std::vector<bool> out(std::size_t(1)<<lines);
  for(std::size_t mask=0;mask<out.size();++mask) {
    // Rows (e | 1) in echelon form; the system is inconsistent, and the
    // sector has a scale, once a row reduces to (0 | nonzero).
    std::map<std::size_t,std::vector<Rational>> echelon;bool consistent=true;
    for(const auto& [support,e]:monomials) {
      if(support&~mask)continue;
      std::vector<Rational> row;for(std::size_t k=0;k<lines;++k)row.emplace_back(long((*e)[k]));row.emplace_back(1);
      for(const auto& [pivot,reduced]:echelon)if(!row[pivot].is_zero()){const auto factor=row[pivot];for(std::size_t j=pivot;j<row.size();++j)row[j]=row[j]-factor*reduced[j];}
      std::size_t pivot=0;while(pivot<lines&&row[pivot].is_zero())++pivot;
      if(pivot==lines){if(!row[lines].is_zero()){consistent=false;break;}continue;}
      const auto factor=row[pivot];for(auto& value:row)value=value/factor;
      echelon.emplace(pivot,std::move(row));
    }
    out[mask]=consistent;
  }
  return out;
}
} // namespace diffexp::scaleless

#pragma once
#include "diffexp/families.hpp"
#include "diffexp/geometry.hpp"
#include "diffexp/ibp.hpp"
#include <flint/acb_mat.h>
#include <optional>

namespace diffexp::causal {
// A supplied contour orientation is physical input, not a proof of its
// homotopy. DerivedFromSymanzik: every real singular point of a level is
// oriented from that level's own second Symanzik polynomial (see Landau).
enum class Assurance { SuppliedPrescription, EuclideanPhysicalDomain, VerifiedHomotopy, DerivedFromSymanzik };
// The side of a level's real singular points: of every one of them for a
// supplied prescription, of the free ones (where the physical branch is
// analytic and both sides agree) for a derived one.
struct LevelPlan { int x_detour_sign=0; };
struct Prescription {
  // -1 is the Feynman propagator: D -> D - i0, hence F -> F - i0.
  int f_rim=-1;
  std::vector<LevelPlan> levels;
  std::string convention="native-D=m2-q2;log-F-rim;v1",provenance;
  Assurance assurance=Assurance::SuppliedPrescription;
  // Orient each real singular point of each level on its own, from the
  // level's Symanzik F; `levels` then only places the free points.
  bool derive_orientations=false;
  void validate(std::size_t depth)const {
    if((f_rim!=1 && f_rim!=-1) || levels.size()!=depth || provenance.empty() || convention.empty())
      throw std::invalid_argument("incomplete causal prescription: F rim, depth and provenance required");
    for(const auto& level:levels)if(level.x_detour_sign!=1 && level.x_detour_sign!=-1)
      throw std::invalid_argument("incomplete causal prescription: each x orientation must be +/-1");
  }
};
inline constexpr const char* derivation_revision="landau-orientation-v1";
// The automatic prescription: nothing is assumed about the kinematics.
inline Prescription derived(std::size_t depth,int f_rim=-1) {
  Prescription p;p.f_rim=f_rim;p.levels.assign(depth,LevelPlan{1});p.derive_orientations=true;
  p.assurance=Assurance::DerivedFromSymanzik;
  p.provenance="Each real singular point is passed on the side Im x = f_rim*sign(dF/dx) at its Landau pinch of the level's Symanzik F; points without a pinch are free.";
  p.validate(depth);return p;
}
inline std::vector<Rational> henn_anchors(){return {Rational("1/5"),Rational("3/10"),Rational("2/5"),Rational("1/2"),Rational("3/5"),Rational("7/10"),Rational("4/5")};}
// Local linearized F orientation only. The caller must establish a real simple
// root and sufficiently small detour; this does not prove a full contour.
inline int simple_root_orientation(int f_rim,const Rational& derivative,bool real_simple_root_proved) {
  if((f_rim!=1 && f_rim!=-1) || !real_simple_root_proved || derivative.is_zero())
    throw std::invalid_argument("simple-root orientation needs a proved real simple root and F rim");
  return f_rim*(derivative>Rational(0)?1:-1);
}
// The side of each zero of a scalar leaf's F(x) strictly inside (0,1), which
// must not vanish at the ends: Im x = f_rim*sign(F'). A repeated zero throws.
struct LeafZero {Rational low,high;int side;};
inline std::vector<LeafZero> leaf_orientations(const Exact& F,unsigned variable,int f_rim=-1,slong bits=256) {
  if(f_rim!=1 && f_rim!=-1)throw std::invalid_argument("leaf orientation needs an F rim");
  std::vector<LeafZero> out;
  for(const auto& zero:simple_real_roots(F.numerator(),variable,Rational(0),Rational(1),bits))
    out.push_back({zero.low,zero.high,f_rim*zero.slope});
  return out;
}
inline bool positive_on_open_orthant(const Exact& polynomial) {
  if(polynomial.is_zero())return false;
  const auto denominator=polynomial.denominator_terms();
  if(denominator.size()!=1)return false;
  for(auto power:denominator.front().powers)if(power)return false;
  for(const auto& term:polynomial.numerator_terms())if(term.coefficient<=Rational(0))return false;
  return denominator.front().coefficient>Rational(0);
}
inline bool euclidean(const feynman::ExampleFamily& definition) {
  auto physical=definition.momenta;physical.lines.resize(definition.physical_count);
  std::vector<std::string> names;for(unsigned i=0;i<definition.physical_count;++i)names.push_back("a"+std::to_string(i));
  ExactField field(names);Exact sample(field);std::vector<Exact> parameters;
  for(unsigned i=0;i<definition.physical_count;++i)parameters.push_back(sample.variable(i));
  auto geometry=feynman::symanzik(physical,parameters);
  return positive_on_open_orthant(geometry.U) && positive_on_open_orthant(geometry.F);
}
inline Prescription euclidean_family(const feynman::ExampleFamily& definition,std::size_t depth) {
  if(!euclidean(definition))
    throw std::invalid_argument("the configured physical U/F are not proved positive: the family is not Euclidean");
  Prescription p;p.f_rim=-1;p.levels.assign(depth,LevelPlan{1});p.assurance=Assurance::EuclideanPhysicalDomain;
  p.provenance="Exact configured physical Symanzik U/F have positive monomials; no open-simplex physical cut. Contour homotopy is not certified.";
  p.validate(depth);return p;
}
// standard_merge_order means the legacy left-associated (first two physical
// denominators) ladder in original ordered family. Caller must verify this.
inline Prescription current_example(const std::string& name,std::size_t depth,
    const std::vector<Rational>& anchors={},bool standard_merge_order=false) {
  Prescription p;
  if(name=="henn_double_pentagon_x0") {
    if(depth!=7 || anchors!=henn_anchors() || !standard_merge_order)
      throw std::invalid_argument("Henn prescription requires original seven-level merge order and prescribed anchors");
    p.f_rim=-1;p.levels.assign(7,LevelPlan{-1});p.levels.front().x_detour_sign=1;
    p.provenance="Examples/FeynmanTrick/HennDoublePentagonBoundary.wl: DeltaPrescriptionSign, LevelDeltaPrescriptionSigns, fixedParameterValues; original ordered left-associated ladder";
  } else {
    auto example=feynman::example_family(name);example.momenta.lines.resize(example.physical_count);
    std::vector<std::string> names;for(unsigned i=0;i<example.physical_count;++i)names.push_back("a"+std::to_string(i));
    ExactField field(names);Exact sample(field);std::vector<Exact> parameters;
    for(unsigned i=0;i<example.physical_count;++i)parameters.push_back(sample.variable(i));
    auto geometry=feynman::symanzik(example.momenta,parameters);
    if(!positive_on_open_orthant(geometry.U) || !positive_on_open_orthant(geometry.F))
      throw std::invalid_argument("example lacks proved positive physical U/F; explicit causal prescription required");
    p.f_rim=-1;p.levels.assign(depth,LevelPlan{1});p.assurance=Assurance::EuclideanPhysicalDomain;
    p.provenance="native families.hpp exact physical Symanzik U/F have strictly positive monomials: no open-positive-simplex physical cut; complex contour homotopy not certified";
  }
  p.validate(depth);return p;
}

// What the level's F says about one real point of its parameter.
//   Free:        the physical branch is analytic there; either side.
//   Oriented:    the contour must pass on `side` (the sign of Im x).
//   Undetermined: neither could be established; `reason` says why.
struct Orientation {
  enum class Kind {Free,Oriented,Undetermined} kind=Kind::Free;
  int side=0;
  std::string reason;
};

// The Landau analysis of one level: F(alpha,x) is the second Symanzik
// polynomial of the level's merged family, with one Feynman parameter alpha
// per propagator, as a polynomial in the level's parameter x. With F - i0
// (f_rim -1) the level's integrals are, for real x, boundary values from the
// side of x on which F acquires a negative imaginary part wherever it
// vanishes on the integration domain; a real point is singular for that
// branch only if the alpha contour is pinched there, which requires a
// solution of the Landau equations on a face of the simplex: alpha_i>0 on the
// face, zero off it, and dF/dalpha_i=0 on it. The statements used:
//
//  (sign)  every monomial coefficient of F in alpha is nonzero at the point
//          and they share one sign: F does not vanish in the open simplex
//          for x near the point. Free.
//  (slope) every monomial coefficient of dF/dx is, at the point, either an
//          identically vanishing polynomial or nonzero, and the nonzero ones
//          share one sign s: Im F has the sign of s*Im x on the whole real
//          simplex, so the real simplex itself is a valid contour for
//          Im x = f_rim*s. Oriented, for any number of loops and masses.
//  (face)  one loop only, where F=alpha^T M(x) alpha and a pinch on a face S
//          is a positive null vector of the principal submatrix M_SS.
//          det M_SS not identically zero and nonzero at the point: no pinch
//          on S. det M_SS identically zero and the rank at the point equal
//          to the generic rank: the null vectors are those present for every
//          x (the soft and collinear regions of massless lines), on which
//          dF/dx vanishes identically; they do not single out the point.
//          Otherwise S is a candidate: if the entries of dM_SS/dx share one
//          sign, any pinch on S has that sign of dF/dx; on an edge with a
//          massless line, where F = a_j (2 b a_k + c a_j), the zero of the
//          second factor enters through the massless vertex with the sign
//          of b'; else the null vector is read from the adjugate, and is
//          either of mixed sign (no pinch), or positive with a definite
//          alpha^T M' alpha, or undecided.
//
// Assigning a side to a point that turns out to have no pinch is harmless,
// so existence of a pinch is never claimed, only excluded. All decisions are
// made in ball arithmetic on the enclosure of the point; anything a ball
// cannot decide is reported as undetermined, never guessed. Levels with more
// than one loop use (sign) and (slope) only: their faces with vanishing U
// need the Symanzik polynomials of contracted subgraphs, which is not
// implemented.
class Landau {
 public:
  using B=kernel::ComplexBall;
  // The merged family of a level; its denominators depend on the parameter
  // named `parameter` of their exact field and on nothing else.
  explicit Landau(const ibp::PropagatorBasis& merged,const std::string& parameter="x")
      :count_(merged.physical_count),loops_(merged.space.loops),field_(names(merged.physical_count,parameter)),u_(field_),f_(field_) {
    const auto& old=merged.space.zero().constant.variables();
    const auto found=std::find(old.begin(),old.end(),parameter);
    if(found==old.end())throw std::invalid_argument("Landau analysis needs the level parameter");
    const std::size_t old_parameter=found-old.begin();
    const auto convert=[&](const Exact& value) {
      for(std::size_t i=0;i<old.size();++i)if(i!=old_parameter && !value.derivative(i).is_zero())
        throw std::invalid_argument("Landau analysis: level kinematics depend on "+old[i]);
      return Exact(field_,value.str());
    };
    const Exact zero(field_),one(field_,1),half=one/(one+one);
    const unsigned L=loops_,E=merged.space.externals();
    std::vector<std::vector<Exact>> Q(L,std::vector<Exact>(L,zero)),R(L,std::vector<Exact>(E,zero)),gram(E,std::vector<Exact>(E,zero));
    for(unsigned e=0;e<E;++e)for(unsigned f=0;f<E;++f)gram[e][f]=convert(merged.space.external_gram[e][f]);
    auto C=zero;
    for(std::size_t n=0;n<count_;++n) {
      const auto alpha=zero.variable(n+1);const auto& d=merged.denominators[n];
      C=C+alpha*convert(d.constant);
      for(std::size_t k=0;k<merged.space.pairs.size();++k) {
        if(d.linear[k].is_zero())continue;
        const auto [i,j]=merged.space.pairs[k];const auto c=alpha*convert(d.linear[k]);
        if(j<L){Q[i][j]=Q[i][j]-c*(i==j?one:half);if(i!=j)Q[j][i]=Q[i][j];}
        else R[i][j-L]=R[i][j-L]-c*half;
      }
    }
    // Symmetric Schur complements, as in feynman::symanzik.
    auto U=one;
    for(unsigned i=0;i<L;++i) {
      if(Q[i][i].is_zero()) {
        unsigned j=i+1;while(j<L && Q[j][j].is_zero())++j;
        if(j==L)throw std::domain_error("degenerate quadratic family: no symmetric pivot");
        std::swap(Q[i],Q[j]);for(auto& row:Q)std::swap(row[i],row[j]);std::swap(R[i],R[j]);
      }
      const auto pivot=Q[i][i];U=U*pivot;
      for(unsigned e=0;e<E;++e)for(unsigned f=0;f<E;++f)C=C+R[i][e]*R[i][f]*gram[e][f]/pivot;
      for(unsigned j=i+1;j<L;++j) {
        for(unsigned e=0;e<E;++e)R[j][e]=R[j][e]-Q[j][i]*R[i][e]/pivot;
        for(unsigned k=i+1;k<L;++k)Q[j][k]=Q[j][k]-Q[j][i]*Q[i][k]/pivot;
      }
    }
    u_=U;f_=U*C;
    if(!f_.denominator().is_rational() || !u_.denominator().is_rational())throw std::logic_error("Symanzik polynomials are not polynomial");
    const auto scale=f_.denominator().rational();
    for(const auto& term:f_.numerator_terms()) {
      std::vector<unsigned long> key(term.powers.begin()+1,term.powers.end());
      auto value=zero.constant(term.coefficient/scale)*zero.variable(0).pow(term.powers[0]);
      auto [at,inserted]=coefficients_.emplace(std::move(key),value);
      if(!inserted)at->second=at->second+value;
    }
    for(auto at=coefficients_.begin();at!=coefficients_.end();)at=at->second.is_zero()?coefficients_.erase(at):std::next(at);
    for(const auto& [key,value]:coefficients_)slopes_.emplace(key,value.derivative(0));
    if(loops_==1 && count_<=12) {
      matrix_.assign(count_,std::vector<Exact>(count_,zero));
      for(const auto& [key,value]:coefficients_) {
        std::vector<std::size_t> where;
        for(std::size_t i=0;i<key.size();++i)for(unsigned long k=0;k<key[i];++k)where.push_back(i);
        if(where.size()!=2)throw std::logic_error("one-loop F is not a quadratic form");
        if(where[0]==where[1])matrix_[where[0]][where[0]]=value;
        else matrix_[where[0]][where[1]]=matrix_[where[1]][where[0]]=value*half;
      }
      for(std::size_t mask=1;mask<(std::size_t(1)<<count_);++mask) {
        Face face(field_);for(std::size_t i=0;i<count_;++i)if(mask>>i&1)face.members.push_back(i);
        elimination(face);faces_.push_back(std::move(face));
      }
    }
  }
  unsigned loops()const{return loops_;}
  std::size_t propagators()const{return count_;}
  // In the field (x, landau_a0, landau_a1, ...).
  const Exact& U()const{return u_;}
  const Exact& F()const{return f_;}
  // A real point given by a rational enclosure. `f_rim` -1 is F - i0.
  Orientation orient(const Rational& low,const Rational& high,int f_rim,slong bits=0)const {
    if((f_rim!=1 && f_rim!=-1) || high<low)throw std::invalid_argument("Landau orientation needs an F rim and an ordered enclosure");
    if(!bits)bits=B::precision();
    B x;arb_set_fmpq(acb_realref(x.raw()),((low+high)/Rational(2)).raw(),bits);
    {B width;arb_set_fmpq(acb_realref(width.raw()),((high-low)/Rational(2)).raw(),bits);arb_add_error(acb_realref(x.raw()),acb_realref(width.raw()));}
    // (sign)
    {
      int common=0;bool decided=!coefficients_.empty();
      for(const auto& [key,value]:coefficients_) {
        const int s=sign(value,x,bits);
        if(s==undecided || (common && s!=common)){decided=false;break;}
        common=s;
      }
      if(decided)return {Orientation::Kind::Free,0,"F does not change sign on the simplex"};
    }
    // (slope)
    {
      int common=0;bool decided=true;
      for(const auto& [key,value]:slopes_) {
        if(value.is_zero())continue;
        const int s=sign(value,x,bits);
        if(s==undecided || (common && s!=common)){decided=false;break;}
        common=s;
      }
      if(decided && common)return {Orientation::Kind::Oriented,f_rim*common,"dF/dx has one sign on the simplex"};
    }
    if(matrix_.empty())
      return {Orientation::Kind::Undetermined,0,loops_==1?"the level has too many propagators for the face analysis":
        "F changes sign on the simplex and dF/dx has no common sign; the face analysis is implemented for one-loop levels only"};
    int votes[2]={0,0};
    for(const auto& face:faces_) {
      const auto verdict=candidate(face,x,bits);
      if(verdict.kind==Orientation::Kind::Undetermined) {
        std::string label;for(auto member:face.members)label+=(label.empty()?"":",")+std::to_string(member);
        return {Orientation::Kind::Undetermined,0,"face {"+label+"}: "+verdict.reason};
      }
      if(verdict.kind==Orientation::Kind::Oriented)++votes[verdict.side>0];
    }
    if(votes[0] && votes[1])return {Orientation::Kind::Undetermined,0,"pinches on different faces ask for opposite sides"};
    if(votes[0] || votes[1])return {Orientation::Kind::Oriented,f_rim*(votes[1]?1:-1),"Landau pinch"};
    return {Orientation::Kind::Free,0,"no Landau pinch on any face"};
  }
  // The one-loop pinch candidates strictly inside (0,1) and what was decided
  // there: a description of the level's physical cuts that does not depend
  // on the master basis, for cache identities.
  std::string signature(int f_rim,slong bits=0)const {
    if(!bits)bits=B::precision();
    std::string out=std::string(derivation_revision)+";loops="+std::to_string(loops_);
    std::map<std::string,std::string> points;
    for(const auto& face:faces_) {
      if(face.determinant.is_zero() || face.determinant.is_rational())continue;
      for(const auto& root:isolated_roots(face.determinant,0,bits)) {
        if(!arb_contains_zero(acb_imagref(root.value.raw())))continue;
        const auto low=rational_bound(acb_realref(root.value.raw()),false),high=rational_bound(acb_realref(root.value.raw()),true);
        if(!(Rational(0)<low) || !(high<Rational(1)))continue;
        const auto verdict=orient(low,high,f_rim,bits);
        points[root.value.real_midpoint(20)]=verdict.kind==Orientation::Kind::Free?"0":
          verdict.kind==Orientation::Kind::Oriented?std::to_string(verdict.side):"?";
      }
    }
    for(const auto& [at,side]:points)out+=";"+at+":"+side;
    return out;
  }
 private:
  struct Face {
    std::vector<std::size_t> members;
    Exact determinant;      // of M_SS, exact in x
    std::size_t rank=0;     // of M_SS over the rational functions of x
    explicit Face(const ExactField& field):determinant(field){}
  };
  static constexpr int undecided=2;
  static std::vector<std::string> names(std::size_t count,const std::string& parameter) {
    std::vector<std::string> out{parameter};for(std::size_t i=0;i<count;++i)out.push_back("landau_a"+std::to_string(i));return out;
  }
  static B value(const Exact& polynomial,const B& x,slong bits) {
    // Horner on the exact coefficients of a polynomial in x alone.
    std::map<unsigned long,Rational> terms;unsigned long degree=0;
    const auto scale=polynomial.denominator().rational();
    for(const auto& term:polynomial.numerator_terms()){terms.emplace(term.powers[0],term.coefficient/scale);degree=std::max(degree,term.powers[0]);}
    B out(0);
    for(unsigned long n=degree+1;n-->0;) {
      acb_mul(out.raw(),out.raw(),x.raw(),bits);
      if(auto at=terms.find(n);at!=terms.end()){B c;arb_set_fmpq(acb_realref(c.raw()),at->second.raw(),bits);acb_add(out.raw(),out.raw(),c.raw(),bits);}
    }
    return out;
  }
  // 0 for the zero polynomial, +/-1 for a certain sign at the point.
  static int sign(const Exact& polynomial,const B& x,slong bits) {
    if(polynomial.is_zero())return 0;
    const auto v=value(polynomial,x,bits);
    if(arb_is_positive(acb_realref(v.raw())))return 1;
    if(arb_is_negative(acb_realref(v.raw())))return -1;
    return undecided;
  }
  // Exact determinant and rank of a principal submatrix by elimination over
  // the rational functions of x.
  void elimination(Face& face)const {
    const auto k=face.members.size();const Exact zero(field_),one(field_,1);
    std::vector<std::vector<Exact>> m(k,std::vector<Exact>(k,zero));
    for(std::size_t i=0;i<k;++i)for(std::size_t j=0;j<k;++j)m[i][j]=matrix_[face.members[i]][face.members[j]];
    auto determinant=one;std::size_t rank=0;
    for(std::size_t column=0;column<k;++column) {
      std::size_t pivot=rank;while(pivot<k && m[pivot][column].is_zero())++pivot;
      if(pivot==k){determinant=zero;continue;}
      if(pivot!=rank){std::swap(m[pivot],m[rank]);determinant=-determinant;}
      determinant=determinant*m[rank][column];
      for(std::size_t i=rank+1;i<k;++i) {
        if(m[i][column].is_zero())continue;
        const auto factor=m[i][column]/m[rank][column];
        for(std::size_t j=column;j<k;++j)m[i][j]=m[i][j]-factor*m[rank][j];
      }
      ++rank;
    }
    face.determinant=rank==k?determinant:zero;face.rank=rank;
  }
  struct BallMatrix {
    acb_mat_t raw;
    BallMatrix(slong rows,slong columns){acb_mat_init(raw,rows,columns);}
    ~BallMatrix(){acb_mat_clear(raw);}
    BallMatrix(const BallMatrix&)=delete;
  };
  Orientation candidate(const Face& face,const B& x,slong bits)const {
    const auto k=static_cast<slong>(face.members.size());
    BallMatrix m(k,k);
    for(slong i=0;i<k;++i)for(slong j=0;j<k;++j)
      acb_set(acb_mat_entry(m.raw,i,j),value(matrix_[face.members[i]][face.members[j]],x,bits).raw());
    if(!face.determinant.is_zero()) {
      if(!value(face.determinant,x,bits).contains_zero())return {};
    } else {
      // A certified lower bound of the rank at the point: eliminate on
      // pivots that are certainly nonzero.
      BallMatrix work(k,k);acb_mat_set(work.raw,m.raw);
      std::vector<bool> used_row(k,false),used_column(k,false);std::size_t rank=0;
      for(;;) {
        slong pr=-1,pc=-1;
        for(slong i=0;i<k && pr<0;++i)if(!used_row[i])for(slong j=0;j<k;++j)
          if(!used_column[j] && !acb_contains_zero(acb_mat_entry(work.raw,i,j))){pr=i;pc=j;break;}
        if(pr<0)break;
        used_row[pr]=used_column[pc]=true;++rank;
        for(slong i=0;i<k;++i)if(!used_row[i]) {
          B factor;acb_div(factor.raw(),acb_mat_entry(work.raw,i,pc),acb_mat_entry(work.raw,pr,pc),bits);
          for(slong j=0;j<k;++j)if(!used_column[j])acb_submul(acb_mat_entry(work.raw,i,j),factor.raw(),acb_mat_entry(work.raw,pr,j),bits);
        }
      }
      if(rank>=face.rank)return {};
    }
    // The face is a candidate. Entries of dM_SS/dx of one sign settle the
    // sign of dF/dx at any pinch in the open face.
    {
      int common=0;bool decided=true;
      for(slong i=0;i<k && decided;++i)for(slong j=i;j<k;++j) {
        const auto derivative=matrix_[face.members[i]][face.members[j]].derivative(0);
        if(derivative.is_zero())continue;
        const int s=sign(derivative,x,bits);
        if(s==undecided || (common && s!=common)){decided=false;break;}
        common=s;
      }
      if(decided && common)return {Orientation::Kind::Oriented,common,""};
    }
    if(face.determinant.is_zero())
      return {Orientation::Kind::Undetermined,0,"the rank of a permanently degenerate face drops at the point"};
    // An edge with a massless line k: F = a_j (2 b a_k + c a_j) on it. The
    // factor a_j is an endpoint factor of the integrand; the zero of the
    // other factor enters the edge through the vertex k where b vanishes
    // (c nonzero there), and next to that vertex Im F has the sign of
    // b' Im x.
    if(k==2)for(slong a=0;a<2;++a) {
      const auto& massless=matrix_[face.members[a]][face.members[a]];
      const auto& other=matrix_[face.members[1-a]][face.members[1-a]];
      if(!massless.is_zero() || other.is_zero() || value(other,x,bits).contains_zero())continue;
      const int s=sign(matrix_[face.members[0]][face.members[1]].derivative(0),x,bits);
      if(s==1 || s==-1)return {Orientation::Kind::Oriented,s,""};
    }
    // A column of the adjugate with a certainly nonzero entry spans the null
    // space of a matrix of corank one.
    std::vector<std::vector<B>> adjugate(k,std::vector<B>(k));
    if(k==1)adjugate[0][0]=B(1);
    else for(slong i=0;i<k;++i)for(slong j=0;j<k;++j) {
      BallMatrix minor(k-1,k-1);
      for(slong r=0,rr=0;r<k;++r){if(r==j)continue;for(slong c=0,cc=0;c<k;++c){if(c==i)continue;acb_set(acb_mat_entry(minor.raw,rr,cc),acb_mat_entry(m.raw,r,c));++cc;}++rr;}
      acb_mat_det(adjugate[i][j].raw(),minor.raw,bits);
      if((i+j)%2)acb_neg(adjugate[i][j].raw(),adjugate[i][j].raw());
    }
    for(slong column=0;column<k;++column) {
      bool spans=false,positive=false,negative=false,doubtful=false;
      for(slong i=0;i<k;++i) {
        const auto* entry=acb_realref(adjugate[i][column].raw());
        if(arb_is_positive(entry)){positive=true;spans=true;}
        else if(arb_is_negative(entry)){negative=true;spans=true;}
        else doubtful=true;
      }
      if(!spans)continue;
      if(positive && negative)return {};
      if(doubtful)return {Orientation::Kind::Undetermined,0,"the null vector may lie on the boundary of the face"};
      // alpha^T M' alpha at the null vector.
      B slope(0);
      for(slong i=0;i<k;++i)for(slong j=0;j<k;++j) {
        const auto derivative=matrix_[face.members[i]][face.members[j]].derivative(0);
        if(derivative.is_zero())continue;
        slope+=adjugate[i][column]*adjugate[j][column]*value(derivative,x,bits);
      }
      if(arb_is_positive(acb_realref(slope.raw())))return {Orientation::Kind::Oriented,1,""};
      if(arb_is_negative(acb_realref(slope.raw())))return {Orientation::Kind::Oriented,-1,""};
      return {Orientation::Kind::Undetermined,0,"dF/dx is not separated from zero at the pinch"};
    }
    return {Orientation::Kind::Undetermined,0,"the corank at the point is not one"};
  }
  std::size_t count_;unsigned loops_;
  ExactField field_;
  Exact u_,f_;
  std::map<std::vector<unsigned long>,Exact> coefficients_,slopes_;
  std::vector<std::vector<Exact>> matrix_;
  std::vector<Face> faces_;
};
} // namespace diffexp::causal

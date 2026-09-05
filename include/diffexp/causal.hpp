#pragma once
#include "diffexp/families.hpp"

namespace diffexp::causal {
// A contour orientation is physical input, not a proof of its homotopy.
enum class Assurance { SuppliedPrescription, EuclideanPhysicalDomain, VerifiedHomotopy };
struct LevelPlan { int x_detour_sign=0; };
struct Prescription {
  int f_rim=0;
  std::vector<LevelPlan> levels;
  std::string convention="native-D=m2-q2;log-F-rim;v1",provenance;
  Assurance assurance=Assurance::SuppliedPrescription;
  void validate(std::size_t depth)const {
    if((f_rim!=1 && f_rim!=-1) || levels.size()!=depth || provenance.empty() || convention.empty())
      throw std::invalid_argument("incomplete causal prescription: F rim, depth and provenance required");
    for(const auto& level:levels)if(level.x_detour_sign!=1 && level.x_detour_sign!=-1)
      throw std::invalid_argument("incomplete causal prescription: each x orientation must be +/-1");
  }
};
inline std::vector<Rational> henn_anchors(){return {Rational("1/5"),Rational("3/10"),Rational("2/5"),Rational("1/2"),Rational("3/5"),Rational("7/10"),Rational("4/5")};}
// Local linearized F orientation only. The caller must establish a real simple
// root and sufficiently small detour; this does not prove a full contour.
inline int simple_root_orientation(int f_rim,const Rational& derivative,bool real_simple_root_proved) {
  if((f_rim!=1 && f_rim!=-1) || !real_simple_root_proved || derivative.is_zero())
    throw std::invalid_argument("simple-root orientation needs a proved real simple root and F rim");
  return f_rim*(derivative>Rational(0)?1:-1);
}
inline bool positive_on_open_orthant(const Exact& polynomial) {
  if(polynomial.is_zero())return false;
  const auto denominator=polynomial.denominator_terms();
  if(denominator.size()!=1)return false;
  for(auto power:denominator.front().powers)if(power)return false;
  for(const auto& term:polynomial.numerator_terms())if(term.coefficient<=Rational(0))return false;
  return denominator.front().coefficient>Rational(0);
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
} // namespace diffexp::causal

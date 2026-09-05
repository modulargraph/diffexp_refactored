#include "diffexp/canonical.hpp"
#include "diffexp/geometry.hpp"
#include "diffexp/rational_transport.hpp"
#include <iostream>
using namespace diffexp;
void require(bool b,const char* why){if(!b)throw std::runtime_error(why);}
int main(){try {
  using B=Jet::Ball;B::set_precision(256);
  const auto complex_error=finite_reference_error(B::from_strings("3","4"),B(0),256);
  require(complex_error>=5 && complex_error<5.000000000001,"reference error includes the complete complex norm");
  require(finite_reference_error(B::from_strings("[1 +/- 0.25]"),B(1),256)>=0.25,
      "reference comparison retains arithmetic radii");
  B indeterminate;acb_indeterminate(indeterminate.raw());
  for(bool invalid_reference:{false,true}) {
    bool failed=false;
    try {
      // A valid earlier component must not hide a later NaN in a maximum.
      double maximum=finite_reference_error(B(1),B(0),256);
      maximum=std::max(maximum,finite_reference_error(invalid_reference?B(0):indeterminate,
          invalid_reference?indeterminate:B(0),256));
      (void)maximum;
    } catch(const std::domain_error&) {failed=true;}
    require(failed,"nonfinite actual/reference coefficient cannot pass a maximum-error comparison");
  }
  B infinite;arb_pos_inf(acb_realref(infinite.raw()));
  bool invalid_bound=false;
  try{(void)finite_reference_error(infinite,B(0),256);}catch(const std::domain_error&){invalid_bound=true;}
  require(invalid_bound,"infinite coefficients cannot pass reference comparisons");
  Jet context(0,12,256);
  B::set_precision(64);
  auto precise=Jet(0,4,384).decimal("0.1234567890123456789012345678901234567890123456789012345678901234567890123456789").at(0);
  require(acb_rel_accuracy_bits(precise.raw())>375 && B::precision()==64,"decimal parsing honors the owned jet precision without changing ambient state");
  B::set_precision(256);
  auto number=evaluate(data::Reader("(* data only *) -2^(-1) + 3 I").read(),context,{}).at(0);
  require((number-B::from_strings("-1/2","3")).is_zero(),"reader operator precedence");
  auto scientific=evaluate(data::Reader("1.25`30*^-3").read(),context,{}).at(0);
  require((scientific-B::from_strings("1/800")).contains_zero(),"annotated scientific input");
  bool rejected=false;try{evaluate(data::Reader("Run[1]").read(),context,{});}catch(const std::invalid_argument&){rejected=true;}
  require(rejected,"data reader must reject executable forms");
  ExactField field({"t","r","I"});Exact t(field,"t"),r(field,"r");
  auto norm=polynomial_norm(t+r,1,t.constant(2));
  require(norm==Exact(field,"t^2-2"),"quadratic-field norm");
  auto roots=polynomial_roots(Exact(field,"(t^2+1)^2*(t-2)"),0);
  require(roots.size()==3,"squarefree complex root isolation");
  auto step=clearance_step(0,roots);
  require(step>0.24 && step<=0.25,"clearance includes complex singularities");
  Jet path(1,12,256);path.set(1,B::from_strings("0","1"));
  auto continued=continue_polynomial_sqrt(path.pow(4),B(0),B(2),B(1));
  require((continued-B::from_strings("-3","4")).contains_zero(),"root-sheet continuation across principal cut");
  // Canonical y=(1-t)^(-eps), with a known analytic epsilon^1 coefficient.
  Jet w(1,44,256);w.set(1,B(-1));
  auto values=canonical_chart({{0,0,0,Rational(-1)}},{w},{{B(1),B(0),B(0)}},B(0),B::from_strings("1/4"),40);
  B ref;acb_log(ref.raw(),B::from_strings("4/3").raw(),256);
  auto diff=values[0][1]-ref;arf_t upper;arf_init(upper);acb_get_abs_ubound_arf(upper,diff.raw(),256);
  require(arf_cmp_2exp_si(upper,-80)<0,"canonical chart logarithm reference");arf_clear(upper);
  // Incompatible singular starts must not be excused by a tiny midpoint.
  Jet singular(0,44,256);singular.set(1,B(1));rejected=false;
  try{canonical_chart({{0,0,0,Rational(1)}},{singular},{{B::from_strings("1e-60"),B(0)}},B(0),B(1),40);}
  catch(const std::runtime_error&){rejected=true;}
  require(rejected,"nonzero residue boundary rejected at every magnitude");
  for(const std::vector<std::string> names:{std::vector<std::string>{"eps","I","x"},std::vector<std::string>{"x"}}) {
    ExactField reordered(names);
    auto endpoint=rational_line({{0,0,0,Exact(reordered,"1/(x-2)")}},{{B::from_strings("1/2")}},8);
    auto error=endpoint[0][0]-B::from_strings("1/4");
    arf_t upper;arf_init(upper);acb_get_abs_ubound_arf(upper,error.raw(),256);
    require(arf_cmp_2exp_si(upper,-230)<0,"physical recurrence geometry uses variable identities, not field positions");arf_clear(upper);
  }
  std::cout<<"Ancillary parsing, algebraic norms, isolated-root clearance and canonical recurrence passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

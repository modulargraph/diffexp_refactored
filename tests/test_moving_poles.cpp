#include "diffexp/moving_poles.hpp"
#include <iostream>
using namespace diffexp;
namespace df=denominator_factors;
namespace mp=moving_poles;
using Matrix=ExactEpsilonMatrix;
void check(bool b,const char* message){if(!b)throw std::runtime_error(message);}
// The connection of J=S G for G'=A G.
Matrix gauged(const Matrix& a,const Matrix& s) {
  auto connection=fuchsify::detail::multiply(fuchsify::detail::multiply(s,a),fuchsify::detail::inverse(s));
  Matrix derivative=s;
  for(auto& row:derivative)for(auto& value:row)value=value.derivative(0);
  const auto extra=fuchsify::detail::multiply(derivative,fuchsify::detail::inverse(s));
  for(std::size_t i=0;i<a.size();++i)for(std::size_t j=0;j<a.size();++j)connection[i][j]=connection[i][j]+extra[i][j];
  return connection;
}
int main(){try {
  ExactField field({"x","eps","I"});
  auto e=[&](const char* s){return Exact(field,s);};
  // Classification and multiplicities of irreducible factors.
  {
    df::Catalogue catalogue(0,1);
    const auto& powers=catalogue.split(e("3*(x-1)*(eps+1)*(5*x-1+2*eps)^2"));
    check(powers.size()==3 && catalogue.factors().size()==3,"three irreducible factors");
    std::size_t parameter=0,epsilon=0,mixed=0;
    for(const auto& [k,power]:powers) {
      const auto kind=catalogue.factors()[k].kind;
      parameter+=kind==df::Kind::parameter;epsilon+=kind==df::Kind::epsilon;
      if(kind==df::Kind::mixed){++mixed;check(power==2,"mixed factor multiplicity");}
    }
    check(parameter==1 && epsilon==1 && mixed==1,"factor kinds");
    // A known factor is divided out, not factored again.
    check(catalogue.split(e("(5*x-1+2*eps)*(x+2)")).size()==2 && catalogue.factors().size()==4,"incremental catalogue");
    check(catalogue.order(e("(x-1)^2/(5*x-1+2*eps)^3"),powers[0].first)+catalogue.order(e("(x-1)^2/(5*x-1+2*eps)^3"),powers[1].first)
      +catalogue.order(e("(x-1)^2/(5*x-1+2*eps)^3"),powers[2].first)==-1,"order of a factor in a rational function");
    const auto weight=df::mixed_weight({{e("1/(5*x-1+2*eps)^2"),e("1/x")},{e("eps/(5*x-1+2*eps)"),e("1/eps")}},catalogue);
    check(weight.factors==1 && weight.total_order==2 && weight.entries==2,"mixed weight");
  }
  // A basis with no moving pole, and the same system in masters that become
  // dependent at 5x-1+2eps=0: J1=G1, J2=f G2+G1.
  const Matrix good{{e("eps/x"),e("0")},{e("1/(x-1)"),e("2*eps/(x-1)")}};
  const auto f=e("5*x-1+2*eps");
  const auto bad=gauged(good,{{e("1"),e("0")},{e("1"),f}});
  {
    df::Catalogue catalogue(0,1);
    check(df::mixed_weight(good,catalogue).total_order==0 && df::mixed_weight(bad,catalogue).total_order==1,"degenerate masters have a moving pole");
    auto result=mp::remove(bad,0,1);
    check(result.changed() && result.steps.size()==1 && result.after.total_order==0 && result.remaining.empty(),"closure row removes the moving pole");
    recursion::verify_physical_basis(result.transform,bad,0);
    check(df::mixed_weight(result.transform.from_basis,catalogue).total_order==0,"old masters stay polynomial in the factor");
    check(df::mixed_weight(result.transform.to_basis,catalogue).total_order==1,"new master carries the pole in its definition");
    for(const auto& row:result.transform.connection)for(const auto& value:row)if(!value.is_zero())
      for(const auto& [k,power]:catalogue.split(value.denominator())) {
        const auto& factor=catalogue.factors()[k].polynomial;
        check(catalogue.factors()[k].kind==df::Kind::epsilon || factor==e("x") || factor==e("x-1"),"no new singular point in x");
      }
    // The same step with a named integral: G2=(J2-J1)/f.
    mp::Options integrals_only;integrals_only.closure_rows=false;
    check(!mp::remove(bad,0,1,{},integrals_only).changed(),"nothing to try without rows");
    const std::vector<mp::Row> rows{{{e("1"),e("0")},"no pole"},{{e("-1")/f,e("1")/f},"integral G2"}};
    auto named=mp::remove(bad,0,1,rows,integrals_only);
    check(named.changed() && named.after.total_order==0 && named.steps[0].source=="integral G2" && named.steps[0].replaced==1,"integral with the pole in its reduction");
    recursion::verify_physical_basis(named.transform,bad,0);
    check(named.transform.connection==good,"the good basis is recovered");
    // Both kinds offered: the smaller connection wins.
    auto both=mp::remove(bad,0,1,rows);
    check(both.changed() && both.characters<=result.characters && both.characters<=named.characters,"smallest connection kept");
    // The specialised transform belongs to the specialised connection.
    auto shifted=mp::shift_epsilon(result.transform,Rational("1/7"));
    Matrix sampled=bad;
    for(auto& row:sampled)for(auto& value:row)if(!value.is_zero())value=value.substitute(std::vector<Exact>{e("x"),e("1/7"),e("I")});
    recursion::verify_physical_basis(shifted,sampled,0);
    for(const auto& row:shifted.connection)for(const auto& value:row)check(value.derivative(1).is_zero(),"specialised basis holds no epsilon");
    auto moved=mp::shift_epsilon(result.transform,Rational("1/7"),false);
    sampled=bad;
    for(auto& row:sampled)for(auto& value:row)if(!value.is_zero())value=value.substitute(std::vector<Exact>{e("x"),e("eps+1/7"),e("I")});
    recursion::verify_physical_basis(moved,sampled,0);
  }
  // A relation of second order, J2=f^2 G2+G1, takes two steps.
  {
    const auto worse=gauged(good,{{e("1"),e("0")},{e("1"),f*f}});
    auto result=mp::remove(worse,0,1);
    check(result.changed() && result.steps.size()==2 && result.after.total_order==0,"second-order relation");
    recursion::verify_physical_basis(result.transform,worse,0);
  }
  // Factors that are not linear in x: the relation holds in the residue field.
  {
    mp::ResidueField field(e("x^2-eps*x-3"),0);
    const auto a=field.residue(e("(x^3+eps)/(2*x-1)"));
    check(a.size()==2 && field.is_zero(field.reduce(field.coefficients((e("2*x-1")*field.value(a)-e("x^3+eps")).numerator()))),"residue of a quotient");
    const auto inverse=field.inverse(a);
    check(field.value(field.multiply(a,inverse))==e("1"),"inverse in the residue field");
    for(const char* factor:{"x^2-eps*x-3","2*x^3-x+eps^2+5"}) {
      const auto g=e(factor);
      const auto curved=gauged(good,{{e("1"),e("0")},{e("x+eps"),g}});
      df::Catalogue catalogue(0,1);
      check(df::mixed_weight(curved,catalogue).total_order==1,"nonlinear moving pole present");
      auto result=mp::remove(curved,0,1);
      check(result.changed() && result.steps.size()==1 && result.after.total_order==0 && result.remaining.empty(),"nonlinear moving pole removed");
      recursion::verify_physical_basis(result.transform,curved,0);
      check(df::mixed_weight(result.transform.from_basis,catalogue).total_order==0,"old masters polynomial in a nonlinear factor");
      for(const auto& row:result.transform.connection)for(const auto& value:row)if(!value.is_zero())
        for(const auto& [k,power]:catalogue.split(value.denominator())) {
          const auto& found=catalogue.factors()[k].polynomial;
          check(catalogue.factors()[k].kind==df::Kind::epsilon || found==e("x") || found==e("x-1"),"no new singular point from a nonlinear factor");
        }
    }
  }
  // A moving singularity that is not apparent, (x-eps)^eps, cannot be removed:
  // the steps never end, and the original basis is returned.
  {
    const Matrix genuine{{e("eps/(x-eps)")}};
    mp::Options bounded;bounded.max_factor_steps=4;
    auto result=mp::remove(genuine,0,1,{},bounded);
    check(!result.changed() && result.remaining.size()==1 && result.transform.connection==genuine &&
      result.transform.to_basis[0][0].is_one(),"genuine moving singularity left alone");
    // Next to an apparent one, only the apparent one goes.
    const auto g=e("3*x-2+eps");
    const auto pair=gauged({{e("eps/(x-eps)"),e("0")},{e("1/x"),e("2*eps/(x-1)")}},{{e("1"),e("0")},{e("1"),g}});
    auto partial=mp::remove(pair,0,1,{},bounded);
    check(partial.changed() && partial.remaining.size()==1 && partial.remaining[0]==e("x-eps").str() && partial.after.total_order==1,"apparent factor removed beside a genuine one");
    recursion::verify_physical_basis(partial.transform,pair,0);
  }
  // Singular points in x alone, on request: the apparent one at 3x-2 goes, the
  // genuine one at 2x-1 (exponent eps) stays, and the ends are never touched.
  {
    const auto h=e("3*x-2");
    const auto fixed=gauged({{e("eps/x+eps/(2*x-1)"),e("0")},{e("1/(x-1)"),e("2*eps/(x-1)")}},{{e("1"),e("0")},{e("1"),h}});
    check(!mp::remove(fixed,0,1).changed(),"fixed points are not touched by default");
    mp::Options apparent;apparent.fixed_points=true;apparent.max_factor_steps=4;
    auto result=mp::remove(fixed,0,1,{},apparent);
    check(result.changed() && result.steps.size()==1 && result.steps[0].factor==h.str() && result.remaining.size()==1 && result.remaining[0]==e("2*x-1").str(),"apparent fixed point");
    recursion::verify_physical_basis(result.transform,fixed,0);
  }
  // Reductions a level already holds become candidate rows; an identity with
  // two integrals outside the basis is not a reduction.
  {
    level::Result level;
    level.ordered_basis={{1,0},{1,1}};
    level.matrix={{e("0"),e("0")},{e("0"),e("0")}};
    level.requested_integrals={{1,1},{2,1}};
    level.target_rows={{e("0"),e("1")},{e("eps"),e("1/x")}};
    level.source_identities={
      ibp::Relation{{{1,2},e("1")},{{1,0},e("-eps")},{{1,1},e("-x")}},
      ibp::Relation{{{2,1},e("1")},{{1,0},e("-eps")},{{1,1},e("-1/x")}},
      ibp::Relation{{{2,2},e("1")},{{3,1},e("-1")},{{1,1},e("-x")}}};
    const auto rows=mp::cached_integrals(level);
    check(rows.size()==2 && rows[0].label=="integral [2,1]" && rows[1].label=="integral [1,2]","cached integral rows");
    check(rows[1].coordinates[0]==e("eps") && rows[1].coordinates[1]==e("x"),"identity read as a reduction");
  }
  std::cout<<"moving poles: all checks passed\n";
  return 0;
} catch(const std::exception& error){std::cerr<<"moving poles: "<<error.what()<<'\n';return 1;}}

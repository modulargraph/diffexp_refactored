#pragma once
#include "diffexp/ibp.hpp"
#include "diffexp/certified_deepest_beta.hpp"
#include "diffexp/frobenius.hpp"

namespace diffexp::feynman {
struct SunrisePlan {
  Exact parameter;
  ibp::DifferentialSystem system;
  std::vector<Exact> scalar_row;
  std::size_t generated_equations;
};
inline SunrisePlan prepare_sunrise() {
  ExactField field({"x","eps","I"});Exact x(field,"x");
  auto raw=ibp::quadratic_family(banana(2,{Rational(1),Rational(1),Rational(1)}),x);
  ibp::Generator generator(ibp::PropagatorBasis(ibp::merge(raw,0,1,x)),Exact(field,"2-2*eps"));
  ibp::ExactReducer reducer(x,300);
  ibp::for_each_seed(2,5,{1,1,50},[&](const ibp::Integral& seed) {
    for(auto& equation:generator.relations(seed))reducer.insert(std::move(equation));
  });
  ibp::BasisReduction basis(reducer,{{1,0,0,0,0},{1,1,0,0,0},{1,1,-1,0,0}},x);
  ibp::DifferentialSystem system{basis.ordered_basis(),{}};
  for(const auto& master:system.ordered_basis)system.matrix.push_back(basis.resolve(generator.derivative(master,0)));
  auto scalar_row=basis.resolve({{ibp::Integral{2,1,0,0,0},x.constant(1)}});
  return {x,std::move(system),std::move(scalar_row),reducer.equation_count()};
}
struct SunriseBoundary {
  int epsilon_low;
  Boundary values;
  bool taylor_tail_certified=false;
};
inline SunriseBoundary sunrise_boundary(const SunrisePlan& plan,const Rational& anchor,
    unsigned epsilon_top,unsigned taylor_order=80) {
  using B=Jet::Ball;
  if(anchor<=Rational(0) || anchor>=Rational(1) || anchor==Rational("1/2"))
    throw std::invalid_argument("sunrise boundary needs a nondegenerate interior anchor");
  const auto& x=plan.parameter;
  auto geometry=symanzik(banana(2,{Rational(1),Rational(1),Rational(1)}),
    {x*x.constant(anchor),x*x.constant(Rational(1)-anchor),x.constant(1)-x});
  CertifiedDeepestBetaOptions proof;proof.working_bits=B::precision();proof.requested_digits=28;
  proof.taylor_order=taylor_order;
  auto scalar11=certified_deepest_beta(geometry,2,2,1,1,epsilon_top,proof);
  auto scalar21=certified_deepest_beta(geometry,2,2,2,1,epsilon_top,proof);
  auto u=anchor*(Rational(1)-anchor);auto scalar10=tadpole(u,u,1,2,2,epsilon_top);
  const int low=std::min({scalar10.epsilon_low,scalar11.epsilon_low,scalar21.epsilon_low});
  const unsigned width=static_cast<int>(epsilon_top)-low+1;
  Jet epsilon(0,width,B::precision());epsilon.set(1,B(1));
  const auto pack=[&](int first,const std::vector<B>& values) {
    auto series=epsilon.constant(0);
    for(unsigned n=0;n<values.size();++n)if(static_cast<int>(n)+first>=low && static_cast<int>(n)+first<=static_cast<int>(epsilon_top))
      series.set(static_cast<int>(n)+first-low,values[n]);
    return series;
  };
  auto j10=pack(scalar10.epsilon_low,scalar10.coefficients),j11=pack(scalar11.epsilon_low,scalar11.coefficients),
       j21=pack(scalar21.epsilon_low,scalar21.coefficients);
  std::vector<Jet> row;
  auto at=epsilon.constant(0);at.set(0,B::from_strings(anchor.str()));
  for(const auto& coefficient:plan.scalar_row)
    row.push_back(evaluate(data::Reader(coefficient.str()).read(),epsilon,{{"x",at},{"eps",epsilon}}));
  auto numerator=(j21-row[0]*j10-row[1]*j11)/row[2];
  Boundary result(3,std::vector<B>(width,B(0)));
  for(unsigned k=0;k<width;++k){result[0][k]=j10.at(k);result[1][k]=j11.at(k);result[2][k]=numerator.at(k);}
  return {low,std::move(result),true};
}

inline std::vector<RationalLineEntry> sunrise_line(const SunrisePlan& plan,const Rational& from,const Rational& to) {
  const auto& x=plan.parameter;auto width=x.constant(to-from);
  std::vector<Exact> point{x.constant(from)+width*x,x.constant(0),x.variable(2)};
  std::vector<RationalLineEntry> entries;
  for(unsigned i=0;i<plan.system.matrix.size();++i)for(unsigned j=0;j<plan.system.matrix.size();++j) {
    const auto& a=plan.system.matrix[i][j];
    if(!a.derivative(1).derivative(1).is_zero())throw std::logic_error("sunrise connection unexpectedly nonaffine in epsilon");
    if(!a.is_zero())entries.push_back({i,j,0,width*a.substitute(point)});
    auto linear=a.derivative(1);
    if(!linear.is_zero())entries.push_back({i,j,1,width*linear.substitute(point)});
  }
  return entries;
}

// Equal-mass sunrise: integrate the outer beta functional using one shared
// three-master solution. Endpoint projection happens exactly before matching,
// so integrable cancellations are not converted into independent divergent
// component errors. The middle contour avoids an apparent pole of the IBP row.
inline LaurentValue sunrise(unsigned epsilon_top=0,unsigned endpoint_order=40,unsigned ordinary_order=80) {
  using B=Jet::Ball;
  if(epsilon_top>20 || endpoint_order<10 || endpoint_order>200)throw std::invalid_argument("sunrise series budget");
  auto plan=prepare_sunrise();const auto& x=plan.parameter;const Rational overlap("1/8");
  auto boundary=sunrise_boundary(plan,overlap,epsilon_top,ordinary_order);
  const auto width=static_cast<int>(epsilon_top)-boundary.epsilon_low+1;
  auto endpoint=FrobeniusSeries::prepare(plan.system.matrix,0,1,endpoint_order,width-1);
  auto constants=endpoint.match(B::from_strings(overlap.str()),boundary.values);
  auto projected=endpoint.project(plan.scalar_row,0,1);
  auto tail=endpoint.integrate_projected(projected,B::from_strings(overlap.str()),constants);
  boundary.values.push_back(std::vector<B>(width,B(0)));
  const auto left=x.constant(overlap),right=x.constant(Rational(1)-overlap);
  const auto imaginary=x.variable(2)*x.constant(Rational("1/10"));
  const std::vector<Exact> vertices{left,left+imaginary,right+imaginary,right};
  for(unsigned leg=0;leg+1<vertices.size();++leg) {
    const auto length=vertices[leg+1]-vertices[leg];
    std::vector<Exact> point{vertices[leg]+length*x,x.constant(0),x.variable(2)};
    std::vector<RationalLineEntry> entries;
    for(unsigned i=0;i<4;++i)for(unsigned j=0;j<3;++j) {
      const auto& coefficient=i==3?plan.scalar_row[j]:plan.system.matrix[i][j];
      if(!coefficient.derivative(1).derivative(1).is_zero())throw std::logic_error("sunrise augmented equation is nonaffine in epsilon");
      for(unsigned e=0;e<2;++e) {
        auto c=e?coefficient.derivative(1):coefficient;
        if(!c.is_zero())entries.push_back({i,j,e,length*c.substitute(point)});
      }
    }
    boundary.values=rational_line(entries,std::move(boundary.values),ordinary_order);
  }
  // Exact symmetry under interchange of the two equal-mass merged lines maps
  // x to 1-x, so both endpoint integrals equal the computed lower one.
  for(unsigned k=0;k<tail.size();++k)boundary.values[3][k]+=B(2)*tail[k];
  return {boundary.epsilon_low,std::move(boundary.values[3]),false};
}
} // namespace diffexp::feynman

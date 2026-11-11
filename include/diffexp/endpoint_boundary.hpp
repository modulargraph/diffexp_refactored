#pragma once
#include "diffexp/affine_frobenius.hpp"
#include "diffexp/epsilon_gauge.hpp"
#include "diffexp/families.hpp"
#include <algorithm>
#include <map>
#include <string>
#include <vector>

// Boundary data of a Feynman-trick level taken at an endpoint of its merge
// parameter instead of at an interior anchor.
//
// A level merges two lines into x*D_l+(1-x)*D_r. At x=0 the merged line is
// D_r: the level's integrals are integrals of the incoming family with line l
// removed (the pinched family), with its original kinematics. Where the
// level's integrals are analytic at that endpoint, the solution of the level's
// differential equation is fixed there by two exact statements: it has no
// component on a local solution with an epsilon-dependent exponent, and its
// analytic part is determined by the values of the masters at the endpoint,
// which are integrals of the pinched family.
//
// Analyticity is a statement about the integrals, not about the equation. It
// holds for lines l and r on one loop momentum in a family whose Symanzik
// polynomials have positive monomials, if the kept line is massive: U then
// depends on the two parameters through their sum only, F(x)-F(0) is bounded
// by a multiple of |x|*b*U (b the parameter of the merged line) and F(0) by
// m^2*b*U from below, so the parametric integrand is analytic in x uniformly
// over the simplex. With a massless kept line the endpoint has collinear
// regions with exponents -eps, -2*eps whose coefficients the pinched family
// does not provide: such an endpoint never qualifies.
namespace diffexp::endpoint_boundary {
using Matrix=fuchsify::Matrix;
using Expansion=AffineFrobeniusSeries::Expansion;

// The kinematic half of the rule. `composition` gives each present line as a
// combination of the original propagators; `kept` is the line that survives
// at the endpoint. `euclidean` is the positivity of the family's U and F.
struct Kinematics {bool qualifies=false;std::string reason;};
inline Kinematics kinematics(const feynman::ExampleFamily& family,bool euclidean,
    const std::vector<std::vector<Rational>>& composition,unsigned left,unsigned right,unsigned kept) {
  if(!euclidean)return {false,"the family's Symanzik polynomials are not proved positive"};
  if(left>=composition.size() || right>=composition.size() || (kept!=left && kept!=right))
    throw std::invalid_argument("endpoint boundary line positions");
  const auto& lines=family.momenta.lines;
  const std::vector<Rational>* routing=nullptr;
  for(const auto position:{left,right})for(std::size_t i=0;i<composition[position].size();++i) {
    if(composition[position][i].is_zero())continue;
    if(i>=lines.size())throw std::invalid_argument("endpoint boundary composition shape");
    const auto& loops=lines[i].loop_coefficients;
    if(!routing){routing=&loops;continue;}
    bool same=true,opposite=true;
    for(std::size_t k=0;k<loops.size();++k){same=same&&loops[k]==(*routing)[k];opposite=opposite&&loops[k]==-(*routing)[k];}
    if(!same && !opposite)return {false,"the merged lines carry different loop momenta"};
  }
  Rational mass(0);
  for(std::size_t i=0;i<composition[kept].size();++i) {
    if(composition[kept][i].sign()<0)return {false,"the kept line is not a positive combination of propagators"};
    if(!composition[kept][i].is_zero())mass=mass+composition[kept][i]*lines[i].mass_squared;
  }
  if(!(mass>Rational(0)))return {false,"the kept line is massless"};
  return {true,{}};
}

// Constants of the local solutions with a fixed exponent, from the values of
// the masters at the endpoint: constants[columns[k]] = sum_j solve[k][j] *
// value[j], every other constant being zero.
struct Determination {
  bool determined=false;
  std::string reason;
  std::vector<unsigned> columns;
  Matrix solve;
};
// `frame` is the fundamental solution in the masters' own coordinates, as an
// expansion at the endpoint (AffineFrobeniusSeries::project of the gauge that
// leads from the series' coordinates to the masters). The solution is
// determined if every fixed-exponent column is a Laurent series without
// logarithm, and if the absence of negative powers together with the values at
// the endpoint fixes all their constants.
inline Determination determine(const AffineFrobeniusSeries& series,const Expansion& frame,std::size_t ei) {
  using affine_frobenius_detail::is_integer;
  Determination out;
  const unsigned d=frame.columns;
  if(frame.rows!=d || d!=series.dimension())throw std::invalid_argument("endpoint boundary frame shape");
  std::vector<int> fixed(d,-1);
  for(const auto& term:frame.terms) {
    const int zero=term.slope.is_zero()?1:0;
    if(fixed[term.column]<0)fixed[term.column]=zero;
    else if(fixed[term.column]!=zero)throw std::logic_error("endpoint boundary column mixes fixed and regulated exponents");
    if(!zero)continue;
    if(!is_integer(term.power)){out.reason="a local solution has the fixed non-integer exponent "+term.power.str();return out;}
    if(term.log_degree){out.reason="a local solution with a fixed exponent has a logarithm";return out;}
  }
  for(unsigned c=0;c<d;++c) {
    if(fixed[c]<0){out.reason="a local solution has no retained term";return out;}
    if(fixed[c])out.columns.push_back(c);
  }
  if(out.columns.empty()){out.reason="no local solution is analytic at the endpoint";return out;}
  Matrix constant;std::vector<AffineFrobeniusSeries::ZeroConstraint> constraints;
  try {
    auto domain=series.dr_domain(frame,false);
    constant=series.dr_endpoint_constant(domain.admissible);
    constraints=std::move(domain.zero_constraints);
  } catch(const std::domain_error& error) {out.reason=error.what();return out;}
  const auto zero=constant[0][0].constant(0);
  const unsigned width=static_cast<unsigned>(out.columns.size());
  // [negative-power rows; value rows | 0; 1], reduced on the constants.
  Matrix system;
  for(const auto& constraint:constraints) {
    std::vector<Exact> row(width+d,zero);bool any=false;
    for(unsigned k=0;k<width;++k){row[k]=constraint.coefficients[out.columns[k]];any=any||!row[k].is_zero();}
    if(any)system.push_back(std::move(row));
  }
  for(unsigned i=0;i<d;++i) {
    std::vector<Exact> row(width+d,zero);
    for(unsigned k=0;k<width;++k)row[k]=constant[i][out.columns[k]];
    row[width+i]=zero.constant(1);system.push_back(std::move(row));
  }
  // Elimination with the pivot of lowest epsilon order: every multiplier is
  // then regular at epsilon=0, and the solution has no deeper pole than the
  // pivots ask for. A left inverse with spurious poles would ask the child
  // for epsilon orders that cancel in the result.
  for(unsigned column=0;column<width;++column) {
    std::optional<unsigned> pivot;std::int64_t lowest=0;
    for(unsigned row=column;row<system.size();++row) {
      const auto valuation=exact_epsilon_valuation(system[row][column],ei);
      if(valuation && (!pivot || *valuation<lowest)){pivot=row;lowest=*valuation;}
    }
    if(!pivot) {
      out.reason="the values at the endpoint do not fix the analytic solution (Taylor coefficients would be needed)";return out;
    }
    std::swap(system[*pivot],system[column]);
    const auto scale=system[column][column];
    for(auto& entry:system[column])if(!entry.is_zero())entry=entry/scale;
    for(unsigned row=0;row<system.size();++row) {
      if(row==column)continue;
      const auto factor=system[row][column];
      if(factor.is_zero())continue;
      for(unsigned j=0;j<system[row].size();++j)if(!system[column][j].is_zero())system[row][j]=system[row][j]-factor*system[column][j];
    }
  }
  out.solve=fuchsify::detail::zeros(width,d,zero);
  for(unsigned k=0;k<width;++k)for(unsigned j=0;j<d;++j)out.solve[k][j]=system[k][width+j];
  out.determined=true;
  return out;
}
// The masters at the distance `point` from the endpoint, as an exact matrix
// applied to their values at the endpoint: the retained Taylor polynomial of
// the determined solution. Its tail is not bounded here.
inline Matrix taylor_map(const Expansion& frame,const Determination& determination,const Rational& point) {
  if(!determination.determined)throw std::invalid_argument("endpoint boundary is not determined");
  const unsigned d=frame.rows,width=static_cast<unsigned>(determination.columns.size());
  const auto zero=determination.solve[0][0].constant(0);
  std::vector<int> position(frame.columns,-1);
  for(unsigned k=0;k<width;++k)position[determination.columns[k]]=static_cast<int>(k);
  auto values=fuchsify::detail::zeros(d,width,zero);
  std::map<long,Rational> powers;
  for(const auto& term:frame.terms) {
    if(position[term.column]<0)continue;
    const long power=affine_frobenius_detail::integer_value(term.power);
    auto found=powers.find(power);
    if(found==powers.end()) {
      Rational value(1);
      for(long k=0;k<std::labs(power);++k)value=power>0?value*point:value/point;
      found=powers.emplace(power,value).first;
    }
    auto& entry=values[term.row][position[term.column]];
    entry=entry+term.coefficient*zero.constant(found->second);
  }
  return fuchsify::detail::multiply(values,determination.solve);
}

// The equation half of the rule, for a level's connection dY/dx = A Y in the
// masters' coordinates: whether the analytic solutions at x=0 (x=1 for
// `upper`) are fixed by the masters' values there.
struct Options {
  fuchsify::Options fuchsification=[]{fuchsify::Options value;value.max_dimension=256;return value;}();
  AffineFrobeniusSeries::Options series=[]{AffineFrobeniusSeries::Options value;value.max_dimension=256;return value;}();
  // Series orders tried in turn: a basis whose gauge has a pole at the
  // endpoint needs more orders to reach the constant term.
  std::vector<unsigned> orders{4,8,16};
};
struct Verdict {
  bool qualifies=false;
  std::string reason;
  std::size_t analytic_dimension=0;
};
inline Verdict examine(const Matrix& connection,std::size_t xi,std::size_t ei,bool upper,const Options& options={}) {
  Verdict out;
  if(connection.empty()){out.reason="empty connection";return out;}
  try {
    auto matrix=connection;
    if(upper) {
      const auto& sample=connection[0][0];
      std::vector<Exact> point;for(std::size_t i=0;i<sample.variable_count();++i)point.push_back(sample.variable(i));
      point.at(xi)=sample.constant(1)-sample.variable(xi);
      for(auto& row:matrix)for(auto& entry:row)if(!entry.is_zero())entry=-entry.substitute(point);
    }
    const auto gauge=epsilon_diagonal_gauge(matrix,ei);
    const auto d=gauge.matrix.size();
    auto diagonal=fuchsify::detail::identity(static_cast<unsigned>(d),gauge.matrix[0][0]);
    for(std::size_t i=0;i<d;++i)diagonal[i][i]=epsilon_gauge_detail::multiply_power(diagonal[i][i],ei,gauge.shifts[i]);
    const auto fuchs=fuchsify::prepare(gauge.matrix,xi,options.fuchsification);
    if(!fuchs.success){out.reason="no Fuchsian form at the endpoint: "+fuchs.reason;return out;}
    const auto physical=fuchsify::detail::multiply(diagonal,fuchs.transform);
    for(const auto order:options.orders) {
      const auto series=AffineFrobeniusSeries::prepare(fuchs.matrix,xi,ei,order,options.series);
      const auto determination=determine(series,series.project(physical),ei);
      if(determination.determined){out.qualifies=true;out.analytic_dimension=determination.columns.size();out.reason.clear();return out;}
      out.reason=determination.reason;
      // Only a series too short to reach the constant term is worth extending.
      if(out.reason.find("increase the endpoint x order")==std::string::npos)return out;
    }
  } catch(const std::exception& error) {out.reason=error.what();}
  return out;
}
} // namespace diffexp::endpoint_boundary

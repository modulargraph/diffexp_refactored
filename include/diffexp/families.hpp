#pragma once
#include "diffexp/feynman.hpp"

namespace diffexp::feynman {

struct ExampleFamily {
  std::string name;
  Family momenta;
  unsigned physical_count;
  int dimension_at_epsilon_zero;
};
inline const std::vector<std::string>& example_names() {
  static const std::vector<std::string> names{
    "bubble","sunrise","banana","banana_unequal","banana4","banana4_unequal",
    "kite","box","pentagon","pentagon_massive","box_bubble","box_triangle",
    "double_box_planar","henn_double_pentagon_x0"};
  return names;
}
inline MomentumLine line(std::initializer_list<long> loops,std::initializer_list<long> externals,
    Rational mass_squared=Rational(0)) {
  MomentumLine result{{},{},std::move(mass_squared)};
  for(auto n:loops)result.loop_coefficients.emplace_back(n);
  for(auto n:externals)result.external_coefficients.emplace_back(n);
  return result;
}

// Exact kinematics and denominator order from Scripts/FTExamples.m and the
// original HennDoublePentagonBoundary.wl family. These are native definitions;
// loading an example never starts a CAS or reads an executable package.
inline ExampleFamily example_family(const std::string& name) {
  Family family{1,{}, {}};int dimension=2;unsigned physical=0;
  if(name=="bubble" || name=="sunrise" || name=="banana" || name=="banana4") {
    const unsigned loops=name=="bubble"?1:name=="sunrise"?2:name=="banana"?3:4;
    family=banana(loops,std::vector<Rational>(loops+1,Rational(1)));
  } else if(name=="banana_unequal") {
    family=banana(3,{Rational(2),Rational("3/2"),Rational("4/3"),Rational(1)});
  } else if(name=="banana4_unequal") {
    family=banana(4,{Rational(2),Rational("3/2"),Rational("4/3"),Rational("5/4"),Rational(1)});
  } else if(name=="kite") {
    family={2,{{Rational(-1)}},{
      line({1,0},{0},Rational(1)),line({1,0},{-1},Rational(1)),
      line({0,1},{0},Rational(1)),line({0,1},{-1},Rational(1)),
      line({1,-1},{0},Rational(1))}};
  } else if(name=="box" || name=="box_bubble" || name=="box_triangle" || name=="double_box_planar") {
    dimension=4;family.external_gram={
      {Rational(0),Rational("-1/2"),Rational("2/3")},
      {Rational("-1/2"),Rational(0),Rational("-1/6")},
      {Rational("2/3"),Rational("-1/6"),Rational(0)}};
    if(name=="box") family.lines={line({1},{0,0,0}),line({1},{1,0,0}),line({1},{1,1,0}),line({1},{1,1,1})};
    else {
      family.loops=2;family.lines={
        line({1,0},{0,0,0}),line({1,0},{1,0,0}),line({1,0},{1,1,0}),line({1,-1},{0,0,0}),
        line({0,1},{0,0,0}),line({0,1},{1,1,0}),line({0,1},{1,1,1})};
      if(name=="box_bubble")family.lines.erase(family.lines.begin(),family.lines.begin()+2);
      else if(name=="box_triangle")family.lines.erase(family.lines.begin());
    }
  } else if(name=="pentagon" || name=="pentagon_massive") {
    dimension=4;
    family.lines={line({1},{0,0,0,0}),line({1},{1,0,0,0}),line({1},{1,1,0,0}),
      line({1},{1,1,1,0}),line({1},{1,1,1,1})};
    if(name=="pentagon") family.external_gram={
      {Rational(0),Rational("-1/2"),Rational(-1),Rational(5)},
      {Rational("-1/2"),Rational(0),Rational(-1),Rational(-1)},
      {Rational(-1),Rational(-1),Rational(0),Rational("-3/2")},
      {Rational(5),Rational(-1),Rational("-3/2"),Rational(0)}};
    else {
      family.external_gram.assign(4,std::vector<Rational>(4,Rational("1/4")));
      for(unsigned i=0;i<4;++i)family.external_gram[i][i]=Rational(-1);
      for(unsigned i=0;i<5;++i)family.lines[i].mass_squared=Rational(i+2)/Rational(i+1);
      family.lines[0].mass_squared=Rational(1);
    }
  } else if(name=="henn_double_pentagon_x0") {
    dimension=4;physical=8;family.loops=2;
    family.external_gram={
      {Rational(0),Rational("3/2"),Rational("-1/2"),Rational("-1/2")},
      {Rational("3/2"),Rational(0),Rational("-1/2"),Rational("-1/2")},
      {Rational("-1/2"),Rational("-1/2"),Rational(0),Rational("1/2")},
      {Rational("-1/2"),Rational("-1/2"),Rational("1/2"),Rational(0)}};
    family.lines={
      line({1,0},{0,0,0,0}),line({1,0},{-1,0,0,0}),line({1,0},{-1,-1,0,0}),
      line({0,1},{0,0,0,0}),line({0,1},{-1,-1,-1,0}),line({0,1},{-1,-1,-1,-1}),
      line({1,-1},{0,0,0,0}),line({1,-1},{0,0,1,0}),
      line({1,0},{-1,-1,-1,-1}),line({0,1},{-1,0,0,0}),line({0,1},{-1,-1,0,0})};
  } else throw std::invalid_argument("unknown native Feynman family: "+name);
  family.validate();if(!physical)physical=family.lines.size();
  return {name,std::move(family),physical,dimension};
}
} // namespace diffexp::feynman

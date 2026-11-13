#pragma once
// Physical points reached from a Euclidean one along a straight line in the
// kinematics (docs/feynman-families.md, "Continuation from a Euclidean point").
//
// With P the requested point and V any point of the same family whose U and F
// have positive monomials, the start E = P + t V (t > 0) has
//
//   F_P - F_E = -t F_V,
//
// every monomial coefficient of F decreases along E -> P. On the line
// lambda -> E + lambda (P - E) the second Symanzik polynomial is linear in
// lambda, so for Im lambda > 0
//
//   Im F(alpha, lambda) = Im(lambda) (F_P - F_E)(alpha) < 0
//
// on the whole open simplex: the integrals are analytic in the upper half
// lambda plane and their boundary value at lambda = 1 from above is F - i0,
// the Feynman prescription. No real singular point of the line needs a side
// of its own; the path only has to stay in the upper half plane.
#include "diffexp/causal.hpp"
#include "diffexp/epsilon_gauge.hpp"
#include "diffexp/laurent_transport.hpp"
#include "diffexp/level_preparation.hpp"
#include "diffexp/geometry.hpp"
#include <set>
#include <boost/json.hpp>
#include <map>
#include <regex>
#include <sys/wait.h>
#include <unistd.h>

namespace diffexp::continuation {

struct Plan {
  feynman::ExampleFamily start;
  Rational scale;
  Rational minimum_scale;
};

inline feynman::ExampleFamily with_kinematics(feynman::ExampleFamily family,
    const std::vector<std::vector<Rational>>& gram,const std::vector<Rational>& masses) {
  if(gram.size()!=family.momenta.external_gram.size())throw std::invalid_argument("continuation external_gram shape");
  for(std::size_t i=0;i<gram.size();++i)if(gram[i].size()!=gram.size())throw std::invalid_argument("continuation external_gram shape");
  if(masses.size()!=family.momenta.lines.size())throw std::invalid_argument("continuation mass_squared needs one entry per propagator");
  family.momenta.external_gram=gram;
  for(std::size_t i=0;i<masses.size();++i)family.momenta.lines[i].mass_squared=masses[i];
  return family;
}

// Monomial coefficients of the second Symanzik polynomial of the physical lines.
inline std::map<std::vector<unsigned long>,Rational> f_coefficients(const feynman::ExampleFamily& definition) {
  auto physical=definition.momenta;physical.lines.resize(definition.physical_count);
  std::vector<std::string> names;for(unsigned i=0;i<definition.physical_count;++i)names.push_back("a"+std::to_string(i));
  ExactField field(names);Exact sample(field);std::vector<Exact> parameters;
  for(unsigned i=0;i<definition.physical_count;++i)parameters.push_back(sample.variable(i));
  const auto F=feynman::symanzik(physical,parameters).F;
  const auto denominator=F.denominator_terms();
  if(denominator.size()!=1)throw std::logic_error("continuation: F is not a polynomial");
  for(auto power:denominator.front().powers)if(power)throw std::logic_error("continuation: F is not a polynomial");
  std::map<std::vector<unsigned long>,Rational> out;
  for(const auto& term:F.numerator_terms())out[term.powers]=term.coefficient/denominator.front().coefficient;
  return out;
}

// A simple rational at least `minimum`: 1, 2 or 5 times a power of ten.
inline Rational simple_above(const Rational& minimum) {
  Rational value(1);
  while(value/Rational(10)>=minimum && value>Rational("1/1000000"))value=value/Rational(10);
  for(;;) {
    for(long f:{1,2,5})if(Rational(f)*value>=minimum)return Rational(f)*value;
    value=value*Rational(10);
  }
}

inline Plan plan(const feynman::ExampleFamily& target,const std::vector<std::vector<Rational>>& gram,
    const std::vector<Rational>& masses,std::optional<Rational> scale) {
  const auto direction=with_kinematics(target,gram,masses);
  // The line keeps every relation of P that V shares; vanishing invariants and
  // masses, and equal masses and squared momenta, define the master basis.
  std::vector<Rational> p,v;
  for(std::size_t i=0;i<gram.size();++i){p.push_back(target.momenta.external_gram[i][i]);v.push_back(gram[i][i]);}
  for(unsigned i=0;i<target.physical_count;++i){p.push_back(target.momenta.lines[i].mass_squared);v.push_back(masses[i]);}
  for(std::size_t i=0;i<p.size();++i) {
    if(p[i].sign()==0 && v[i].sign()!=0)throw std::invalid_argument("continuation direction must keep the vanishing squared momenta and masses of the family");
    for(std::size_t j=i+1;j<p.size();++j)if(p[i]==p[j] && !(v[i]==v[j]))
      throw std::invalid_argument("continuation direction must keep the equal squared momenta and masses of the family");
  }
  for(std::size_t i=0;i<gram.size();++i)for(std::size_t j=0;j<gram.size();++j)
    if(target.momenta.external_gram[i][j].sign()==0 && gram[i][j].sign()!=0)
      throw std::invalid_argument("continuation direction must keep the vanishing scalar products of the family");
  for(std::size_t i=target.physical_count;i<masses.size();++i)if(!(masses[i]==target.momenta.lines[i].mass_squared) && masses[i].sign()!=0)
    throw std::invalid_argument("continuation direction: numerator slots carry no mass along the line; give them 0");
  if(!causal::euclidean(direction))
    throw std::invalid_argument("continuation euclidean_direction is not Euclidean: its U and F must have positive monomials");
  const auto fp=f_coefficients(target),fv=f_coefficients(direction);
  Rational minimum(0);
  for(const auto& [monomial,c]:fp) {
    if(c.sign()>0)continue;
    auto found=fv.find(monomial);
    if(found==fv.end())throw std::invalid_argument("continuation direction lacks an F monomial that is not positive at the target point");
    auto needed=Rational(0)-c/found->second;if(needed>minimum)minimum=needed;
  }
  Plan out{target,Rational(0),minimum};
  out.scale=scale?*scale:simple_above(minimum*Rational(2)+Rational("1/1000"));
  if(!(out.scale>minimum))throw std::invalid_argument("continuation scale "+out.scale.str()+" must exceed "+minimum.str()+" for a Euclidean start");
  auto g=target.momenta.external_gram;
  for(std::size_t i=0;i<g.size();++i)for(std::size_t j=0;j<g.size();++j)g[i][j]=g[i][j]+out.scale*gram[i][j];
  std::vector<Rational> m;
  for(std::size_t i=0;i<masses.size();++i)m.push_back(target.momenta.lines[i].mass_squared+(i<target.physical_count?out.scale*masses[i]:Rational(0)));
  out.start=with_kinematics(target,g,m);
  out.start.name=target.name+"__euclidean_start";
  if(!causal::euclidean(out.start))throw std::logic_error("continuation start point is not Euclidean");
  return out;
}

// Runs this program on another family file and returns its standard output;
// its progress goes to the same standard error.
inline std::string run_capture(const std::vector<std::string>& arguments) {
  int channel[2];if(pipe(channel)!=0)throw std::runtime_error("continuation: pipe failed");
  const pid_t child=fork();
  if(child<0)throw std::runtime_error("continuation: fork failed");
  if(child==0) {
    close(channel[0]);dup2(channel[1],STDOUT_FILENO);close(channel[1]);
    std::vector<char*> argv;for(const auto& a:arguments)argv.push_back(const_cast<char*>(a.c_str()));argv.push_back(nullptr);
    execvp(argv[0],argv.data());_exit(127);
  }
  close(channel[1]);std::string output;char block[65536];
  for(;;){auto got=read(channel[0],block,sizeof block);if(got<=0)break;output.append(block,static_cast<std::size_t>(got));}
  close(channel[0]);int status=0;waitpid(child,&status,0);
  if(output.empty())throw std::runtime_error("continuation: the boundary run at the Euclidean start produced no response (exit status "+std::to_string(WIFEXITED(status)?WEXITSTATUS(status):-1)+")");
  return output;
}

// One Arb decimal interval string, "[mid +/- rad]", "[+/- rad]" or "mid".
inline void read_arb(arb_t out,const std::string& text,slong bits) {
  static const std::regex ball(R"(^\s*\[\s*([^\s\]+/]*?)\s*(?:\+/-\s*([^\]\s]+))?\s*\]\s*$)");
  std::smatch m;std::string mid=text,rad;
  if(std::regex_match(text,m,ball)){mid=m[1].str();rad=m[2].str();}
  if(mid.empty())mid="0";
  if(arb_set_str(out,mid.c_str(),bits)!=0)throw std::invalid_argument("continuation: cannot read the number "+text);
  if(!rad.empty()) {
    arb_t r;arb_init(r);
    if(arb_set_str(r,rad.c_str(),bits)!=0){arb_clear(r);throw std::invalid_argument("continuation: cannot read the radius in "+text);}
    mag_t e;mag_init(e);arb_get_mag(e,r);arb_add_error_mag(out,e);mag_clear(e);arb_clear(r);
  }
}
inline LaurentBoundary read_rows(const boost::json::object& response,slong bits) {
  LaurentBoundary out{static_cast<int>(response.at("epsilon_low").as_int64()),{}};
  for(const auto& row:response.at("coefficients").as_array()) {
    out.values.emplace_back();
    for(const auto& entry:row.as_array()) {
      kernel::ComplexBall value;const auto& o=entry.as_object();
      read_arb(acb_realref(value.raw()),std::string(o.at("real").as_string()),bits);
      read_arb(acb_imagref(value.raw()),std::string(o.at("imaginary").as_string()),bits);
      out.values.back().push_back(std::move(value));
    }
  }
  out.high();
  return out;
}

// Transports the masters Y (Y_i = eps^shift_i Z_i, the gauged Z have a
// holomorphic connection) from the start along 0 -> i h -> 1 + i h -> 1 and
// returns the target rows at the endpoint through epsilon^top.
inline LaurentBoundary transport_to_target(const ExactEpsilonMatrix& gauged,const std::vector<std::int64_t>& shifts,
    const LaurentBoundary& start,const std::vector<std::vector<Exact>>& targets,int top,const Rational& height,
    unsigned taylor_order,slong bits,std::ostream* progress=nullptr) {
  using B=kernel::ComplexBall;
  const auto n=gauged.size();
  LaurentBoundary result{0,{}};
  if(!n) {
    result.values.assign(targets.size(),std::vector<B>(static_cast<std::size_t>(top+1),B(0)));
    return result;
  }
  const auto smax=*std::max_element(shifts.begin(),shifts.end()),smin=*std::min_element(shifts.begin(),shifts.end());
  const int lowZ=start.low-static_cast<int>(smax),highZ=start.high()-static_cast<int>(smax);
  LaurentBoundary z{lowZ,{}};
  for(std::size_t i=0;i<n;++i) {
    z.values.emplace_back();
    for(int k=lowZ;k<=highZ;++k) {
      const int source=k+static_cast<int>(shifts[i]);
      z.values.back().push_back(source>=start.low?start.values[i][static_cast<std::size_t>(source-start.low)]:B(0));
    }
  }
  const auto& sample=gauged[0][0];const auto h=height.str();
  std::vector<Exact> vertices{sample.constant(0),sample.parse("("+h+")*I"),sample.parse("1+("+h+")*I"),sample.constant(1)};
  // The singular points of the line next to its two ends. The path leaves
  // and reaches the real axis there; a point at an end is refused.
  {
    std::set<std::string> seen;std::vector<B> roots;
    std::vector<Exact> at_zero_epsilon{sample.variable(0),sample.constant(0),sample.variable(2)};
    for(const auto& row:gauged)for(const auto& entry:row) {
      if(entry.is_zero())continue;
      auto denominator=entry.denominator().substitute(at_zero_epsilon);
      if(denominator.is_rational() || !seen.insert(denominator.str()).second)continue;
      for(auto& root:polynomial_roots(denominator,0,bits))roots.push_back(std::move(root));
    }
    double nearest[2]={1e300,1e300};std::size_t inside=0;
    for(const auto& root:roots) {
      for(int end=0;end<2;++end) {
        B d;acb_sub_si(d.raw(),root.raw(),end,bits);
        if(acb_contains_zero(d.raw()))throw std::domain_error(std::string("continuation: the ")+(end?"requested point":"Euclidean start")+
          " is a singular point of the line's differential equation; choose another euclidean_direction or scale");
        arb_t a;arb_init(a);acb_abs(a,d.raw(),bits);nearest[end]=std::min(nearest[end],arf_get_d(arb_midref(a),ARF_RND_NEAR));arb_clear(a);
      }
      if(arb_contains_zero(acb_imagref(root.raw())) && arb_is_positive(acb_realref(root.raw())) && arf_cmp_si(arb_midref(acb_realref(root.raw())),1)<0)++inside;
    }
    if(progress)*progress<<"Continuation line: "<<roots.size()<<" singular points, "<<inside<<" of them real in (0,1); nearest to the start at distance "
      <<nearest[0]<<", to the requested point at "<<nearest[1]<<".\n"<<std::flush;
  }
  LaurentBoundary end{0,{}};
  for(slong precision=bits;;precision*=2) {
    try {
      kernel::ComplexBall::set_precision(precision);
      end=transport_laurent(gauged,z,vertices,taylor_order);
      break;
    } catch(const std::domain_error& error) {
      if(precision>=8192)throw;
      if(progress)*progress<<"Continuation: transport at "<<precision<<" bits failed ("<<error.what()<<"); repeating at "<<2*precision<<" bits.\n"<<std::flush;
    }
  }
  // Laurent coefficients of the target rows.
  std::int64_t vmin=0;
  struct Series {std::int64_t low;std::vector<B> values;};
  std::vector<std::vector<std::optional<Series>>> rows(targets.size());
  for(std::size_t r=0;r<targets.size();++r)for(std::size_t j=0;j<n;++j) {
    const auto& entry=targets[r].at(j);auto v=exact_epsilon_valuation(entry,1);
    if(!v){rows[r].emplace_back();continue;}
    vmin=std::min(vmin,*v);
    const int needed=top-(lowZ+static_cast<int>(shifts[j]))-static_cast<int>(*v);
    if(needed<0){rows[r].emplace_back();continue;}
    auto regular=epsilon_gauge_detail::multiply_power(entry,1,-*v);
    auto series=feynman::scalar_functional_detail::epsilon_series(regular,1,static_cast<unsigned>(needed));
    Series s{*v,{}};
    for(const auto& c:series) {
      B b;const auto q=c.rational();acb_set_fmpq(b.raw(),q.raw(),bits);s.values.push_back(std::move(b));
    }
    rows[r].push_back(std::move(s));
  }
  int low=lowZ+static_cast<int>(smin)+static_cast<int>(vmin);if(low>top)low=top;
  result.low=low;
  for(std::size_t r=0;r<targets.size();++r) {
    std::vector<B> out(static_cast<std::size_t>(top-low+1),B(0));
    for(std::size_t j=0;j<n;++j) {
      if(!rows[r][j])continue;const auto& s=*rows[r][j];
      const int ylow=lowZ+static_cast<int>(shifts[j]);
      for(std::size_t m=0;m<s.values.size();++m) {
        if(acb_is_zero(s.values[m].raw()))continue;
        const int tp=static_cast<int>(s.low)+static_cast<int>(m);
        for(int k=low;k<=top;++k) {
          const int yp=k-tp;
          if(yp<ylow || yp>highZ+static_cast<int>(shifts[j]))continue;
          B term;acb_mul(term.raw(),s.values[m].raw(),end.values[j][static_cast<std::size_t>(yp-ylow)].raw(),bits);
          acb_add(out[static_cast<std::size_t>(k-low)].raw(),out[static_cast<std::size_t>(k-low)].raw(),term.raw(),bits);
        }
      }
    }
    result.values.push_back(std::move(out));
  }
  // Leading orders that are zero within their radii in every row.
  while(result.low<top) {
    bool zero=true;for(const auto& row:result.values)zero=zero && acb_contains_zero(row.front().raw());
    if(!zero)break;
    for(auto& row:result.values)row.erase(row.begin());
    ++result.low;
  }
  return result;
}
} // namespace diffexp::continuation

namespace diffexp::ibp {
// The family with Gram matrix and masses start + x (end - start).
inline QuadraticFamily line_family(const feynman::Family& start,const feynman::Family& end,const Exact& x,std::size_t physical_count) {
  std::vector<std::vector<Exact>> gram;
  for(std::size_t i=0;i<start.external_gram.size();++i) {
    gram.emplace_back();
    for(std::size_t j=0;j<start.external_gram.size();++j)
      gram.back().push_back(x.constant(start.external_gram[i][j])+x*x.constant(end.external_gram[i][j]-start.external_gram[i][j]));
  }
  ScalarProducts sp(start.loops,std::move(gram),x);
  QuadraticFamily result{sp,{},{}};
  result.external_coefficients.emplace();
  for(std::size_t n=0;n<start.lines.size();++n) {
    const auto& line=start.lines[n];
    result.external_coefficients->emplace_back();
    for(const auto& c:line.external_coefficients)result.external_coefficients->back().push_back(x.constant(c));
    auto d=sp.zero();
    d.constant=x.constant(line.mass_squared)+x*x.constant(end.lines[n].mass_squared-line.mass_squared);
    auto c=line.loop_coefficients;c.insert(c.end(),line.external_coefficients.begin(),line.external_coefficients.end());
    for(unsigned i=0;i<c.size();++i)for(unsigned j=0;j<c.size();++j)d=plus(d,scaled(sp.dot(i,j),x.constant(Rational(0)-c[i]*c[j])));
    if(result.physical.size()<physical_count)result.physical.push_back(std::move(d));
    else result.auxiliary.push_back(std::move(d));
  }
  return result;
}
} // namespace diffexp::ibp

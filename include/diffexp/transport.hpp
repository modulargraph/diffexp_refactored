#pragma once
#include "diffexp/cli.hpp"
#include "diffexp/geometry.hpp"
#include "diffexp/frobenius.hpp"
#include "diffexp/canonical.hpp"
#include <chrono>
#include <memory>

namespace diffexp::transport {
namespace json=boost::json;
using B=kernel::ComplexBall;
inline std::string string(const json::value& v) {return std::string(v.as_string());}
inline void known_keys(const json::object& o,std::initializer_list<std::string_view> keys,const char* scope) {
  for(const auto& kv:o)if(std::find(keys.begin(),keys.end(),std::string_view(kv.key()))==keys.end())throw std::invalid_argument(std::string("unknown ")+scope+" field: "+std::string(kv.key()));
}
inline long integer(const json::object& o,const char* key,long fallback,long low,long high) {
  auto p=o.if_contains(key);
  if(!p){if(fallback<low || fallback>high)throw std::invalid_argument(std::string("missing transport ")+key);return fallback;}
  if(p->is_uint64() && high>=0 && p->as_uint64()<=static_cast<std::uint64_t>(high) && (low<=0 || p->as_uint64()>=static_cast<std::uint64_t>(low)))return static_cast<long>(p->as_uint64());
  if(!p->is_int64() || p->as_int64()<low || p->as_int64()>high)throw std::invalid_argument(std::string("invalid transport ")+key);
  return p->as_int64();
}
inline std::string expression_key(const data::Expr& e) {
  std::string s=e.head+"[";for(const auto& a:e.args)s+=expression_key(a)+",";return s+"]";
}
struct Algebra {
  std::vector<data::Expr> radicands;
  std::map<std::string,unsigned> root_ids;
  data::Expr extract(data::Expr e) {
    for(auto& a:e.args)a=extract(std::move(a));
    if(e.head=="Sqrt" && e.args.size()==1) {
      auto key=expression_key(e.args[0]);auto [it,inserted]=root_ids.try_emplace(key,radicands.size());
      if(inserted){if(radicands.size()>=32)throw std::invalid_argument("transport algebraic root budget exceeded");radicands.push_back(e.args[0]);}
      return {"r"+std::to_string(it->second),{}};
    }
    return e;
  }
};
struct Entry {unsigned row,column,epsilon;Exact coefficient;};
struct Compiled {
  ExactField field;
  std::vector<Exact> squares,root_derivatives;
  std::vector<Entry> entries;
  std::vector<Exact> letters,basis_prefactors;
  std::vector<CanonicalEntry> canonical_entries;
  bool canonical=true;
  std::vector<B> singularities;
  std::vector<Jet> square_polynomials;
  unsigned dimension,epsilon_order;
  Compiled(std::vector<std::string> names,unsigned d,unsigned k):field(std::move(names)),dimension(d),epsilon_order(k){}
  Exact reduced(const Exact& q) const {
    auto result=q;
    for(unsigned r=squares.size();r-->0;)result=reduce_square(result.numerator(),2+r,squares[r])/reduce_square(result.denominator(),2+r,squares[r]);
    return reduce_square(result.numerator(),1,result.constant(-1))/reduce_square(result.denominator(),1,result.constant(-1));
  }
  Exact derivative(const Exact& q) const {
    auto out=q.derivative(0);for(unsigned r=0;r<root_derivatives.size();++r)out=out+q.derivative(2+r)*root_derivatives[r];return reduced(out);
  }
  Exact norm(Exact p) const {
    for(unsigned r=squares.size();r-->0;)p=polynomial_norm(p.numerator(),2+r,squares[r]);
    return polynomial_norm(p.numerator(),1,p.constant(-1));
  }
};
inline Compiled compile(const json::object& request,unsigned d,unsigned k,slong bits) {
  Algebra algebra;std::map<std::string,data::Expr> paths;
  for(const auto& kv:request.at("paths").as_object()) {
    auto name=std::string(kv.key());
    if(name=="x" || name=="I" || (name.size()>1 && name[0]=='r' && std::all_of(name.begin()+1,name.end(),[](unsigned char c){return std::isdigit(c);})))throw std::invalid_argument("path variable collides with a reserved native symbol: "+name);
  }
  for(const auto& kv:request.at("paths").as_object())paths.emplace(std::string(kv.key()),algebra.extract(data::Reader(string(kv.value())).read()));
  std::vector<data::Expr> expressions;
  for(const auto& v:request.at("entries").as_array())expressions.push_back(algebra.extract(data::Reader(string(v.as_object().at("expression"))).read()));
  std::vector<data::Expr> prefactors;
  if(auto supplied=request.if_contains("basis_prefactors")) {
    if(!supplied->is_array() || supplied->as_array().size()!=d)throw std::invalid_argument("basis_prefactors requires one expression string per component");
    for(const auto& value:supplied->as_array()) {
      if(!value.is_string())throw std::invalid_argument("basis_prefactors entries must be expression strings");
      prefactors.push_back(algebra.extract(data::Reader(string(value)).read()));
    }
  }
  std::vector<std::string> names{"x","I"};for(unsigned r=0;r<algebra.radicands.size();++r)names.push_back("r"+std::to_string(r));
  Compiled out(names,d,k);Exact x(out.field,"x"),imag(out.field,"I");
  std::map<std::string,Exact> vars{{"x",x},{"I",imag}};
  for(unsigned r=0;r<algebra.radicands.size();++r)vars.emplace("r"+std::to_string(r),x.variable(2+r));
  for(const auto& [name,e]:paths)vars.insert_or_assign(name,evaluate_exact(e,x,vars));
  for(const auto& e:algebra.radicands) {
    auto square=evaluate_exact(e,x,vars);out.squares.push_back(square);
    auto derivative=out.derivative(square)/(x.constant(2)*x.variable(out.squares.size()+1));out.root_derivatives.push_back(std::move(derivative));
  }
  for(const auto& prefactor:prefactors)out.basis_prefactors.push_back(out.reduced(evaluate_exact(prefactor,x,vars)));
  unsigned index=0;std::map<std::pair<std::string,std::string>,Exact> prepared_coefficients;std::map<std::string,unsigned> letter_ids;
  for(const auto& v:request.at("entries").as_array()) {
    const auto& o=v.as_object();known_keys(o,{"row","column","epsilon","variable","expression","coefficient"},"matrix entry");unsigned row=integer(o,"row",-1,0,d-1),col=integer(o,"column",-1,0,d-1),eps=integer(o,"epsilon",0,0,100);
    if(row>=d || col>=d)throw std::invalid_argument("matrix position exceeds dimension");
    const auto variable=string(o.at("variable"));auto key=std::make_pair(variable,expression_key(expressions[index]));
    auto found=prepared_coefficients.find(key);Exact coefficient(x.constant(0));
    if(found!=prepared_coefficients.end())coefficient=found->second;
    else {
      coefficient=evaluate_exact(expressions[index],x,vars);
      if(variable=="dlog")coefficient=out.derivative(coefficient)/coefficient;
      else {auto it=vars.find(variable);if(it==vars.end())throw std::invalid_argument("matrix variable absent from path: "+variable);coefficient=coefficient*out.derivative(it->second);}
      coefficient=out.reduced(coefficient);prepared_coefficients.emplace(std::move(key),coefficient);
    }
    if(variable=="dlog" && eps==1) {
      const auto letter_key=expression_key(expressions[index]);auto [it,added]=letter_ids.try_emplace(letter_key,out.letters.size());
      if(added)out.letters.push_back(evaluate_exact(expressions[index],x,vars));
      out.canonical_entries.push_back({row,col,it->second,o.if_contains("coefficient")?Rational(string(o.at("coefficient"))):Rational(1)});
    } else out.canonical=false;
    ++index;
    if(auto p=o.if_contains("coefficient"))coefficient=coefficient*Exact(out.field,string(*p));
    if(!coefficient.is_zero() && eps<=k)out.entries.push_back({row,col,eps,std::move(coefficient)});
  }
  if(!out.canonical) {
    std::map<std::tuple<unsigned,unsigned,unsigned>,Exact> combined;
    for(const auto& e:out.entries){auto [it,inserted]=combined.try_emplace({e.row,e.column,e.epsilon},x.constant(0));it->second=it->second+e.coefficient;}
    out.entries.clear();for(auto& [key,q]:combined){q=out.reduced(q);if(!q.is_zero()){auto [row,col,epsilon]=key;out.entries.push_back({row,col,epsilon,std::move(q)});}}
  }
  std::set<std::string> polynomials,raw_polynomials;
  const auto add_roots=[&](const Exact& q) {
    if(q.is_rational() || !raw_polynomials.insert(q.str()).second)return;
    auto p=out.norm(q);if(p.is_rational())return;
    p=p/p.constant(p.numerator_terms().front().coefficient);if(!polynomials.insert(p.str()).second)return;
    auto roots=polynomial_roots(p,0,bits);out.singularities.insert(out.singularities.end(),roots.begin(),roots.end());
  };
  for(const auto& square:out.squares){add_roots(square.numerator());add_roots(square.denominator());}
  for(const auto& e:out.entries)add_roots(e.coefficient.denominator());
  // Independent polynomial radicands cover the original algebraic alphabets.
  // Reject unsupported towers instead of resetting principal branches per chart.
  for(const auto& square:out.squares) {
    if(!square.denominator().is_rational())throw std::invalid_argument("transport currently requires polynomial square-root radicands");
    unsigned degree=0;
    for(const auto& term:square.numerator_terms()) {
      for(unsigned r=2;r<term.powers.size();++r)if(term.powers[r])throw std::invalid_argument("nested algebraic transport roots are not yet supported");
      degree=std::max(degree,static_cast<unsigned>(term.powers[0]));
    }
    Jet t(0,std::max(2u,degree+1),bits);t.set(1,B(1));
    out.square_polynomials.push_back(evaluate(data::Reader(square.str()).read(),t,{{"x",t}}));
  }
  // Norms also contain zeros on conjugate algebraic sheets. A real-path
  // candidate is discarded only when every actual denominator and radicand
  // excludes zero on the continuously transported sheet at that candidate.
  out.singularities.erase(std::remove_if(out.singularities.begin(),out.singularities.end(),[&](const B& candidate) {
    if(!arb_contains_zero(acb_imagref(candidate.raw())))return false;
    // Root isolation may round a real endpoint root to a ball straddling 0 or
    // 1. Check its intersection with the path, not whole-ball containment.
    B point,domain=B::from_strings("1/2");arb_add_error_2exp_si(acb_realref(domain.raw()),-1);
    if(!arb_intersection(acb_realref(point.raw()),acb_realref(candidate.raw()),acb_realref(domain.raw()),bits))return false;
    Jet context(0,1,bits);context.set(0,point);std::map<std::string,Jet> env{{"x",context}};
    for(unsigned root=0;root<out.square_polynomials.size();++root) {
      auto square=out.square_polynomials[root].evaluate_polynomial(point);if(square.contains_zero())return false;
      auto start=out.square_polynomials[root].evaluate_polynomial(B(0));if(start.contains_zero())return false;
      B value;acb_sqrt(value.raw(),start.raw(),bits);
      try{value=continue_polynomial_sqrt(out.square_polynomials[root],B(0),point,value);}catch(const std::exception&){return false;}
      auto j=context.constant(0);j.set(0,value);env.emplace("r"+std::to_string(root),j);
    }
    for(const auto& entry:out.entries)
      if(evaluate(data::Reader(entry.coefficient.denominator().str()).read(),context,env).at(0).contains_zero())return false;
    return true;
  }),out.singularities.end());
  return out;
}
inline std::string decimal_mid(const arb_t a,long digits) {
  arb_t point;arb_init(point);arb_set_arf(point,arb_midref(a));char* p=arb_get_str(point,digits,ARB_STR_NO_RADIUS);
  std::string s(p);flint_free(p);arb_clear(point);return s;
}
inline json::object number(const B& b,long digits) {
  if(!b.is_finite())throw std::runtime_error("non-finite transport output");
  arb_t radius;arb_init(radius);mag_t mag;mag_init(mag);mag_add(mag,arb_radref(acb_realref(b.raw())),arb_radref(acb_imagref(b.raw())));arf_set_mag(arb_midref(radius),mag);
  auto s=decimal_mid(radius,digits);arb_clear(radius);mag_clear(mag);
  return {{"real_midpoint",decimal_mid(acb_realref(b.raw()),digits)},{"imaginary_midpoint",decimal_mid(acb_imagref(b.raw()),digits)},{"radius",s}};
}
inline json::array matrix_json(const Boundary& values,long digits) {
  json::array rows;for(const auto& row:values){json::array cells;for(const auto& b:row)cells.push_back(number(b,digits));rows.push_back(std::move(cells));}return rows;
}
inline B magnitude(const B& v) {B b;acb_abs(acb_realref(b.raw()),v.raw(),B::precision());return b;}
inline B arithmetic_error(const B& v) {B b;mag_t m;mag_init(m);mag_add(m,arb_radref(acb_realref(v.raw())),arb_radref(acb_imagref(v.raw())));arf_set_mag(arb_midref(acb_realref(b.raw())),m);mag_clear(m);return b;}
// Endpoint normalization is supplied mathematical basis data. It is never
// inferred from a component's value, a reference result, or a family name.
inline std::vector<B> principal_roots_at(const Compiled& c,const B& point) {
  std::vector<Exact> replacements;
  if(point.is_finite() && acb_is_exact(point.raw()) && !c.squares.empty()) {
    const auto rational_string=[](const arf_t coordinate) {
      fmpq_t q;fmpq_init(q);arf_get_fmpq(q,coordinate);
      char* raw=fmpq_get_str(nullptr,10,q);std::string text(raw);flint_free(raw);fmpq_clear(q);return text;
    };
    const auto& model=c.squares.front();
    for(unsigned i=0;i<model.variables().size();++i)replacements.push_back(model.variable(i));
    replacements[0]=Exact(c.field,rational_string(arb_midref(acb_realref(point.raw()))))+
      model.variable(1)*Exact(c.field,rational_string(arb_midref(acb_imagref(point.raw()))));
  }
  std::vector<B> roots;for(unsigned r=0;r<c.square_polynomials.size();++r) {
    B square;
    if(!replacements.empty()) {
      // Cancel endpoint polynomial coefficients exactly before rounding. Tiny
      // artificial imaginary intervals on a negative real radicand straddle
      // the principal square-root cut and otherwise destroy its precision.
      auto exact=c.reduced(c.squares[r].substitute(replacements));Jet context(0,1,B::precision());
      square=evaluate(data::Reader(exact.str()).read(),context,{}).at(0);
    } else square=c.square_polynomials[r].evaluate_polynomial(point);
    B root;acb_sqrt(root.raw(),square.raw(),B::precision());
    if(!root.is_finite())throw std::invalid_argument("non-finite principal endpoint root");roots.push_back(std::move(root));
  }return roots;
}
inline std::vector<B> basis_prefactor_ratios(const Compiled& c,const B& endpoint,const std::vector<B>& continued_roots) {
  if(c.basis_prefactors.empty())return {};
  if(continued_roots.size()!=c.square_polynomials.size())throw std::invalid_argument("basis prefactor root dimensions");
  Jet point(0,1,B::precision());point.set(0,endpoint);std::map<std::string,Jet> continued{{"x",point}},principal{{"x",point}};
  auto principal_roots=principal_roots_at(c,endpoint);
  for(unsigned r=0;r<continued_roots.size();++r) {
    if(!continued_roots[r].is_finite())throw std::invalid_argument("non-finite continued endpoint root");
    auto a=point.constant(0),b=point.constant(0);a.set(0,continued_roots[r]);b.set(0,principal_roots[r]);
    continued.emplace("r"+std::to_string(r),std::move(a));principal.emplace("r"+std::to_string(r),std::move(b));
  }
  std::vector<B> ratios;for(const auto& prefactor:c.basis_prefactors) {
    auto expression=data::Reader(prefactor.str()).read();auto denominator=evaluate(expression,point,continued).at(0),numerator=evaluate(expression,point,principal).at(0);
    if(!denominator.is_finite() || denominator.contains_zero() || !numerator.is_finite())
      throw std::invalid_argument("basis prefactor has a zero-containing or non-finite endpoint divisor/value");
    auto ratio=numerator/denominator;if(!ratio.is_finite())throw std::invalid_argument("non-finite basis prefactor ratio");ratios.push_back(std::move(ratio));
  }return ratios;
}
// Reusable for a completed transport: reconstruct the continued roots from its
// exact path, then apply this helper without repeating the differential equation.
inline void apply_basis_prefactors(const Compiled& c,const B& endpoint,const std::vector<B>& continued_roots,
    Boundary& values,Boundary& errors) {
  if(c.basis_prefactors.empty())return;
  if(values.size()!=c.dimension || errors.size()!=c.dimension)throw std::invalid_argument("basis prefactor boundary dimensions");
  for(unsigned i=0;i<c.dimension;++i)if(values[i].size()!=c.epsilon_order+1 || errors[i].size()!=c.epsilon_order+1)
    throw std::invalid_argument("basis prefactor epsilon window");
  auto ratios=basis_prefactor_ratios(c,endpoint,continued_roots);
  // Stage changes so a rejected divisor or output never partially converts a boundary.
  auto converted=values,estimates=errors;
  for(unsigned i=0;i<c.dimension;++i)for(unsigned e=0;e<=c.epsilon_order;++e) {
    if(!errors[i][e].is_finite() || !arb_is_zero(acb_imagref(errors[i][e].raw())) || !arb_is_nonnegative(acb_realref(errors[i][e].raw())))
      throw std::invalid_argument("basis prefactor errors must be finite real nonnegative bounds");
    converted[i][e]=ratios[i]*values[i][e];
    if(!converted[i][e].is_finite())throw std::invalid_argument("non-finite converted basis value");
    estimates[i][e]=magnitude(ratios[i])*errors[i][e]+arithmetic_error(converted[i][e]);
    if(!estimates[i][e].is_finite())throw std::invalid_argument("non-finite converted basis error");
  }
  values=std::move(converted);errors=std::move(estimates);
}
inline Boundary boundary(const json::object& r,unsigned d,unsigned k,slong bits) {
  Boundary out(d,std::vector<B>(k+1,B(0)));const auto& input=r.at("boundary").as_array();
  if(input.size()!=d)throw std::invalid_argument("transport boundary dimension");Jet ctx(0,1,bits);
  for(unsigned i=0;i<d;++i) {
    if(input[i].as_array().size()!=k+1)throw std::invalid_argument("transport boundary epsilon window");
    for(unsigned e=0;e<=k;++e)out[i][e]=evaluate(data::Reader(string(input[i].as_array()[e])).read(),ctx,{}).at(0);
  }
  if(auto p=r.if_contains("boundary_errors")) {
    if(p->as_array().size()!=d)throw std::invalid_argument("boundary error dimension");
    for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e) {
      auto error=evaluate(data::Reader(string(p->as_array().at(i).as_array().at(e))).read(),ctx,{}).at(0);
      if(!arb_is_zero(acb_imagref(error.raw())) || !arb_is_nonnegative(acb_realref(error.raw())))throw std::invalid_argument("boundary errors must be real and nonnegative");
      arb_add_error(acb_realref(out[i][e].raw()),acb_realref(error.raw()));arb_add_error(acb_imagref(out[i][e].raw()),acb_realref(error.raw()));
    }
  }
  return out;
}
struct NumericalEntry {unsigned epsilon;std::vector<std::pair<unsigned,unsigned>> positions;data::Expr coefficient;};
inline std::vector<NumericalEntry> numerical_entries(const Compiled& c) {
  std::vector<NumericalEntry> out;std::map<std::pair<unsigned,std::string>,unsigned> ids;
  for(const auto& e:c.entries) {
    auto key=std::make_pair(e.epsilon,e.coefficient.str());auto [it,inserted]=ids.try_emplace(key,out.size());
    if(inserted)out.push_back({e.epsilon,{},data::Reader(key.second).read()});out[it->second].positions.emplace_back(e.row,e.column);
  }
  return out;
}
// Only feedback in the epsilon-zero dependency graph repeatedly recirculates
// the same uncertainty. Triangular iterated-integral systems keep the direct
// enclosing recurrence; a dense homogeneous map is unnecessary work there.
inline bool zero_order_feedback(const Compiled& c) {
  std::vector<std::vector<unsigned>> outgoing(c.dimension);std::vector<unsigned> incoming(c.dimension,0),ready;
  for(const auto& e:c.entries)if(e.epsilon==0) {
    if(e.row==e.column)return true;
    outgoing[e.column].push_back(e.row);++incoming[e.row];
  }
  for(unsigned i=0;i<c.dimension;++i)if(incoming[i]==0)ready.push_back(i);
  for(std::size_t at=0;at<ready.size();++at)for(auto row:outgoing[ready[at]])if(--incoming[row]==0)ready.push_back(row);
  return ready.size()!=c.dimension;
}
struct AcbArray {
  slong size;acb_ptr p;explicit AcbArray(slong n):size(n),p(_acb_vec_init(n)){}~AcbArray(){_acb_vec_clear(p,size);}AcbArray(const AcbArray&)=delete;
};
struct Chart {Boundary values,errors,truncation_errors;json::object saved;bool centered=false;};
inline Chart chart(const Compiled& c,const std::vector<NumericalEntry>& entries,const Boundary& initial,
    const std::vector<B>& root_values,const B& center,const B& step,unsigned order,bool save,long digits,
    bool center_uncertainty=true,std::size_t map_cell_budget=10000000,bool share_canonical=true) {
  const auto bits=B::precision();const unsigned d=c.dimension,w=c.epsilon_order+1;
  // Separate inherited uncertainty before the recurrence can repeatedly use it.
  // This selector changes only the enclosing arithmetic, never accuracy controls.
  bool centered=false;B floor(1);acb_mul_2exp_si(floor.raw(),floor.raw(),-bits/2);
  const bool coupled_at_zero=!c.canonical && zero_order_feedback(c);
  if(center_uncertainty && !c.canonical && coupled_at_zero)for(const auto& row:initial)for(const auto& value:row)
    if(!arb_le(acb_realref(arithmetic_error(value).raw()),acb_realref((floor*(magnitude(value)+B(1))).raw())))centered=true;
  const std::size_t base_cells=static_cast<std::size_t>(order+1)*d*w;
  // A large system retains the original interval recurrence when a local map
  // would exceed this finite workspace budget. No input radius is removed.
  if(centered && (base_cells>map_cell_budget || d+1>map_cell_budget/base_cells))centered=false;
  const unsigned columns=centered?d+1:1;
  Boundary delta(d,std::vector<B>(w,B(0)));
  Jet x(0,order+2,bits);x.set(0,center);x.set(1,B(1));std::map<std::string,Jet> vars{{"x",x}};
  for(unsigned r=0;r<c.squares.size();++r) {
    auto square=evaluate(data::Reader(c.squares[r].str()).read(),x,vars),root=square.sqrt();
    auto factor=root.constant(0);factor.set(0,root_values[r]/root.at(0));root=root*factor;vars.emplace("r"+std::to_string(r),std::move(root));
  }
  AcbArray coeffs(static_cast<slong>(c.canonical?c.letters.size():entries.size())*order),values(static_cast<slong>(base_cells)*columns);
  if(c.canonical) {
    for(unsigned l=0;l<c.letters.size();++l){auto letter=evaluate(data::Reader(c.letters[l].str()).read(),x,vars);auto jet=letter.derivative()/letter;
      for(unsigned n=0;n<order;++n){auto value=jet.at(n);acb_swap(coeffs.p+l*order+n,value.raw());}}
  } else for(unsigned e=0;e<entries.size();++e) {auto jet=evaluate(entries[e].coefficient,x,vars);for(unsigned n=0;n<order;++n){auto value=jet.at(n);acb_swap(coeffs.p+e*order+n,value.raw());}}
  const slong stride=static_cast<slong>(columns)*d*w;
  const auto at=[&](unsigned n,unsigned column,unsigned i,unsigned e){return values.p+n*stride+(static_cast<slong>(column)*d+i)*w+e;};
  for(unsigned i=0;i<d;++i)for(unsigned e=0;e<w;++e) {
    if(centered) {
      acb_get_mid(at(0,0,i,e),initial[i][e].raw());
      acb_sub(delta[i][e].raw(),initial[i][e].raw(),at(0,0,i,e),bits);
    } else acb_set(at(0,0,i,e),initial[i][e].raw());
  }
  if(centered)for(unsigned i=0;i<d;++i)acb_one(at(0,i+1,i,0));
  B sum;std::vector<B> weights;for(const auto& e:c.canonical_entries)weights.push_back(B::from_strings(e.coefficient.str()));
  // Rows often reuse the same letter/source convolution. Cache that identical
  // dot result at each Taylor/epsilon order; preserve every row accumulation
  // and its order, including its exact rational weight.
  std::map<std::pair<unsigned,unsigned>,unsigned> pair_ids;std::vector<unsigned> entry_pairs;
  if(c.canonical && share_canonical)for(const auto& e:c.canonical_entries) {
    auto [it,inserted]=pair_ids.try_emplace({e.letter,e.column},pair_ids.size());entry_pairs.push_back(it->second);
  }
  AcbArray convolutions(pair_ids.size()*w);std::vector<int> convolution_order(pair_ids.size()*w,-1);
  for(unsigned n=0;n<order;++n)for(unsigned column=0;column<columns;++column) {
    if(c.canonical) {
      for(unsigned p=0;p<c.canonical_entries.size();++p) {const auto& entry=c.canonical_entries[p];for(unsigned e=1;e<w;++e){
        acb_srcptr convolution;
        if(share_canonical) {
          const auto id=entry_pairs[p]*w+e;
          if(convolution_order[id]!=static_cast<int>(n)) {
            acb_dot(convolutions.p+id,nullptr,0,coeffs.p+entry.letter*order,1,at(n,column,entry.column,e-1),-stride,n+1,bits);
            convolution_order[id]=n;
          }
          convolution=convolutions.p+id;
        } else {
          acb_dot(sum.raw(),nullptr,0,coeffs.p+entry.letter*order,1,at(n,column,entry.column,e-1),-stride,n+1,bits);
          convolution=sum.raw();
        }
        acb_addmul(at(n+1,column,entry.row,e),convolution,weights[p].raw(),bits);
      }}
    } else for(unsigned p=0;p<entries.size();++p)for(const auto& [i,j]:entries[p].positions)for(unsigned e=entries[p].epsilon;e<w;++e) {
      auto target=at(n+1,column,i,e);acb_dot(sum.raw(),target,0,coeffs.p+p*order,1,at(n,column,j,e-entries[p].epsilon),-stride,n+1,bits);acb_swap(target,sum.raw());
    }
    for(unsigned i=0;i<d;++i)for(unsigned e=0;e<w;++e)acb_div_ui(at(n+1,column,i,e),at(n+1,column,i,e),n+1,bits);
  }
  Chart result{Boundary(d,std::vector<B>(w,B(0))),Boundary(d,std::vector<B>(w,B(0))),Boundary(d,std::vector<B>(w,B(0))),{},centered};
  const auto coefficient=[&](unsigned n,unsigned i,unsigned e) {
    B value;acb_set(value.raw(),at(n,0,i,e));
    if(centered)for(unsigned j=0;j<d;++j)for(unsigned l=0;l<=e;++l)
      acb_addmul(value.raw(),at(n,j+1,i,l),delta[j][e-l].raw(),bits);
    return value;
  };
  const auto evaluate_column=[&](unsigned column,unsigned i,unsigned e) {
    B value;for(unsigned n=order+1;n-->0;){acb_mul(value.raw(),value.raw(),step.raw(),bits);acb_add(value.raw(),value.raw(),at(n,column,i,e),bits);}return value;
  };
  std::vector<B> step_powers(order+1,B(1));for(unsigned n=1;n<=order;++n)step_powers[n]=step_powers[n-1]*step;
  for(unsigned i=0;i<d;++i)for(unsigned e=0;e<w;++e) {
    auto& v=result.values[i][e];v=evaluate_column(0,i,e);
    // Evaluate the map first: cancellations between Taylor coefficients must
    // occur before the same input uncertainty is applied, once per chart.
    if(centered)for(unsigned j=0;j<d;++j)for(unsigned l=0;l<=e;++l) {
      auto map=evaluate_column(j+1,i,l);acb_addmul(v.raw(),map.raw(),delta[j][e-l].raw(),bits);
    }
    if(!v.is_finite())throw std::runtime_error("transport chart produced a non-finite value");
    // Preserve the existing last-four-term estimate for the entire uncertain
    // input, not just its midpoint. This is still not a remainder certificate.
    for(unsigned n=order>3?order-3:1;n<=order;++n) {
      auto term=coefficient(n,i,e)*step_powers[n];result.errors[i][e]+=magnitude(term);
    }
    result.truncation_errors[i][e]=result.errors[i][e];
    result.errors[i][e]+=arithmetic_error(v);
  }
  if(save) {
    json::array all;for(unsigned n=0;n<=order;++n){Boundary m(d,std::vector<B>(w));for(unsigned i=0;i<d;++i)for(unsigned e=0;e<w;++e)m[i][e]=coefficient(n,i,e);all.push_back(matrix_json(m,digits));}
    result.saved={{"center",decimal_mid(acb_realref(center.raw()),digits)},{"start",decimal_mid(acb_realref(center.raw()),digits)},
      {"end",decimal_mid(acb_realref((center+step).raw()),digits)},{"coefficients",std::move(all)}};
  }
  return result;
}
// Solve asymptotic coefficient constraints in the exact Frobenius frame. Each
// supplied row fixes every power/log coefficient below its explicit cutoff.
inline Boundary asymptotic_constants(const FrobeniusSeries& frame,const json::object& asymptotic) {
  const auto d=frame.dimension(),kmax=frame.epsilon_order(),size=d*(kmax+1);Jet ctx(0,1,B::precision());
  using Key=std::tuple<unsigned,unsigned,Rational,unsigned>;
  std::map<unsigned,Rational> cutoffs;for(const auto& v:asymptotic.at("cutoffs").as_array()) {
    auto& o=v.as_object();auto row=integer(o,"row",-1,0,d-1);cutoffs.insert_or_assign(row,Rational(string(o.at("power"))));
  }
  std::map<Key,B> expected;
  for(const auto& v:asymptotic.at("constraints").as_array()) {
    auto& o=v.as_object();unsigned row=integer(o,"row",-1,0,d-1),e=integer(o,"epsilon",0,0,kmax),l=integer(o,"log_degree",0,0,1000);
    auto key=Key{row,e,Rational(string(o.at("power"))),l};auto [it,inserted]=expected.try_emplace(key,B(0));
    it->second+=evaluate(data::Reader(string(o.at("value"))).read(),ctx,{}).at(0);
  }
  std::map<Key,std::vector<Rational>> equations;
  for(const auto& m:frame.monomials()) {
    auto cutoff=cutoffs.find(m.row);
    for(unsigned e=m.epsilon;e<=kmax;++e) {
      auto key=Key{m.row,e,m.power,m.log_degree};
      if(!(cutoff!=cutoffs.end() && m.power<cutoff->second) && !expected.contains(key))continue;
      auto [it,inserted]=equations.try_emplace(key,size,Rational(0));
      it->second[m.column*(kmax+1)+e-m.epsilon]+=m.coefficient;
    }
  }
  for(const auto& [key,v]:expected)equations.try_emplace(key,size,Rational(0));
  // Exact elimination selects independent constraints without a numerical rank
  // tolerance. Numerical right-hand sides carry only supplied/rounding error.
  std::vector<std::vector<Rational>> pivots(size);std::vector<B> rhs(size);unsigned rank=0;
  for(auto& [key,row]:equations) {
    B value=expected.contains(key)?expected.at(key):B(0);
    for(unsigned col=0;col<size;++col)if(!row[col].is_zero()) {
      auto q=row[col];
      if(pivots[col].empty()) {
        for(auto& x:row)x=x/q;value=value/B::from_strings(q.str());pivots[col]=row;rhs[col]=value;++rank;break;
      }
      for(unsigned j=col;j<size;++j)row[j]-=q*pivots[col][j];value-=B::from_strings(q.str())*rhs[col];
    }
    if(std::all_of(row.begin(),row.end(),[](const auto& q){return q.is_zero();}) && !value.contains_zero())
      throw std::invalid_argument("asymptotic boundary constraints are inconsistent");
  }
  if(rank!=size)throw std::invalid_argument("asymptotic boundary leaves "+std::to_string(size-rank)+" integration constants undetermined");
  std::vector<B> solution(size,B(0));for(unsigned col=size;col-->0;) {auto value=rhs[col];for(unsigned j=col+1;j<size;++j)value-=B::from_strings(pivots[col][j].str())*solution[j];solution[col]=value;}
  Boundary out(d,std::vector<B>(kmax+1));for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=kmax;++e)out[i][e]=solution[i*(kmax+1)+e];return out;
}
inline std::optional<std::pair<Boundary,Boundary>> finite_start(const Compiled& c,const json::object& asymptotic,const B& endpoint,unsigned order) {
  if(!c.squares.empty())return std::nullopt;
  const unsigned d=c.dimension,w=c.epsilon_order+1;std::set<unsigned> rows;
  for(const auto& v:asymptotic.at("cutoffs").as_array()) {const auto& o=v.as_object();if(Rational(string(o.at("power")))<=Rational(0) || Rational(string(o.at("power")))>Rational(1))return std::nullopt;rows.insert(integer(o,"row",-1,0,d-1));}
  if(rows.size()!=d)return std::nullopt;
  Boundary initial(d,std::vector<B>(w,B(0)));Jet context(0,order+3,B::precision());context.set(1,B(1));
  for(const auto& v:asymptotic.at("constraints").as_array()) {const auto& o=v.as_object();
    if(Rational(string(o.at("power")))!=Rational(0) || integer(o,"log_degree",0,0,1000)!=0)return std::nullopt;
    initial[integer(o,"row",-1,0,d-1)][integer(o,"epsilon",0,0,w-1)]+=evaluate(data::Reader(string(o.at("value"))).read(),context,{}).at(0);
  }
  struct Local {unsigned row,col,eps;int power;std::vector<B> coefficients;};std::vector<Local> matrix;
  for(const auto& e:c.entries) {
    auto numerator=evaluate(data::Reader(e.coefficient.numerator().str()).read(),context,{{"x",context}}),denominator=evaluate(data::Reader(e.coefficient.denominator().str()).read(),context,{{"x",context}});
    unsigned nv=0,dv=0;while(nv<numerator.length()&&numerator.at(nv).is_zero())++nv;while(dv<denominator.length()&&denominator.at(dv).is_zero())++dv;
    if(nv==numerator.length())continue;if(dv==denominator.length() || static_cast<int>(nv)-static_cast<int>(dv)<-1)throw std::invalid_argument("finite singular start is not regular singular");
    auto q=numerator.shifted_down(nv)/denominator.shifted_down(dv);Local local{e.row,e.column,e.epsilon,static_cast<int>(nv)-static_cast<int>(dv),{}};
    for(unsigned n=0;n<=order+1;++n)local.coefficients.push_back(q.at(n));matrix.push_back(std::move(local));
  }
  for(unsigned k=0;k<w;++k)for(unsigned row=0;row<d;++row) {B residue(0);for(const auto& e:matrix)if(e.row==row&&e.power==-1&&e.eps<=k)residue+=e.coefficients[0]*initial[e.col][k-e.eps];if(!residue.contains_zero())throw std::invalid_argument("boundary is incompatible with a finite singular start");}
  struct Matrices {acb_mat_t a,inv;Matrices(unsigned n){acb_mat_init(a,n,n);acb_mat_init(inv,n,n);}~Matrices(){acb_mat_clear(a);acb_mat_clear(inv);}} m(d);
  std::vector<Boundary> values(order+1,Boundary(d,std::vector<B>(w,B(0))));values[0]=initial;
  // At each order invert n I - R_0 once, sharing it over the epsilon window.
  for(unsigned n=1;n<=order;++n) {
    acb_mat_zero(m.a);for(unsigned i=0;i<d;++i)acb_set_ui(acb_mat_entry(m.a,i,i),n);
    for(const auto& e:matrix)if(e.power==-1 && e.eps==0)acb_sub(acb_mat_entry(m.a,e.row,e.col),acb_mat_entry(m.a,e.row,e.col),e.coefficients[0].raw(),B::precision());
    if(!acb_mat_inv(m.inv,m.a,B::precision()))throw std::invalid_argument("finite boundary does not fix resonant analytic constants; supply more asymptotic terms");
    for(unsigned k=0;k<w;++k) {
      std::vector<B> rhs(d,B(0));
      for(const auto& e:matrix)if(e.eps<=k)for(unsigned j=0;j<e.coefficients.size();++j) {
        int power=e.power+static_cast<int>(j),source=static_cast<int>(n)-1-power;if(source<0)break;
        if(power==-1 && e.eps==0)continue;
        acb_addmul(rhs[e.row].raw(),e.coefficients[j].raw(),values[source][e.col][k-e.eps].raw(),B::precision());
      }
      for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j)acb_addmul(values[n][i][k].raw(),acb_mat_entry(m.inv,i,j),rhs[j].raw(),B::precision());
    }
  }
  Boundary output(d,std::vector<B>(w,B(0))),errors=output;B power(1);std::vector<B> powers(order+1,B(1));for(unsigned n=1;n<=order;++n)powers[n]=powers[n-1]*endpoint;
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<w;++k) {
    for(unsigned n=order+1;n-->0;)output[i][k]=output[i][k]*endpoint+values[n][i][k];
    for(unsigned n=order>3?order-3:1;n<=order;++n)errors[i][k]+=magnitude(values[n][i][k]*powers[n]);errors[i][k]+=arithmetic_error(output[i][k]);
  }
  return std::make_pair(output,errors);
}

inline json::value run(const json::value& input) {
  const auto start=std::chrono::steady_clock::now();const auto& r=input.as_object();
  known_keys(r,{"schema","dimension","epsilon_order","taylor_order","working_bits","accuracy_goal","division_order","save_segments","paths","entries","boundary","boundary_errors","asymptotic","initial_only","basis_prefactors"},"transport");
  if(auto a=r.if_contains("asymptotic")) {
    known_keys(a->as_object(),{"constraints","cutoffs"},"asymptotic");
    for(const auto& v:a->as_object().at("constraints").as_array())known_keys(v.as_object(),{"row","epsilon","power","log_degree","value"},"asymptotic constraint");
    for(const auto& v:a->as_object().at("cutoffs").as_array())known_keys(v.as_object(),{"row","power"},"asymptotic cutoff");
  }
  if(auto p=r.if_contains("schema");p && *p!="DiffExp.Transport/v1")throw std::invalid_argument("unsupported transport schema");
  unsigned d=integer(r,"dimension",0,1,1000),k=integer(r,"epsilon_order",4,0,100),order=integer(r,"taylor_order",50,8,1000);
  if(!d || static_cast<std::uint64_t>(d)*(k+1)*(order+1)>10000000)throw std::invalid_argument("transport coefficient budget exceeded");
  auto bits=integer(r,"working_bits",384,64,100000),division=integer(r,"division_order",4,2,100),goal=integer(r,"accuracy_goal",0,0,20000);
  if(goal*3.32193>bits-32)throw std::invalid_argument("working precision has insufficient reserve for AccuracyGoal");
  B::set_precision(bits);const long digits=bits*30103/100000+5;bool save=r.if_contains("save_segments") && r.at("save_segments").as_bool();
  std::cerr<<"Preparing generic transport: "<<d<<" components, epsilon 0.."<<k<<"\n";
  auto c=compile(r,d,k,bits);const auto prepared=std::chrono::steady_clock::now();
  Boundary current,errors(d,std::vector<B>(k+1,B(0)));double center=0;json::array segments;unsigned charts=0;std::uint64_t saved_budget_used=0;
  if(auto a=r.if_contains("asymptotic")) {
    if(!c.squares.empty())throw std::invalid_argument("asymptotic matching currently requires a rational connection");
    std::vector<B> finite=c.singularities;finite.erase(std::remove_if(finite.begin(),finite.end(),[](const B& p){return p.contains_zero();}),finite.end());
    double next=std::min(0.125,clearance_endpoint(0,finite));next/=4;
    B endpoint;acb_set_d(endpoint.raw(),next);
    // A short frame fixes the constants. Increase the radius precision by using
    // a small initial point; subsequent charts use the physical basis directly.
    const auto evaluate_start=[&]() {
    if(auto finite_result=finite_start(c,a->as_object(),endpoint,std::min(order,40u))) {current=std::move(finite_result->first);errors=std::move(finite_result->second);}
    else {
      ExactField f({"x","eps"});Exact x(f,"x"),eps(f,"eps");FrobeniusSeries::ExactMatrix matrix(d,std::vector<Exact>(d,Exact(f,0)));
      for(const auto& e:c.entries)matrix[e.row][e.column]=matrix[e.row][e.column]+Exact(f,e.coefficient.str())*eps.pow(e.epsilon);
      auto frame=FrobeniusSeries::prepare(matrix,0,1,std::min(order,32u),k);
      auto constants=asymptotic_constants(frame,a->as_object());current=frame.solution(endpoint,constants);
      auto shorter=FrobeniusSeries::prepare(matrix,0,1,std::min(order,32u)-4,k).solution(endpoint,constants);
      for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e)errors[i][e]=magnitude(current[i][e]-shorter[i][e])+arithmetic_error(current[i][e]);
    }
    };
    evaluate_start();
    if(goal>0) {
      B tolerance=B::from_strings("1e-"+std::to_string(goal+4));unsigned attempts=0;
      const auto acceptable=[&](){for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e){auto bound=(magnitude(current[i][e])+B(1))*tolerance;if(!arb_le(acb_realref(errors[i][e].raw()),acb_realref(bound.raw())))return false;}return true;};
      while(!acceptable()){if(++attempts>30)throw std::runtime_error("asymptotic start cannot meet AccuracyGoal with the supplied boundary precision");next/=2;acb_set_d(endpoint.raw(),next);evaluate_start();}
    }
    for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e) {
      arb_add_error(acb_realref(current[i][e].raw()),acb_realref(errors[i][e].raw()));
      arb_add_error(acb_imagref(current[i][e].raw()),acb_realref(errors[i][e].raw()));
    }
    center=next;++charts;
    if(r.if_contains("initial_only") && r.at("initial_only").as_bool()) {
      apply_basis_prefactors(c,endpoint,{},current,errors);
      const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
      return json::object{{"schema","DiffExp.TransportResult/v1"},{"parameter",decimal_mid(acb_realref(endpoint.raw()),digits)},
        {"values",matrix_json(current,digits)},{"errors",matrix_json(errors,digits)},{"charts",charts},{"segments",std::move(segments)},
        {"timings",json::object{{"total_seconds",elapsed},{"preparation_seconds",std::chrono::duration<double>(prepared-start).count()}}},{"omitted_tails_certified",false}};
    }
  } else current=boundary(r,d,k,bits);
  std::vector<B> root_values;B initial_center;acb_set_d(initial_center.raw(),center);
  for(const auto& p:c.square_polynomials){B root;auto value=p.evaluate_polynomial(initial_center);acb_sqrt(root.raw(),value.raw(),bits);root_values.push_back(root);}
  auto entries=numerical_entries(c);
  if(center==0 && std::any_of(c.singularities.begin(),c.singularities.end(),[](const B& p){return p.contains_zero();})) {
    if(!c.canonical)throw std::invalid_argument("a singular initial point requires asymptotic boundary conditions");
    auto finite=c.singularities;finite.erase(std::remove_if(finite.begin(),finite.end(),[](const B& p){return p.contains_zero();}),finite.end());
    auto next=std::min(1.0,clearance_endpoint(0,finite)*4.0/division);B end;acb_set_d(end.raw(),next);
    Jet x(0,order+4,bits);x.set(1,B(1));std::map<std::string,Jet> vars{{"x",x}};
    for(unsigned root=0;root<c.squares.size();++root) {
      auto square=evaluate(data::Reader(c.squares[root].str()).read(),x,vars);auto value=square.sqrt();auto factor=x.constant(0);factor.set(0,root_values[root]/value.at(0));
      vars.emplace("r"+std::to_string(root),value*factor);
    }
    std::vector<Jet> letters;for(const auto& letter:c.letters)letters.push_back(evaluate(data::Reader(letter.str()).read(),x,vars));
    Boundary values,shorter;
    const auto evaluate_initial=[&](){values=canonical_chart(c.canonical_entries,letters,current,B(0),end,order);shorter=canonical_chart(c.canonical_entries,letters,current,B(0),end,order-4);
      for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e)errors[i][e]=magnitude(values[i][e]-shorter[i][e])+arithmetic_error(values[i][e]);};
    evaluate_initial();
    if(goal>0) {
      auto tolerance=B::from_strings("1e-"+std::to_string(goal+4));unsigned attempts=0;
      const auto acceptable=[&](){for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e){auto bound=(magnitude(values[i][e])+B(1))*tolerance;if(!arb_le(acb_realref(errors[i][e].raw()),acb_realref(bound.raw())))return false;}return true;};
      while(!acceptable()){if(++attempts>30)throw std::runtime_error("finite singular start cannot meet AccuracyGoal with the supplied boundary precision");next/=2;acb_set_d(end.raw(),next);evaluate_initial();}
    }
    for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e) {
      arb_add_error(acb_realref(values[i][e].raw()),acb_realref(errors[i][e].raw()));
      arb_add_error(acb_imagref(values[i][e].raw()),acb_realref(errors[i][e].raw()));
    }
    current=std::move(values);
    for(unsigned root=0;root<root_values.size();++root)root_values[root]=continue_polynomial_sqrt(c.square_polynomials[root],B(0),end,root_values[root]);
    center=next;++charts;
  }
  while(center<1) {
    if(++charts>20000)throw std::runtime_error("transport chart budget exhausted");
    auto next=clearance_endpoint(center,c.singularities);next=center+(next-center)*4.0/division;next=std::min(next,1.0);
    if(!(next>center))throw std::runtime_error("transport has no nonzero continuation step");
    B at,end;acb_set_d(at.raw(),center);acb_set_d(end.raw(),next);
    if(save) {
      const std::uint64_t estimate=static_cast<std::uint64_t>(order+1)*d*(k+1)*(3*(digits+32)+80);
      if(estimate>64*1024*1024 || saved_budget_used>64*1024*1024-estimate)throw std::length_error("saved transport segments exceed the 64 MiB output budget");
      saved_budget_used+=estimate;
    }
    auto result=chart(c,entries,current,root_values,at,end-at,order,save,digits);
    if(goal>0) {
      // Bounded local refinement at a fixed expansion order. No retry of the
      // complete transport, and no relaxation of the requested accuracy.
      B tolerance=B::from_strings("1e-"+std::to_string(goal+8));unsigned attempts=0;
      const auto acceptable=[&]() {for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e){auto scale=magnitude(result.values[i][e])+B(1);if(!arb_le(acb_realref(result.truncation_errors[i][e].raw()),acb_realref((scale*tolerance).raw())))return false;}return true;};
      while(!acceptable()) {
        if(++attempts>20)throw std::runtime_error("chart could not meet AccuracyGoal within 20 local refinements");
        next=(center+next)/2;acb_set_d(end.raw(),next);result=chart(c,entries,current,root_values,at,end-at,order,save,digits);
      }
    }
    for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e) {
      errors[i][e]=result.errors[i][e];
      // Carry the estimated omitted terms into the next chart as input radii.
      // This remains an error estimate, not a proof of the infinite remainder.
      auto tail=result.truncation_errors[i][e];
      arb_add_error(acb_realref(result.values[i][e].raw()),acb_realref(tail.raw()));
      arb_add_error(acb_imagref(result.values[i][e].raw()),acb_realref(tail.raw()));
    }
    if(save)segments.push_back(std::move(result.saved));current=std::move(result.values);
    for(unsigned root=0;root<root_values.size();++root)root_values[root]=continue_polynomial_sqrt(c.square_polynomials[root],at,end,root_values[root]);
    center=next;
    if(charts%20==0)std::cerr<<"Transport chart "<<charts<<", path parameter "<<center<<"\n";
  }
  apply_basis_prefactors(c,B(1),root_values,current,errors);
  if(goal>0) {
    auto tolerance=B::from_strings("1e-"+std::to_string(goal));
    for(unsigned i=0;i<d;++i)for(unsigned e=0;e<=k;++e) {
      auto bound=(magnitude(current[i][e])+B(1))*tolerance;
      auto error=arithmetic_error(current[i][e]);
      if(!arb_le(acb_realref(error.raw()),acb_realref(bound.raw())))throw std::runtime_error("propagated error estimate exceeds AccuracyGoal; increase ExpansionOrder or boundary precision");
    }
  }
  const auto finished=std::chrono::steady_clock::now();
  json::object output{{"schema","DiffExp.TransportResult/v1"},{"parameter","1"},{"values",matrix_json(current,digits)},
    {"errors",matrix_json(errors,digits)},{"charts",charts},{"segments",std::move(segments)},{"omitted_tails_certified",false},
    {"timings",json::object{{"preparation_seconds",std::chrono::duration<double>(prepared-start).count()},
      {"numerical_seconds",std::chrono::duration<double>(finished-prepared).count()},{"total_seconds",std::chrono::duration<double>(finished-start).count()}}}};
  if(!c.basis_prefactors.empty()) {
    output["basis_convention"]="principal_endpoint";
    output["segments_basis_convention"]="continued";
    output["endpoint_roots_continued"]=matrix_json(Boundary{root_values},digits).at(0);
    output["endpoint_roots_principal"]=matrix_json(Boundary{principal_roots_at(c,B(1))},digits).at(0);
  }
  return output;
}
inline int run_file(const std::string& path) {
  std::ifstream file;if(path!="-"){file.open(path);if(!file)throw std::runtime_error("cannot open transport request");}
  auto& in=path=="-"?std::cin:file;std::string text;char block[8192];while(in.read(block,sizeof(block))||in.gcount()){text.append(block,in.gcount());if(text.size()>64*1024*1024)throw std::invalid_argument("transport request exceeds 64 MiB");}
  std::cout<<json::serialize(run(json::parse(text)))<<'\n';return 0;
}
} // namespace diffexp::transport

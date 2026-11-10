#pragma once

#include "diffexp/laurent_transport.hpp"
#include <flint/acb_poly.h>
#include <limits>
#include <map>
#include <memory>

// The ordinary-chart executor for epsilon-regular rational connections
// y' = A(x,eps) y. Every entry N/Q of A becomes a stream w = y_column/Q,
// solved by the finite-lag recurrence Q w = y, and a finite sum
// (n+1) y_{n+1} = sum N w: a Taylor order costs O(degree), not O(order).
// expand() computes the Taylor coefficients of one chart, evaluate() sums
// them at any step chosen afterwards.
// Storage index k means physical epsilon degree lower+k, with a COMMON certified
// finite lower bound supplied by the caller. Thus negative physical degrees are
// supported; missing coefficients above the supplied window are never read.
// Only arithmetic enclosure of the retained Taylor polynomial is claimed here.
namespace diffexp::rational_circuit {
using B = Jet::Ball;
struct Options {
  unsigned max_dimension = 5000, max_degree = 256, max_epsilon_degree = 256;
  std::size_t max_terms = 200000, max_cells = 20000000,
              max_operations = 200000000;
};
struct UnsupportedDomain : std::domain_error {
  using std::domain_error::domain_error;
};
struct InvalidPivot : std::domain_error {
  using std::domain_error::domain_error;
};
struct Scalar { Rational real{0}, imaginary{0}; };
using Polynomial = std::vector<Scalar>;
struct Block { unsigned epsilon, polynomial; };
// A root of the epsilon-zero denominator polynomial, with a multiplicity not
// below its true one. Extra roots or multiplicity only loosen the majorant.
struct Pole { B root; unsigned multiplicity; };
struct Denominator {
  std::vector<Block> blocks;
  unsigned lag = 0;
  bool one = false;
  // Enclosures of every root of the epsilon-zero block. Empty with
  // majorized=false when they could not be isolated.
  std::vector<Pole> poles;
  bool majorized = false;
};
struct Stream {
  unsigned column, denominator;
  // Includes BOTH denominator recurrence and all numerator consumers.
  unsigned lag = 0;
};
struct Entry { unsigned row, stream; std::vector<Block> numerator; };
struct PlanData {
  unsigned dimension = 0, epsilon_high = 0, expected_order = 0;
  std::vector<Polynomial> polynomials;
  std::vector<Denominator> denominators;
  std::vector<Stream> streams;
  std::vector<Entry> entries;
  std::vector<std::vector<unsigned>> row_entries;
  std::size_t visited_terms = 0, coefficient_cells = 0;
  Options options;
};
// Handles own immutable data; a prepared plan remains valid after its original
// compiled handle has been destroyed. No numerical RHS payload lives here.
class Compiled {
 public:
  explicit Compiled(PlanData data)
      : data_(std::make_shared<const PlanData>(std::move(data))) {}
  const PlanData &data() const { return *data_; }
 private:
  std::shared_ptr<const PlanData> data_;
};
struct Statistics {
  // addmul_operations counts the nonzero terms of the pooled dot products.
  std::size_t addmul_operations = 0, scalar_operations = 0, dot_products = 0;
  std::size_t live_cells = 0, auxiliary_cells = 0, coefficient_cells = 0;
  // TaylorTailEstimate of the evaluated chart (see rational_transport.hpp);
  // -infinity when every late term is an exact zero.
  double truncation_log2 = -std::numeric_limits<double>::infinity();
};
namespace detail {
inline bool zero(const Scalar &a) {
  return a.real.is_zero() && a.imaginary.is_zero();
}
inline std::size_t add(std::size_t a, std::size_t b, std::size_t limit,
                       const char *message) {
  if (a > limit || b > limit-a) throw std::length_error(message);
  return a+b;
}
inline std::size_t mul(std::size_t a, std::size_t b, std::size_t limit,
                       const char *message) {
  if (b && a > limit/b) throw std::length_error(message);
  return a*b;
}
inline std::size_t cells(std::size_t a, std::size_t b, std::size_t limit) {
  return mul(a,b,limit,"rational circuit storage budget exhausted");
}
using Support = std::map<unsigned, std::map<unsigned, Scalar>>;
struct Compiler {
  PlanData out;
  std::map<std::string, unsigned> polynomials;
  // I^2=-1 participates in zero tests, equality and interning.
  // Exact has already cancelled factors over its declared rational field.
  // Further polynomial GCD over Q(i) is not attempted: removable chart poles
  // that would require that specialization are rejected by prepare().
  Support support(const Exact &value, const std::vector<Exact::Term> &terms,
                  std::size_t xi, std::optional<std::size_t> ei) {
    Support result;
    for (const auto &term : terms) {
      if (++out.visited_terms > out.options.max_terms)
        throw std::length_error("rational circuit exact term budget exhausted");
      if (term.powers[xi] > out.options.max_degree ||
          (ei && term.powers[*ei] > out.options.max_epsilon_degree))
        throw std::length_error("rational circuit degree budget exhausted");
      unsigned ipower = 0;
      for (std::size_t i=0; i<term.powers.size(); ++i)
        if (i != xi && (!ei || i != *ei) && term.powers[i]) {
          if (value.variables()[i] != "I")
            throw UnsupportedDomain("rational circuit has an unsubstituted parameter");
          ipower = term.powers[i] % 4;
        }
      auto q = ipower >= 2 ? -term.coefficient : term.coefficient;
      auto &a = result[ei ? unsigned(term.powers[*ei]) : 0][unsigned(term.powers[xi])];
      if (ipower % 2) a.imaginary += q;
      else a.real += q;
    }
    for (auto e=result.begin(); e!=result.end();) {
      for (auto x=e->second.begin(); x!=e->second.end();)
        if (zero(x->second)) x=e->second.erase(x); else ++x;
      if (e->second.empty()) e=result.erase(e); else ++e;
    }
    return result;
  }
  std::vector<Block> blocks(const Support &support, unsigned subtract,
                             unsigned offset) {
    std::vector<Block> result;
    for (const auto &[ep, coefficients] : support) {
      const std::uint64_t epsilon = std::uint64_t(ep-subtract)+offset;
      if (epsilon > out.epsilon_high) continue; // Cannot reach this window.
      const auto size = std::size_t(coefficients.rbegin()->first)+1;
      // Budget dense shifted support, before allocation.
      if (size > out.options.max_cells)
        throw std::length_error("rational circuit polynomial storage budget exhausted");
      Polynomial polynomial(size);
      for (const auto &[x,a] : coefficients) polynomial[x]=a;
      std::string key;
      for (const auto &a : polynomial)
        key += a.real.str()+","+a.imaginary.str()+";";
      auto found=polynomials.find(key);
      unsigned id;
      if (found == polynomials.end()) {
        out.coefficient_cells=add(out.coefficient_cells,size,out.options.max_cells,
                                  "rational circuit coefficient storage budget exhausted");
        id=unsigned(out.polynomials.size());
        polynomials.emplace(std::move(key),id);
        out.polynomials.push_back(std::move(polynomial));
      } else id=found->second;
      result.push_back({unsigned(epsilon),id});
    }
    return result;
  }
};
inline unsigned lag(const std::vector<Block> &blocks, const PlanData &data) {
  unsigned result=0;
  for (auto b:blocks)
    result=std::max(result,unsigned(data.polynomials[b.polynomial].size()-1));
  return result;
}
// Integers are converted without rounding; other rationals at bits.
inline void set_exact(arb_t out, const Rational &value, slong bits) {
  if (fmpz_is_one(fmpq_denref(value.raw()))) arb_set_fmpz(out, fmpq_numref(value.raw()));
  else arb_set_fmpq(out, value.raw(), bits);
}
// Roots of an exact Gaussian polynomial d(x), isolated at a fixed precision
// that does not depend on the caller's working precision. A real polynomial
// yields exact multiplicities. Otherwise the roots of the real polynomial
// d*conj(d) are used: its roots are those of d and their conjugates, and a
// conjugate that d certainly does not vanish at is discarded.
inline std::vector<Pole> poles(const Polynomial &d) {
  constexpr slong bits = 128;
  struct Poly { fmpq_poly_t q; Poly() { fmpq_poly_init(q); } ~Poly() { fmpq_poly_clear(q); } } norm;
  bool real = true;
  for (const auto &a : d) real = real && a.imaginary.is_zero();
  if (real) for (std::size_t i=0;i<d.size();++i) fmpq_poly_set_coeff_fmpq(norm.q,i,d[i].real.raw());
  else for (std::size_t i=0;i<d.size();++i) for (std::size_t j=0;j<d.size();++j) {
    Rational sum = d[i].real*d[j].real + d[i].imaginary*d[j].imaginary;
    if (sum.is_zero()) continue;
    fmpq_t old; fmpq_init(old); fmpq_poly_get_coeff_fmpq(old,norm.q,i+j);
    fmpq_add(old,old,sum.raw()); fmpq_poly_set_coeff_fmpq(norm.q,i+j,old); fmpq_clear(old);
  }
  std::vector<Pole> result;
  for (auto &root : isolated_roots(norm.q,bits)) {
    if (!real) {
      B value(0), coefficient;
      for (std::size_t i=d.size();i-->0;) {
        acb_mul(value.raw(),value.raw(),root.value.raw(),bits);
        arb_set_fmpq(acb_realref(coefficient.raw()),d[i].real.raw(),bits);
        arb_set_fmpq(acb_imagref(coefficient.raw()),d[i].imaginary.raw(),bits);
        acb_add(value.raw(),value.raw(),coefficient.raw(),bits);
      }
      if (!acb_contains_zero(value.raw())) continue;
    }
    result.push_back({std::move(root.value),root.multiplicity});
  }
  return result;
}
// Owned array of FLINT magnitudes, all zero initially.
struct Magnitudes {
  std::vector<mag_struct> values;
  explicit Magnitudes(std::size_t n = 0) : values(n) { for (auto &v : values) mag_init(&v); }
  Magnitudes(Magnitudes &&other) noexcept : values(std::move(other.values)) { other.values.clear(); }
  Magnitudes &operator=(Magnitudes &&other) noexcept { std::swap(values,other.values); return *this; }
  Magnitudes(const Magnitudes &other) : values(other.values.size()) {
    for (std::size_t i=0;i<values.size();++i) { mag_init(&values[i]); mag_set(&values[i],&other.values[i]); }
  }
  ~Magnitudes() { for (auto &v : values) mag_clear(&v); }
  mag_struct *operator[](std::size_t i) { return &values[i]; }
  const mag_struct *operator[](std::size_t i) const { return &values[i]; }
  std::size_t size() const { return values.size(); }
};
struct AcbPolynomial {
  acb_poly_t raw;
  AcbPolynomial() { acb_poly_init(raw); }
  ~AcbPolynomial() { acb_poly_clear(raw); }
  AcbPolynomial(const AcbPolynomial &) = delete;
};
} // namespace detail

inline Compiled compile(const std::vector<RationalLineEntry> &entries,
                        unsigned dimension, std::size_t xi,
                        std::optional<std::size_t> ei, unsigned epsilon_high,
                        unsigned expected_order, const Options &options = {}) {
  if (!dimension || dimension>options.max_dimension || !expected_order ||
      expected_order>1000 || epsilon_high>options.max_epsilon_degree ||
      epsilon_high==std::numeric_limits<unsigned>::max() ||
      entries.size()>options.max_terms || !options.max_cells ||
      !options.max_operations || !options.max_terms || (ei && xi==*ei))
    throw std::invalid_argument("rational circuit compilation dimensions or budgets");
  detail::Compiler compiler;
  auto &out=compiler.out;
  out.dimension=dimension; out.epsilon_high=epsilon_high;
  out.expected_order=expected_order; out.options=options;
  out.row_entries.resize(dimension);
  using Key=std::vector<std::pair<unsigned,unsigned>>;
  std::map<Key,unsigned> denominators;
  std::map<std::pair<unsigned,unsigned>,unsigned> streams;
  for (const auto &entry:entries) {
    const auto &a=entry.coefficient;
    if (entry.row>=dimension || entry.column>=dimension || xi>=a.variable_count() ||
        (ei && *ei>=a.variable_count()) || a.variables()[xi]=="I" ||
        (ei && a.variables()[*ei]=="I"))
      throw std::invalid_argument("rational circuit sparse or variable index");
    auto numerator=compiler.support(a,a.numerator_terms(),xi,ei);
    auto denominator=compiler.support(a,a.denominator_terms(),xi,ei);
    if (denominator.empty())
      throw UnsupportedDomain("rational circuit denominator is zero after I specialization");
    if (numerator.empty()) continue;
    const unsigned nv=numerator.begin()->first, dv=denominator.begin()->first;
    const std::uint64_t positive=std::uint64_t(entry.epsilon)+nv;
    if (positive<dv)
      throw UnsupportedDomain("rational circuit requires an epsilon-regular connection");
    const std::uint64_t valuation=positive-dv;
    if (valuation>epsilon_high) continue; // Proven nonnegative valuation.
    // Exact keeps integer polynomials: they stay integers here, so that their
    // shifts to a dyadic chart centre are short exact numbers (see prepare).
    auto qs=compiler.blocks(denominator,dv,0);
    auto ps=compiler.blocks(numerator,nv,unsigned(valuation));
    Key key;
    for (auto b:qs) key.emplace_back(b.epsilon,b.polynomial);
    auto [denpos,newden]=denominators.emplace(std::move(key),out.denominators.size());
    if (newden) {
      const auto &p=out.polynomials[qs[0].polynomial];
      const bool one=qs.size()==1 && qs[0].epsilon==0 && p.size()==1 &&
                     p[0].real==Rational(1) && p[0].imaginary.is_zero();
      const auto qlag=detail::lag(qs,out);
      Denominator denominator{std::move(qs),qlag,one,{},false};
      // The interval recurrence for y/Q wraps: its radii grow at the rate of
      // the root of prod(1+z/|root|)^multiplicity=2, a third or less of the
      // pole distance for a triple pole. With the roots known, the quotient
      // is instead enclosed by its midpoint recurrence plus a majorant.
      try { denominator.poles=detail::poles(p); denominator.majorized=true; }
      catch (const std::runtime_error &) {}
      out.denominators.push_back(std::move(denominator));
    }
    const unsigned did=denpos->second;
    auto [position,inserted]=streams.emplace(std::make_pair(entry.column,did),out.streams.size());
    if (inserted) out.streams.push_back({entry.column,did,out.denominators[did].lag});
    out.streams[position->second].lag=std::max(out.streams[position->second].lag,detail::lag(ps,out));
    out.row_entries[entry.row].push_back(unsigned(out.entries.size()));
    out.entries.push_back({entry.row,position->second,std::move(ps)});
  }
  return Compiled(std::move(out));
}
inline Compiled compile(const std::vector<RationalLineEntry> &entries,
                        unsigned dimension, unsigned expected_order,
                        unsigned epsilon_high, const Options &options = {}) {
  if (entries.empty())
    return compile(entries,dimension,0,std::nullopt,epsilon_high,expected_order,options);
  const auto &names=entries[0].coefficient.variables();
  auto x=std::find(names.begin(),names.end(),"x"), e=std::find(names.begin(),names.end(),"eps");
  if (x==names.end()) throw std::invalid_argument("rational circuit requires x variable");
  return compile(entries,dimension,x-names.begin(),
                 e==names.end()?std::nullopt:std::optional<std::size_t>(e-names.begin()),
                 epsilon_high,expected_order,options);
}
inline Compiled compile(const ExactEpsilonMatrix &matrix,std::size_t xi,std::size_t ei,
                        unsigned expected_order,unsigned epsilon_high,const Options &options={}) {
  std::vector<RationalLineEntry> entries;
  for (unsigned i=0;i<matrix.size();++i) {
    if (matrix[i].size()!=matrix.size()) throw std::invalid_argument("rational circuit square matrix");
    for (unsigned j=0;j<matrix.size();++j)
      if (!matrix[i][j].is_zero()) entries.push_back({i,j,0,matrix[i][j]});
  }
  return compile(entries,unsigned(matrix.size()),xi,ei,epsilon_high,expected_order,options);
}

// Largest number of terms pooled into one dot product for this order/window.
inline std::size_t dot_capacity(const PlanData& plan,unsigned order,unsigned width) {
  if(!order || !width || width>plan.epsilon_high+1)
    throw std::invalid_argument("rational circuit dot schedule order/window");
  const auto limit=plan.options.max_cells;
  std::size_t capacity=0;
  for(const auto& row:plan.row_entries) {
    std::size_t count=0;
    for(auto id:row)for(auto b:plan.entries[id].numerator)if(b.epsilon<width)
      count=detail::add(count,std::min<std::size_t>(order,plan.polynomials[b.polynomial].size()),limit,
                        "rational circuit dot scratch budget exhausted");
    capacity=std::max(capacity,count);
  }
  for(const auto& q:plan.denominators)if(!q.one) {
    std::size_t count=0;
    for(auto b:q.blocks)if(b.epsilon<width)
      count=detail::add(count,std::min<std::size_t>(order,plan.polynomials[b.polynomial].size()),limit,
                        "rational circuit dot scratch budget exhausted");
    capacity=std::max(capacity,count);
  }
  return capacity;
}

// Conservative count of executor arithmetic terms, including zero terms, for
// one chart of this plan. Taylor-shift internals are preparation work and are
// not chart operations. Throws std::length_error beyond plan.options.max_operations,
// so orchestration can test a candidate batch before any numerical preparation.
inline std::size_t operation_upper_bound(const PlanData& plan,unsigned order,unsigned max_width) {
  if(!order || !max_width || max_width>plan.epsilon_high+1)
    throw std::invalid_argument("rational circuit operation estimate order/window");
  std::size_t total=0;
  auto spend=[&](std::size_t count) {
    total=detail::add(total,count,plan.options.max_operations,
                      "rational circuit predicted operation budget exhausted");
  };
  spend(detail::mul(detail::mul(order,plan.dimension,plan.options.max_operations,"rational circuit operation budget"),
                    max_width,plan.options.max_operations,"rational circuit operation budget"));
  spend(detail::mul(detail::mul(2*(order+1),plan.dimension,plan.options.max_operations,"rational circuit operation budget"),
                    max_width,plan.options.max_operations,"rational circuit operation budget"));
  for (unsigned n=0;n<order;++n) {
    for (const auto &stream:plan.streams) {
      const auto &q=plan.denominators[stream.denominator];
      if (q.one) continue;
      spend(max_width); // cached inverse pivot multiplication
      for (auto b:q.blocks) if (b.epsilon<max_width) {
        std::size_t length=std::min<std::size_t>(n+1,plan.polynomials[b.polynomial].size());
        if (b.epsilon==0) --length;
        spend(detail::mul(length,max_width-b.epsilon,plan.options.max_operations,"rational circuit operation budget"));
      }
    }
    for (const auto &entry:plan.entries)
      for (auto b:entry.numerator) if (b.epsilon<max_width)
        spend(detail::mul(std::min<std::size_t>(n+1,plan.polynomials[b.polynomial].size()),
                          max_width-b.epsilon,plan.options.max_operations,"rational circuit operation budget"));
  }
  return total;
}

struct PreparedData {
  Compiled compiled;
  B center;
  slong precision;
  unsigned order, max_width;
  std::vector<std::vector<B>> polynomials;
  std::vector<B> inverse_pivots;
  // Per denominator, an upper bound of 1/|root-center| for every pole,
  // repeated by multiplicity: the factors of the quotient majorant.
  std::vector<detail::Magnitudes> inverse_distances;
  std::vector<unsigned> ring_lengths;
  std::size_t coefficient_cells=0, auxiliary_cells=0, workspace_cells=0,
              live_cells=0, preparation_cells=0, operation_upper_bound=0, dot_capacity=0;
};
class Prepared {
 public:
  explicit Prepared(PreparedData data)
      : data_(std::make_shared<const PreparedData>(std::move(data))) {}
  const PreparedData &data() const { return *data_; }
 private:
  std::shared_ptr<const PreparedData> data_;
};
// max_width is an epsilon coefficient count, not a physical epsilon order.
// The cache identity is (immutable plan, center ball, order, width, precision).
inline Prepared prepare(const Compiled &compiled,const B &center,unsigned order,
                        unsigned max_width=0) {
  const auto &plan=compiled.data();
  if (!max_width) max_width=plan.epsilon_high+1;
  if (!order || order>1000 || max_width>plan.epsilon_high+1 || !center.is_finite())
    throw std::invalid_argument("rational circuit preparation order, width or center");
  PreparedData out{compiled,center,B::precision(),order,max_width,{},{},{}};
  const auto limit=plan.options.max_cells;
  out.coefficient_cells=plan.coefficient_cells;
  for (const auto &stream:plan.streams) {
    const bool alias=plan.denominators[stream.denominator].one;
    const unsigned length=alias?0:std::min(order,stream.lag+1);
    out.ring_lengths.push_back(length);
    out.auxiliary_cells=detail::add(out.auxiliary_cells,detail::cells(length,max_width,limit),limit,
                                  "rational circuit auxiliary storage budget exhausted");
  }
  out.dot_capacity=dot_capacity(plan,order,max_width);
  const auto ycells=detail::cells(detail::cells(order+1,plan.dimension,limit),max_width,limit);
  const auto resultcells=detail::cells(plan.dimension,max_width,limit);
  out.workspace_cells=detail::add(ycells,out.auxiliary_cells,limit,"rational circuit workspace budget exhausted");
  out.workspace_cells=detail::add(out.workspace_cells,resultcells,limit,"rational circuit result storage budget exhausted");
  // Borrowed acb_struct views contain no owned limbs, but their storage is
  // nevertheless charged as two logical scalar slots per possible dot term.
  out.workspace_cells=detail::add(out.workspace_cells,
    detail::add(detail::cells(2,out.dot_capacity,limit),1,limit,"rational circuit dot scratch budget exhausted"),
    limit,"rational circuit dot workspace budget exhausted");
  out.live_cells=detail::add(out.workspace_cells,out.coefficient_cells,limit,"rational circuit live storage budget exhausted");
  out.live_cells=detail::add(out.live_cells,plan.denominators.size()+2,limit,"rational circuit pivot storage budget exhausted");
  // Count owned ball cells, including the two visible FLINT polynomial shift
  // buffers. FLINT internal workspaces and arbitrary-precision limb bytes are
  // not inferred from these logical cell counts.
  std::size_t largest_polynomial=0;
  for (const auto &polynomial:plan.polynomials)
    largest_polynomial=std::max(largest_polynomial,polynomial.size());
  out.preparation_cells=detail::add(out.coefficient_cells,
      detail::cells(2,largest_polynomial,limit),limit,"rational circuit preparation storage budget exhausted");
  out.preparation_cells=detail::add(out.preparation_cells,plan.denominators.size()+3,limit,
      "rational circuit preparation pivot storage budget exhausted");
  out.operation_upper_bound=operation_upper_bound(plan,order,max_width);
  out.polynomials.reserve(plan.polynomials.size());
  for (const auto &p:plan.polynomials) {
    // Exact integer coefficients shifted to a dyadic centre are exact dyadic
    // numbers. Evaluate the shift without rounding and round only what does
    // not fit the working precision: short exact coefficients make every
    // recurrence product cheaper and free of coefficient radii.
    std::vector<B> coefficients(p.size());slong exact=64;
    for (unsigned n=0;n<p.size();++n) {
      detail::set_exact(acb_realref(coefficients[n].raw()),p[n].real,out.precision);
      detail::set_exact(acb_imagref(coefficients[n].raw()),p[n].imaginary,out.precision);
      exact=std::max<slong>(exact,acb_bits(coefficients[n].raw())+64);
    }
    exact+=static_cast<slong>(p.size())*(acb_bits(center.raw())+2);
    // Horner shift p(center+t): every intermediate fits in `exact` bits.
    for (std::size_t i=p.size();i-->1;)
      for (std::size_t j=i-1;j+1<p.size();++j)
        acb_addmul(coefficients[j].raw(),center.raw(),coefficients[j+1].raw(),exact);
    for (auto &coefficient:coefficients)
      acb_set_round(coefficient.raw(),coefficient.raw(),out.precision);
    out.polynomials.push_back(std::move(coefficients));
  }
  out.inverse_pivots.reserve(plan.denominators.size());
  for (const auto &denominator:plan.denominators) {
    B pivot(0);
    for (auto b:denominator.blocks) if (b.epsilon==0)
      pivot+=out.polynomials[b.polynomial][0];
    if (!pivot.is_finite() || pivot.contains_zero())
      throw InvalidPivot("rational circuit Q(center,epsilon=0) does not exclude zero; nonunit or unresolved removable pole");
    B inverse;
    acb_inv(inverse.raw(),pivot.raw(),out.precision);
    if (!inverse.is_finite()) throw InvalidPivot("rational circuit inverse pivot is not finite");
    out.inverse_pivots.push_back(std::move(inverse));
    std::size_t levels=0;
    for (const auto &pole:denominator.poles) levels+=pole.multiplicity;
    detail::Magnitudes inverses(levels);
    levels=0;
    for (const auto &pole:denominator.poles) {
      B distance;
      acb_sub(distance.raw(),pole.root.raw(),center.raw(),128);
      detail::Magnitudes lower(1);
      acb_get_mag_lower(lower[0],distance.raw());
      if (mag_is_zero(lower[0])) throw InvalidPivot("rational circuit center is not separated from a denominator root");
      for (unsigned copy=0;copy<pole.multiplicity;++copy) mag_inv(inverses[levels++],lower[0]);
    }
    out.inverse_distances.push_back(std::move(inverses));
  }
  return Prepared(std::move(out));
}

// Taylor coefficients of one chart for one right-hand side. They do not
// depend on the step, which can therefore be chosen (and shortened) after the
// recurrence from the retained tail, without recomputing anything.
struct Series {
  unsigned dimension=0,width=0,order=0;
  slong precision=0;
  // Structurally reachable rows, in storage order; every other row is zero.
  std::vector<unsigned> rows;
  // y[(n*rows.size()+r)*width+k] is Taylor order n, row rows[r], epsilon index k.
  std::vector<B> y;
  const B& at(unsigned n,std::size_t r,unsigned k)const {return y[(std::size_t(n)*rows.size()+r)*width+k];}
};
inline Series expand(const Prepared &prepared,const Boundary &boundary,Statistics *statistics=nullptr) {
  const auto &p=prepared.data();
  const auto &plan=p.compiled.data();
  if (B::precision()!=p.precision)
    throw std::invalid_argument("rational circuit prepared precision mismatch");
  if (boundary.size()!=plan.dimension || boundary.empty() || boundary[0].empty() ||
      boundary[0].size()>p.max_width)
    throw std::invalid_argument("rational circuit boundary shape");
  const unsigned width=unsigned(boundary[0].size()), d=plan.dimension, order=p.order;
  for (const auto &row:boundary) {
    if (row.size()!=width) throw std::invalid_argument("rational circuit ragged boundary");
    for (const auto &value:row)
      if (!value.is_finite()) throw std::invalid_argument("rational circuit nonfinite boundary");
  }
  // Exact structural reachability is computed separately for each RHS. A zero
  // ball proves a zero input; merely small balls never suppress dependencies.
  std::vector<std::vector<unsigned>> dependents(d);
  for (const auto &entry:plan.entries)
    if (!entry.numerator.empty() && entry.numerator.front().epsilon<width)
      dependents[plan.streams[entry.stream].column].push_back(entry.row);
  std::vector<bool> active(d,false), active_stream(plan.streams.size(),false);
  Series series{d,width,order,p.precision,{},{}};
  auto &queue=series.rows;
  for (unsigned i=0;i<d;++i)
    for (const auto &value:boundary[i]) if (!value.is_zero()) {
      active[i]=true; queue.push_back(i); break;
    }
  for (std::size_t next=0;next<queue.size();++next)
    for (auto row:dependents[queue[next]]) if (!active[row]) {
      active[row]=true; queue.push_back(row);
    }
  std::vector<unsigned> row_index(d);
  for (unsigned i=0;i<queue.size();++i) row_index[queue[i]]=i;
  for (const auto &entry:plan.entries)
    if (active[plan.streams[entry.stream].column] &&
        !entry.numerator.empty() && entry.numerator.front().epsilon<width)
      active_stream[entry.stream]=true;
  Statistics stats;
  stats.coefficient_cells=p.coefficient_cells;
  auto &y=series.y;
  y.resize(std::size_t(order+1)*queue.size()*width);
  std::vector<std::size_t> offsets(plan.streams.size());
  for (unsigned id=0;id<plan.streams.size();++id) {
    offsets[id]=stats.auxiliary_cells;
    if (active_stream[id]) stats.auxiliary_cells+=std::size_t(p.ring_lengths[id])*width;
  }
  std::vector<B> w(stats.auxiliary_cells);
  stats.live_cells=y.size()+w.size()+std::size_t(d)*width+p.coefficient_cells+p.inverse_pivots.size()+2+
    2*p.dot_capacity+1;
  // These are non-owning copies of read-only FLINT structs. Never clear them:
  // their limb pointers belong to the prepared coefficients and Y/W arrays.
  // acb_dot reads these views and writes only the distinct accumulator below.
  std::vector<acb_struct> dot_left(p.dot_capacity),dot_right(p.dot_capacity);
  B dot_result;
  std::size_t dot_count=0;
  // midpoint=true enters only the exact midpoint of b (its radius is
  // accounted for separately by the caller).
  auto collect=[&](const B& a,const B& b,bool midpoint=false) {
    if(dot_count>=p.dot_capacity)throw std::logic_error("rational circuit dot schedule exceeded");
    dot_left[dot_count]=*a.raw();dot_right[dot_count]=*b.raw();
    if(midpoint) {
      auto* view=&dot_right[dot_count];
      for(auto* radius:{arb_radref(acb_realref(view)),arb_radref(acb_imagref(view))}){MAG_EXP(radius)=0;MAG_MAN(radius)=0;}
    }
    ++dot_count;
  };
  auto flush=[&](B& value,bool subtract) {
    if(!dot_count)return;
    acb_dot(dot_result.raw(),value.raw(),subtract,dot_left.data(),1,dot_right.data(),1,
            static_cast<slong>(dot_count),p.precision);
    acb_swap(value.raw(),dot_result.raw());++stats.dot_products;stats.addmul_operations+=dot_count;dot_count=0;
  };
  auto yat=[&](unsigned n,unsigned i,unsigned k)->B& {
    return y[(std::size_t(n)*queue.size()+row_index[i])*width+k];
  };
  // slot[id] is the ring row of Taylor order n; earlier orders precede it cyclically.
  std::vector<unsigned> slot(plan.streams.size(),0);
  auto wat=[&](unsigned id,unsigned n,unsigned back,unsigned k)->B& {
    if (!p.ring_lengths[id]) return yat(n-back,plan.streams[id].column,k);
    const unsigned row=slot[id]>=back?slot[id]-back:slot[id]+p.ring_lengths[id]-back;
    return w[offsets[id]+std::size_t(row)*width+k];
  };
  for (auto i:queue)
    for (unsigned k=0;k<width;++k) yat(0,i,k)=boundary[i][k];
  // Majorant state per active stream: chain[(level*width)+k] is the bound of
  // the previous Taylor order after `level+1` of the pole factors.
  std::vector<detail::Magnitudes> chains(plan.streams.size());
  for (unsigned id=0;id<plan.streams.size();++id)
    if (active_stream[id] && p.ring_lengths[id] && plan.denominators[plan.streams[id].denominator].majorized)
      chains[id]=detail::Magnitudes(p.inverse_distances[plan.streams[id].denominator].size()*width);
  detail::Magnitudes scratch(2);
  for (unsigned n=0;n<order;++n) {
    // QW=Y is triangular in (Taylor degree, epsilon storage index).
    for (unsigned id=0;id<plan.streams.size();++id) {
      if (!active_stream[id] || !p.ring_lengths[id]) continue;
      slot[id]=n%p.ring_lengths[id];
      const auto &stream=plan.streams[id];
      const auto &denominator=plan.denominators[stream.denominator];
      const auto &inverse=p.inverse_distances[stream.denominator];
      for (unsigned k=0;k<width;++k) {
        auto &value=wat(id,n,0,k);
        value=yat(n,stream.column,k); // overwrite the expired ring row
        for (auto b:denominator.blocks) if (b.epsilon<=k) {
          const auto &q=p.polynomials[b.polynomial];
          for (unsigned a=0;a<q.size() && a<=n;++a) {
            if ((!a && !b.epsilon) || q[a].is_zero()) continue;
            const auto &previous=wat(id,n,a,k-b.epsilon);
            if (!previous.is_zero()) collect(q[a],previous,denominator.majorized && !b.epsilon);
          }
        }
        flush(value,true);
        ++stats.scalar_operations;
        acb_mul(value.raw(),value.raw(),p.inverse_pivots[stream.denominator].raw(),p.precision);
        if (!denominator.majorized) continue;
        // With Q0 the epsilon-zero block and m the midpoints of earlier
        // quotients, `value` encloses (z_n - sum_{a>0} q_a m_{n-a})/q_0 for
        // every admissible right-hand side z, so its radius bounds the defect
        // d_n of Q0*m = z - q_0*d. Hence w-m = (q_0*d)/Q0, and the Taylor
        // coefficients of q_0/Q0(center+t) = prod(1-t/root_i)^-1 are bounded
        // by those of prod(1-t/|root_i|)^-1: a chain of positive first-order
        // recurrences, free of the cancellation an interval recurrence loses.
        auto *bound=scratch[0];
        auto *real=arb_radref(acb_realref(value.raw())), *imaginary=arb_radref(acb_imagref(value.raw()));
        mag_hypot(bound,real,imaginary);
        for (std::size_t level=0;level<inverse.size();++level) {
          auto *previous=chains[id][level*width+k];
          mag_addmul(bound,previous,inverse[level]);
          mag_set(previous,bound);
        }
        mag_set(real,bound); mag_set(imaginary,bound);
      }
    }
    for(auto row:queue)for(unsigned k=0;k<width;++k) {
      for(auto id:plan.row_entries[row]) {
        const auto& entry=plan.entries[id];
        if(!active_stream[entry.stream])continue;
        for(auto b:entry.numerator)if(b.epsilon<=k) {
          const auto& numerator=p.polynomials[b.polynomial];
          for(unsigned a=0;a<numerator.size() && a<=n;++a) {
            if(numerator[a].is_zero())continue;
            const auto& source=wat(entry.stream,n,a,k-b.epsilon);
            if(!source.is_zero())collect(numerator[a],source);
          }
        }
      }
      auto& value=yat(n+1,row,k);
      flush(value,false);
      ++stats.scalar_operations;
      acb_div_ui(value.raw(),value.raw(),n+1,p.precision);
    }
  }
  if (statistics) *statistics=stats;
  return series;
}
// Horner value of the retained polynomial at step; truncation_log2 receives
// its TaylorTailEstimate (see rational_transport.hpp) when requested.
inline Boundary evaluate(const Series &series,const B &step,double *truncation_log2=nullptr) {
  if (B::precision()!=series.precision || !step.is_finite())
    throw std::invalid_argument("rational circuit series precision or step");
  Boundary result;
  result.reserve(series.dimension);
  for (unsigned i=0;i<series.dimension;++i) result.emplace_back(series.width);
  for (std::size_t r=0;r<series.rows.size();++r)
    for (unsigned k=0;k<series.width;++k) {
      auto &value=result[series.rows[r]][k];
      for (unsigned n=series.order+1;n-->0;) {
        acb_mul(value.raw(),value.raw(),step.raw(),series.precision);
        acb_add(value.raw(),value.raw(),series.at(n,r,k).raw(),series.precision);
      }
    }
  if (truncation_log2) {
    TaylorTailEstimate tail(step,series.order);
    for (std::size_t r=0;r<series.rows.size();++r)
      for (unsigned k=0;k<series.width;++k)
        tail.component(result[series.rows[r]][k].raw(),[&](unsigned n){return series.at(n,r,k).raw();});
    *truncation_log2=tail.log2_relative();
  }
  return result;
}
inline Boundary chart(const Prepared &prepared,const Boundary &boundary,const B &step,
                      Statistics *statistics=nullptr) {
  Statistics stats;
  auto result=evaluate(expand(prepared,boundary,&stats),step,&stats.truncation_log2);
  if (statistics) *statistics=stats;
  return result;
}
inline Boundary chart(const Compiled &compiled,const Boundary &boundary,const B &center,
                      const B &step,unsigned order,Statistics *statistics=nullptr) {
  if (boundary.empty() || boundary[0].empty())
    throw std::invalid_argument("rational circuit empty boundary");
  return chart(prepare(compiled,center,order,unsigned(boundary[0].size())),boundary,step,statistics);
}
} // namespace diffexp::rational_circuit

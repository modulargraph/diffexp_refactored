#pragma once

#include "diffexp/laurent_transport.hpp"
#include <flint/acb_poly.h>
#include <limits>
#include <map>
#include <memory>

// A compact ordinary-chart executor for epsilon-regular rational connections.
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
  // Pool finite convolution terms in FLINT dot products. False retains the
  // scalar accumulation oracle at the same precision, order and support.
  bool grouped_dot = true;
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
struct Denominator {
  std::vector<Block> blocks;
  unsigned lag = 0;
  bool one = false;
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
  // addmul_operations counts nonzero multiplication-accumulation terms,
  // including terms executed together by a FLINT dot product.
  std::size_t addmul_operations = 0, scalar_operations = 0, dot_products = 0;
  std::size_t live_cells = 0, auxiliary_cells = 0, coefficient_cells = 0;
};
namespace detail {
inline bool zero(const Scalar &a) {
  return a.real.is_zero() && a.imaginary.is_zero();
}
inline Scalar divide(const Scalar &a, const Scalar &b) {
  const auto norm = b.real*b.real + b.imaginary*b.imaginary;
  if (norm.is_zero()) throw UnsupportedDomain("zero Gaussian denominator");
  return {(a.real*b.real+a.imaginary*b.imaginary)/norm,
          (a.imaginary*b.real-a.real*b.imaginary)/norm};
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
  // I^2=-1 participates in zero tests, equality, normalization and interning.
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
                             unsigned offset, const Scalar &normalizer) {
    std::vector<Block> result;
    for (const auto &[ep, coefficients] : support) {
      const std::uint64_t epsilon = std::uint64_t(ep-subtract)+offset;
      if (epsilon > out.epsilon_high) continue; // Cannot reach this window.
      const auto size = std::size_t(coefficients.rbegin()->first)+1;
      // Budget dense shifted support, before allocation.
      if (size > out.options.max_cells)
        throw std::length_error("rational circuit polynomial storage budget exhausted");
      Polynomial polynomial(size);
      for (const auto &[x,a] : coefficients) polynomial[x]=divide(a,normalizer);
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
inline B ball(const Scalar &value) {
  B result;
  arb_set_fmpq(acb_realref(result.raw()), value.real.raw(), B::precision());
  arb_set_fmpq(acb_imagref(result.raw()), value.imaginary.raw(), B::precision());
  return result;
}
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
    const auto normalizer=denominator.begin()->second.begin()->second;
    auto qs=compiler.blocks(denominator,dv,0,normalizer);
    auto ps=compiler.blocks(numerator,nv,unsigned(valuation),normalizer);
    Key key;
    for (auto b:qs) key.emplace_back(b.epsilon,b.polynomial);
    auto [denpos,newden]=denominators.emplace(std::move(key),out.denominators.size());
    if (newden) {
      const auto &p=out.polynomials[qs[0].polynomial];
      const bool one=qs.size()==1 && qs[0].epsilon==0 && p.size()==1 &&
                     p[0].real==Rational(1) && p[0].imaginary.is_zero();
      const auto qlag=detail::lag(qs,out);
      out.denominators.push_back({std::move(qs),qlag,one});
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

// Maximum borrowed term slots for the selected order/window. This is exposed
// so orchestration can account for complete per-worker scratch before batching.
inline std::size_t dot_capacity(const PlanData& plan,unsigned order,unsigned width) {
  if(!order || !width || width>plan.epsilon_high+1)
    throw std::invalid_argument("rational circuit dot schedule order/window");
  const auto limit=plan.options.max_cells;
  std::size_t capacity=0;
  if(plan.options.grouped_dot) {
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
  }
  return capacity;
}

struct PreparedData {
  Compiled compiled;
  B center;
  slong precision;
  unsigned order, max_width;
  std::vector<std::vector<B>> polynomials;
  std::vector<B> inverse_pivots;
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
  if(plan.options.grouped_dot)
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
  // Conservative count of executor arithmetic terms, including zero terms.
  // Taylor-shift internals are preparation work and are not chart operations.
  auto spend=[&](std::size_t count) {
    out.operation_upper_bound=detail::add(out.operation_upper_bound,count,plan.options.max_operations,
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
  out.polynomials.reserve(plan.polynomials.size());
  for (const auto &p:plan.polynomials) {
    detail::AcbPolynomial input,shifted;
    for (unsigned n=0;n<p.size();++n) {
      const auto a=detail::ball(p[n]);
      acb_poly_set_coeff_acb(input.raw,n,a.raw());
    }
    acb_poly_taylor_shift(shifted.raw,input.raw,center.raw(),out.precision);
    std::vector<B> coefficients(p.size());
    for (unsigned n=0;n<p.size();++n)
      acb_poly_get_coeff_acb(coefficients[n].raw(),shifted.raw,n);
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
  }
  return Prepared(std::move(out));
}

inline Boundary chart(const Prepared &prepared,const Boundary &boundary,const B &step,
                      Statistics *statistics=nullptr) {
  const auto &p=prepared.data();
  const auto &plan=p.compiled.data();
  if (B::precision()!=p.precision)
    throw std::invalid_argument("rational circuit prepared precision mismatch");
  if (boundary.size()!=plan.dimension || boundary.empty() || boundary[0].empty() ||
      boundary[0].size()>p.max_width || !step.is_finite())
    throw std::invalid_argument("rational circuit boundary shape or step");
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
  std::vector<unsigned> queue;
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
  std::vector<B> y(std::size_t(order+1)*queue.size()*width);
  std::vector<std::size_t> offsets(plan.streams.size());
  for (unsigned id=0;id<plan.streams.size();++id) {
    offsets[id]=stats.auxiliary_cells;
    if (active_stream[id]) stats.auxiliary_cells+=std::size_t(p.ring_lengths[id])*width;
  }
  std::vector<B> w(stats.auxiliary_cells);
  stats.live_cells=y.size()+w.size()+std::size_t(d)*width+p.coefficient_cells+p.inverse_pivots.size()+2+
    (plan.options.grouped_dot?2*p.dot_capacity+1:0);
  // These are non-owning copies of read-only FLINT structs. Never clear them:
  // their limb pointers belong to the prepared coefficients and Y/W arrays.
  // acb_dot reads these views and writes only the distinct accumulator below.
  std::vector<acb_struct> dot_left(p.dot_capacity),dot_right(p.dot_capacity);
  std::optional<B> dot_result;
  if(plan.options.grouped_dot)dot_result.emplace();
  std::size_t dot_count=0;
  // Report measured counts even when an arithmetic budget interrupts execution.
  struct Publish { Statistics *target; Statistics &value;
    ~Publish() { if (target) *target=value; }
  } publish{statistics,stats};
  auto spend=[&](bool addmul) {
    auto &count=addmul?stats.addmul_operations:stats.scalar_operations;
    ++count;
    if (stats.addmul_operations+stats.scalar_operations>plan.options.max_operations)
      throw std::length_error("rational circuit actual operation budget exhausted");
  };
  auto collect=[&](const B& a,const B& b) {
    if(dot_count>=p.dot_capacity)throw std::logic_error("rational circuit dot schedule exceeded");
    dot_left[dot_count]=*a.raw();dot_right[dot_count]=*b.raw();++dot_count;
    spend(true);
  };
  auto flush=[&](B& value,bool subtract) {
    if(!dot_count)return;
    acb_dot(dot_result->raw(),value.raw(),subtract,dot_left.data(),1,dot_right.data(),1,
            static_cast<slong>(dot_count),p.precision);
    acb_swap(value.raw(),dot_result->raw());++stats.dot_products;dot_count=0;
  };
  auto yat=[&](unsigned n,unsigned i,unsigned k)->B& {
    return y[(std::size_t(n)*queue.size()+row_index[i])*width+k];
  };
  auto wat=[&](unsigned id,unsigned n,unsigned k)->B& {
    const auto &stream=plan.streams[id];
    if (!p.ring_lengths[id]) return yat(n,stream.column,k);
    return w[offsets[id]+std::size_t(n%p.ring_lengths[id])*width+k];
  };
  for (auto i:queue)
    for (unsigned k=0;k<width;++k) yat(0,i,k)=boundary[i][k];
  for (unsigned n=0;n<order;++n) {
    // QW=Y is triangular in (Taylor degree, epsilon storage index).
    for (unsigned id=0;id<plan.streams.size();++id) {
      if (!active_stream[id] || !p.ring_lengths[id]) continue;
      const auto &stream=plan.streams[id];
      const auto &denominator=plan.denominators[stream.denominator];
      for (unsigned k=0;k<width;++k) {
        auto &value=wat(id,n,k);
        value=yat(n,stream.column,k); // overwrite the expired ring row
        for (auto b:denominator.blocks) if (b.epsilon<=k) {
          const auto &q=p.polynomials[b.polynomial];
          for (unsigned a=0;a<q.size() && a<=n;++a) {
            if ((!a && !b.epsilon) || q[a].is_zero()) continue;
            const auto &previous=wat(id,n-a,k-b.epsilon);
            if (previous.is_zero()) continue;
            if(plan.options.grouped_dot)collect(q[a],previous);
            else {spend(true);acb_submul(value.raw(),q[a].raw(),previous.raw(),p.precision);}
          }
        }
        if(plan.options.grouped_dot)flush(value,true);
        spend(false);
        acb_mul(value.raw(),value.raw(),p.inverse_pivots[stream.denominator].raw(),p.precision);
      }
    }
    if(plan.options.grouped_dot) {
      for(auto row:queue)for(unsigned k=0;k<width;++k) {
        for(auto id:plan.row_entries[row]) {
          const auto& entry=plan.entries[id];
          if(!active_stream[entry.stream])continue;
          for(auto b:entry.numerator)if(b.epsilon<=k) {
            const auto& numerator=p.polynomials[b.polynomial];
            for(unsigned a=0;a<numerator.size() && a<=n;++a) {
              if(numerator[a].is_zero())continue;
              const auto& source=wat(entry.stream,n-a,k-b.epsilon);
              if(!source.is_zero())collect(numerator[a],source);
            }
          }
        }
        flush(yat(n+1,row,k),false);
      }
    } else for (const auto &entry:plan.entries) if (active_stream[entry.stream])
      for (auto b:entry.numerator) if (b.epsilon<width) {
        const auto &numerator=p.polynomials[b.polynomial];
        for (unsigned a=0;a<numerator.size() && a<=n;++a) {
          if (numerator[a].is_zero()) continue;
          for (unsigned k=b.epsilon;k<width;++k) {
            const auto &source=wat(entry.stream,n-a,k-b.epsilon);
            if (source.is_zero()) continue;
            spend(true);
            acb_addmul(yat(n+1,entry.row,k).raw(),numerator[a].raw(),source.raw(),p.precision);
          }
        }
      }
    for (auto i:queue)
      for (unsigned k=0;k<width;++k) {
        spend(false);
        acb_div_ui(yat(n+1,i,k).raw(),yat(n+1,i,k).raw(),n+1,p.precision);
      }
  }
  Boundary result;
  result.reserve(d);
  for (unsigned i=0;i<d;++i) result.emplace_back(width);
  for (auto i:queue)
    for (unsigned k=0;k<width;++k)
      for (unsigned n=order+1;n-->0;) {
        spend(false); acb_mul(result[i][k].raw(),result[i][k].raw(),step.raw(),p.precision);
        spend(false); acb_add(result[i][k].raw(),result[i][k].raw(),yat(n,i,k).raw(),p.precision);
      }
  return result;
}
inline Boundary chart(const Compiled &compiled,const Boundary &boundary,const B &center,
                      const B &step,unsigned order,Statistics *statistics=nullptr) {
  if (boundary.empty() || boundary[0].empty())
    throw std::invalid_argument("rational circuit empty boundary");
  return chart(prepare(compiled,center,order,unsigned(boundary[0].size())),boundary,step,statistics);
}
} // namespace diffexp::rational_circuit

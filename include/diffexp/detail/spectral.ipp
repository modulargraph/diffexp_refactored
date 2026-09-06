#pragma once
#include <array>
#include <charconv>
#include <cmath>
#include <numeric>
#include <optional>

// Ordinary sequential-epsilon Chebyshev collocation, following the CHESS
// method of Yuanche Liu and Yang Zhang. Refinement below estimates the omitted
// spectral tail; it does not turn finite ball arithmetic into a tail proof.
namespace diffexp::transport {
struct SpectralOptions {
  unsigned accuracy_goal=30;
  unsigned max_nodes=128;
  slong matrix_bits=0;
  double seconds_budget=30;
  std::uint64_t work_budget=100000000;
  std::size_t cell_budget=2000000;
  unsigned guard_digits=1;
};
struct SpectralDiagnostics {
  std::vector<unsigned> nodes_tried;
  unsigned selected_nodes=0;
  unsigned sampled_nodes=0;
  unsigned absolute_stability_components=0;
  double preparation_seconds=0;
  double numerical_seconds=0;
  std::string rejection_reason;
};
struct SpectralResult {
  Boundary values; // Incoming boxes and finite arithmetic only.
  Boundary errors; // Estimated omitted spectral tails, not injected here.
  SpectralDiagnostics diagnostics;
  std::vector<B> endpoint_roots;
};
namespace spectral_detail {
// A cost/sampling guard, not a theorem about the solution's convergence rate.
inline double node_forecast(const Compiled& c,unsigned goal) {
  double log_rho=std::numeric_limits<double>::infinity();
  for(const auto& singularity:c.singularities) {
    const std::complex<double> z(2*arf_get_d(arb_midref(acb_realref(singularity.raw())),ARF_RND_NEAR)-1,
        2*arf_get_d(arb_midref(acb_imagref(singularity.raw())),ARF_RND_NEAR));
    const auto root=std::sqrt(z*z-std::complex<double>(1));
    const auto rho=std::max(std::abs(z+root),std::abs(z-root));
    log_rho=std::min(log_rho,std::log(rho));
  }
  return log_rho>0?(goal+4)*std::log(10.)/log_rho:std::numeric_limits<double>::infinity();
}
struct Rejected:std::runtime_error {using std::runtime_error::runtime_error;};
struct PrecisionScope {
  slong previous;
  explicit PrecisionScope(slong bits):previous(B::precision()){B::set_precision(bits);}
  ~PrecisionScope(){B::set_precision(previous);}
};
struct Matrix {
  acb_mat_t value;
  Matrix(unsigned rows,unsigned columns){acb_mat_init(value,rows,columns);}
  ~Matrix(){acb_mat_clear(value);}
  Matrix(const Matrix&)=delete;
  Matrix& operator=(const Matrix&)=delete;
};
inline B upper(const B& nonnegative) {
  B result;mag_t m;mag_init(m);arb_get_mag(m,acb_realref(nonnegative.raw()));
  arf_set_mag(arb_midref(acb_realref(result.raw())),m);mag_clear(m);return result;
}
inline B difference(const B& a,const B& b) {
  B x,y;acb_get_mid(x.raw(),a.raw());acb_get_mid(y.raw(),b.raw());
  return upper(magnitude(x-y));
}
inline B midpoint_scale(const B& value) {
  B midpoint;acb_get_mid(midpoint.raw(),value.raw());return B(1)+magnitude(midpoint);
}
inline bool le(const B& a,const B& b){return arb_le(acb_realref(a.raw()),acb_realref(b.raw()));}
inline double ratio_upper(const B& numerator,const B& denominator) {
  B ratio=numerator/denominator;arf_t a;arf_init(a);
  arb_get_ubound_arf(a,acb_realref(ratio.raw()),B::precision());
  double result=arf_get_d(a,ARF_RND_CEIL);arf_clear(a);return result;
}
// If errors behave as C*q^N, two differences with node increments a,b obey
// r=q^a*(1-q^b)/(1-q^a). Solve for q rather than assuming equal increments.
inline std::optional<double> tail_factor(double ratio,unsigned a,unsigned b) {
  if(!(ratio>=0) || !std::isfinite(ratio) || ratio>=0.75)return std::nullopt;
  if(ratio==0)return 0.0;
  double low=0,high=1;
  for(unsigned iteration=0;iteration<64;++iteration) {
    const double q=(low+high)/2,logq=std::log(q);
    const double model=std::exp(a*logq)*std::expm1(b*logq)/std::expm1(a*logq);
    if(model<ratio)low=q;else high=q;
  }
  // Keep positive underflows conservative even at thousands of requested digits.
  const double remaining=std::max(std::pow(high,b),std::numeric_limits<double>::min());
  if(!(remaining<0.75))return std::nullopt;
  return 16.0*remaining/(1.0-remaining);
}
// The same nonzero-image-disk proof as geometry's continuation, with checks
// inside recursive subdivisions so a difficult root cannot consume the budget.
template<class Check>
B continue_root(const Jet& polynomial,const B& left,const B& right,const B& root,
    Check& check,unsigned depth=0) {
  check();if(depth>24)throw Rejected("spectral root continuation subdivision limit");
  B interval;arb_union(acb_realref(interval.raw()),acb_realref(left.raw()),acb_realref(right.raw()),polynomial.bits());
  auto start=polynomial.evaluate_polynomial(left);
  if(start.contains_zero())throw Rejected("spectral root continuation meets a branch point");
  auto deviation=polynomial.evaluate_polynomial(interval)/start-B(1);
  arf_t bound;arf_init(bound);acb_get_abs_ubound_arf(bound,deviation.raw(),polynomial.bits());
  bool safe=arf_cmp_ui(bound,1)<0;arf_clear(bound);
  if(safe){B factor;auto quotient=polynomial.evaluate_polynomial(right)/start;acb_sqrt(factor.raw(),quotient.raw(),polynomial.bits());return root*factor;}
  auto middle=(left+right)/B(2);
  auto intermediate=continue_root(polynomial,left,middle,root,check,depth+1);
  return continue_root(polynomial,middle,right,intermediate,check,depth+1);
}
using NodeKey=std::pair<unsigned,unsigned>;
struct NodeLess {bool operator()(NodeKey a,NodeKey b)const{return a.first*b.second<b.first*a.second;}};
struct Node {B point;std::vector<B> roots,forms;};
struct Target {unsigned row;B weight;bool integer=false;slong integer_weight=0;};
struct Product {unsigned column,letter;std::vector<Target> targets;};
} // namespace spectral_detail

inline std::optional<SpectralResult> spectral_try(const Compiled& c,const Boundary& initial,
    const SpectralOptions& options,SpectralDiagnostics& diagnostics) {
  using namespace spectral_detail;
  diagnostics={};const auto started=std::chrono::steady_clock::now();auto prepared=started;bool preparation_done=false;
  const auto elapsed=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();};
  const auto finish_times=[&]{auto now=std::chrono::steady_clock::now();diagnostics.preparation_seconds=std::chrono::duration<double>((preparation_done?prepared:now)-started).count();diagnostics.numerical_seconds=preparation_done?std::chrono::duration<double>(now-prepared).count():0;};
  const auto check=[&]{if(elapsed()>options.seconds_budget)throw Rejected("spectral time budget exhausted");};
  const slong bits=B::precision(),matrix_bits=options.matrix_bits?options.matrix_bits:bits;
  try {
    if(!c.canonical)throw Rejected("spectral transport requires an epsilon-linear shared connection");
    if(c.letter_forms.size()!=c.letters.size())throw Rejected("spectral shared-letter metadata is inconsistent");
    if(!options.accuracy_goal)throw Rejected("spectral refinement requires a positive accuracy goal");
    if(matrix_bits<bits)throw Rejected("spectral matrix precision is below working precision");
    if(matrix_bits>100000)throw Rejected("spectral matrix precision budget exceeded");
    if(options.max_nodes<8 || options.max_nodes>128)throw Rejected("spectral node budget must be in 8..128");
    if(!(options.seconds_budget>0) || !std::isfinite(options.seconds_budget))throw Rejected("spectral time budget must be finite and positive");
    if((options.accuracy_goal+static_cast<double>(options.guard_digits))*3.321928094887362>bits-32)
      throw Rejected("insufficient working precision for spectral accuracy goal");
    if(initial.size()!=c.dimension)throw Rejected("spectral boundary dimension mismatch");
    for(const auto& row:initial){if(row.size()!=c.epsilon_order+1)throw Rejected("spectral boundary epsilon dimension mismatch");for(const auto& value:row)if(!value.is_finite())throw Rejected("non-finite spectral boundary");}
    B domain=B::from_strings("1/2");arb_add_error_2exp_si(acb_realref(domain.raw()),-1);
    for(const auto& point:c.singularities)
      if(arb_contains_zero(acb_imagref(point.raw())) && arb_overlaps(acb_realref(point.raw()),acb_realref(domain.raw())))
        throw Rejected("spectral interval intersects a singularity or branch point");
    if(node_forecast(c,options.accuracy_goal)>512)throw Rejected("nearby singularities require local transport");
    check();
    const unsigned d=c.dimension,k=c.epsilon_order;
    const B tolerance=B::from_strings("1e-"+std::to_string(options.accuracy_goal+options.guard_digits));
    B roundoff_floor(1);acb_mul_2exp_si(roundoff_floor.raw(),roundoff_floor.raw(),16-bits);
    const auto acceptable=[&](const B& value,const B& tail){return le(arithmetic_error(value)+B(2)*tail,tolerance*(magnitude(value)+B(1)));};
    for(const auto& row:initial)if(!acceptable(row[0],B(0)))throw Rejected("constant epsilon boundary boxes exceed spectral accuracy budget");
    std::vector<data::Expr> letters,squares;
    for(const auto& letter:c.letters){check();letters.push_back(data::Reader(letter.str()).read());}
    for(const auto& square:c.squares)squares.push_back(data::Reader(square.str()).read());
    std::vector<Product> products;std::map<std::pair<unsigned,unsigned>,unsigned> product_ids;
    for(const auto& entry:c.canonical_entries) {
      auto [it,added]=product_ids.try_emplace({entry.column,entry.letter},products.size());
      if(added)products.push_back({entry.column,entry.letter,{}});
      const auto weight=entry.coefficient.str();Target target{entry.row,B::from_strings(weight)};
      auto parsed=std::from_chars(weight.data(),weight.data()+weight.size(),target.integer_weight);
      target.integer=parsed.ec==std::errc{} && parsed.ptr==weight.data()+weight.size();
      products[it->second].targets.push_back(std::move(target));
    }
    std::vector<Jet> root_polynomials;
    std::map<NodeKey,Node,NodeLess> cache;
    {
      PrecisionScope precision(matrix_bits);Node zero;zero.point=B(0);
      for(unsigned r=0;r<squares.size();++r) {
        auto degree=c.square_polynomials[r].length();Jet x(0,degree,matrix_bits);x.set(1,B(1));
        auto polynomial=evaluate(squares[r],x,{{"x",x}});auto value=polynomial.evaluate_polynomial(B(0));B root;
        acb_sqrt(root.raw(),value.raw(),matrix_bits);
        if(!root.is_finite() || root.contains_zero())throw Rejected("spectral initial root is singular");
        root_polynomials.push_back(std::move(polynomial));zero.roots.push_back(root);
      }
      cache.emplace(NodeKey{0,1},std::move(zero));
    }
    const auto node=[&](unsigned j,unsigned n)->const Node& {
      const unsigned divisor=std::gcd(j,n);NodeKey key{j/divisor,n/divisor};
      if(auto found=cache.find(key);found!=cache.end())return found->second;
      check();if((cache.size()+1)*(letters.size()+squares.size()+1)>options.cell_budget)throw Rejected("spectral node-sample cell budget exhausted");
      PrecisionScope precision(matrix_bits);Node value;B angle=B::from_strings(std::to_string(key.first)+"/"+std::to_string(key.second)),cosine;
      acb_cos_pi(cosine.raw(),angle.raw(),matrix_bits);value.point=(B(1)-cosine)/B(2);
      auto left=std::prev(cache.lower_bound(key));
      Jet x(0,2,matrix_bits);x.set(0,value.point);x.set(1,B(1));std::map<std::string,Jet> vars{{"x",x}};
      for(unsigned r=0;r<squares.size();++r) {
        auto root=continue_root(root_polynomials[r],left->second.point,value.point,left->second.roots[r],check);
        if(!root.is_finite() || root.contains_zero())throw Rejected("spectral sample root is singular");
        value.roots.push_back(root);auto square=evaluate(squares[r],x,vars);auto jet=x.constant(0);
        jet.set(0,root);jet.set(1,square.at(1)/(B(2)*root));vars.emplace("r"+std::to_string(r),std::move(jet));
      }
      for(unsigned l=0;l<letters.size();++l) {
        if(l%32==0)check();auto jet=evaluate(letters[l],x,vars);
        auto scalar=c.letter_forms[l]?jet.at(0):jet.at(1)/jet.at(0);
        if(!scalar.is_finite())throw Rejected("spectral matrix has a non-finite sample");value.forms.push_back(std::move(scalar));
      }
      auto inserted=cache.emplace(key,std::move(value));diagnostics.sampled_nodes=cache.size()-1;return inserted.first->second;
    };
    prepared=std::chrono::steady_clock::now();preparation_done=true;
    std::uint64_t work=0;Boundary previous,previous_difference;unsigned previous_nodes=0,previous_increment=0;
    std::string last_reason="three refinement levels are required";
    for(unsigned n:std::array<unsigned,9>{8,12,16,24,32,48,64,96,128}) {
      if(n>options.max_nodes)break;check();const unsigned m=n+1;
      const std::uint64_t candidate=static_cast<std::uint64_t>(m)*m*m+static_cast<std::uint64_t>(k)*n*(c.canonical_entries.size()+products.size()+static_cast<std::uint64_t>(m)*d);
      if(candidate>options.work_budget || work>options.work_budget-candidate)throw Rejected("spectral arithmetic work budget exhausted");work+=candidate;
      if(static_cast<std::uint64_t>(m)*(3ull*d+2ull*m)>options.cell_budget)throw Rejected("spectral collocation cell budget exhausted");
      diagnostics.nodes_tried.push_back(n);
      std::vector<const Node*> nodes(m);nodes[0]=&cache.begin()->second;for(unsigned j=1;j<m;++j)nodes[j]=&node(j,n);
      Matrix derivative(m,m),integral(m,m),state(m,d),rhs(m,d),next(m,d);
      std::vector<B> weights(m);for(unsigned j=0;j<m;++j){weights[j]=B(j%2?-1:1);if(j==0||j==n)weights[j]=weights[j]/B(2);}
      acb_one(acb_mat_entry(derivative.value,0,0));
      for(unsigned row=1;row<m;++row){check();B diagonal(0);for(unsigned column=0;column<m;++column)if(column!=row){auto entry=weights[column]/weights[row]/(nodes[row]->point-nodes[column]->point);acb_set(acb_mat_entry(derivative.value,row,column),entry.raw());diagonal=diagonal-entry;}acb_set(acb_mat_entry(derivative.value,row,row),diagonal.raw());}
      if(!acb_mat_inv(integral.value,derivative.value,bits))throw Rejected("spectral scalar collocation inverse is unresolved");check();
      Boundary values(d,std::vector<B>(k+1,B(0)));
      for(unsigned row=0;row<d;++row){values[row][0]=initial[row][0];for(unsigned j=0;j<m;++j)acb_set(acb_mat_entry(state.value,j,row),initial[row][0].raw());}
      for(unsigned epsilon=1;epsilon<=k;++epsilon) {
        acb_mat_zero(rhs.value);B product;
        for(unsigned j=1;j<m;++j){check();unsigned group=0;for(const auto& p:products){if((++group%256)==0)check();acb_mul(product.raw(),nodes[j]->forms[p.letter].raw(),acb_mat_entry(state.value,j,p.column),bits);for(const auto& target:p.targets){auto destination=acb_mat_entry(rhs.value,j,target.row);if(target.integer)acb_addmul_si(destination,product.raw(),target.integer_weight,bits);else acb_addmul(destination,target.weight.raw(),product.raw(),bits);}}}
        acb_mat_mul(next.value,integral.value,rhs.value,bits);check();
        // The inverse's boundary column is exactly the constant function one;
        // add it directly, avoiding artificial loss of boundary-box correlation.
        for(unsigned j=0;j<m;++j)for(unsigned row=0;row<d;++row)acb_add(acb_mat_entry(next.value,j,row),acb_mat_entry(next.value,j,row),initial[row][epsilon].raw(),bits);
        for(unsigned row=0;row<d;++row){acb_set(values[row][epsilon].raw(),acb_mat_entry(next.value,n,row));if(!values[row][epsilon].is_finite())throw Rejected("non-finite spectral endpoint");}
        acb_mat_swap(state.value,next.value);
      }
      if(previous.empty()){previous=std::move(values);previous_nodes=n;continue;}
      Boundary current_difference(d,std::vector<B>(k+1,B(0)));
      for(unsigned row=0;row<d;++row)for(unsigned epsilon=1;epsilon<=k;++epsilon)current_difference[row][epsilon]=difference(values[row][epsilon],previous[row][epsilon]);
      const unsigned increment=n-previous_nodes;
      if(!previous_difference.empty()) {
        Boundary tails(d,std::vector<B>(k+1,B(0)));bool converged=true,boxes_ok=true;unsigned absolute_stability=0;
        for(unsigned row=0;row<d;++row)for(unsigned epsilon=0;epsilon<=k;++epsilon) {
          const auto& value=values[row][epsilon];
          if(!acceptable(value,B(0)))boxes_ok=false;
          if(epsilon==0)continue;
          const auto floor=roundoff_floor*(midpoint_scale(value)+midpoint_scale(previous[row][epsilon]));
          const auto& before=previous_difference[row][epsilon];const auto& after=current_difference[row][epsilon];
          std::optional<double> factor;
          if(le(before+after,floor))factor=0;
          else if(!before.is_zero())factor=tail_factor(ratio_upper(after,before),previous_increment,increment);
          if(!factor) {
            ++absolute_stability;
            // Relative decay is meaningless for tiny cancellation/noise terms.
            // Demand two small absolute changes with a safety factor instead.
            tails[row][epsilon]=upper(B(16)*(before+after)+floor);
          } else {
            B multiplier;acb_set_d(multiplier.raw(),*factor);tails[row][epsilon]=upper(after*multiplier+floor);
          }
          if(!acceptable(value,tails[row][epsilon])){converged=false;last_reason="estimated spectral tail exceeds accuracy budget";}
        }
        if(!boxes_ok)throw Rejected("spectral arithmetic or input boxes exceed accuracy budget");
        if(converged){diagnostics.absolute_stability_components=absolute_stability;diagnostics.selected_nodes=n;finish_times();return SpectralResult{std::move(values),std::move(tails),diagnostics,nodes.back()->roots};}
      }
      previous=std::move(values);previous_difference=std::move(current_difference);previous_nodes=n;previous_increment=increment;
    }
    throw Rejected("spectral node budget exhausted: "+last_reason);
  } catch(const Rejected& error){diagnostics.rejection_reason=error.what();}
    catch(const std::exception& error){diagnostics.rejection_reason=std::string("spectral evaluation failed: ")+error.what();}
  finish_times();return std::nullopt;
}
} // namespace diffexp::transport

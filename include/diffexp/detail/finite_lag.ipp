#pragma once
// Included after transport's Compiled and Chart types. All plans are derived
// from the supplied connection; there is no family-name dispatch.
namespace diffexp::transport {
struct FiniteLagProduct {
  unsigned column,epsilon;
  unsigned coefficient;
  std::vector<std::pair<unsigned,B>> targets;
};
struct FiniteLagPlan {
  std::vector<unsigned> gauge;
  std::vector<FiniteLagProduct> products;
  std::vector<std::pair<std::vector<B>,std::vector<B>>> coefficients;
  unsigned max_degree=0;
};
namespace finite_lag_detail {
inline std::map<unsigned,Exact> decompose(const transport::Compiled& c,Exact q){
 auto num=q.numerator(),den=q.denominator();
 for(unsigned r=c.squares.size();r-->0;){
  bool present=false;for(const auto& t:den.numerator_terms())if(t.powers[2+r])present=true;
  if(!present)continue;
  std::vector<Exact> subs;for(const auto& name:q.variables())subs.emplace_back(c.field,name);
  subs[2+r]=-subs[2+r];auto conjugate=den.substitute(subs);
  num=c.reduced(num*conjugate);den=c.reduced(den*conjugate);
  auto reduced=c.reduced(num/den);num=reduced.numerator();den=reduced.denominator();
 }
 for(const auto& t:den.numerator_terms())for(unsigned r=0;r<c.squares.size();++r)if(t.powers[2+r])throw std::runtime_error("unrationalized denominator");
 std::map<unsigned,Exact> out;
 for(const auto& t:num.numerator_terms()){
  unsigned mask=0;Exact term(c.field,t.coefficient.str());
  for(unsigned v=0;v<t.powers.size();++v){if(v<2)term=term*term.variable(v).pow(t.powers[v]);else{if(t.powers[v]>1)throw std::runtime_error("unreduced root power");if(t.powers[v])mask|=1u<<(v-2);}}
  auto [it,added]=out.try_emplace(mask,Exact(c.field,0));it->second=it->second+term/den;
 }
 Exact reconstructed(c.field,0);for(const auto& [m,a]:out){auto term=a;for(unsigned r=0;r<c.squares.size();++r)if(m&(1u<<r))term=term*term.variable(r+2);reconstructed=reconstructed+term;}
 if(!c.reduced(reconstructed-q).is_zero())throw std::runtime_error("exact algebraic reconstruction failed");
 return out;
}
inline std::vector<B> polynomial(const Exact& a) {
  if(!a.denominator().is_rational())throw std::runtime_error("nonpolynomial finite-lag coefficient");
  std::vector<B> out(1,B(0));
  for(const auto& t:a.numerator_terms()) {
    for(unsigned v=2;v<t.powers.size();++v)if(t.powers[v])throw std::runtime_error("remaining algebraic coefficient");
    if(t.powers[0]>256)throw std::runtime_error("finite-lag polynomial degree budget");
    const unsigned n=t.powers[0],ip=t.powers[1]%4;
    out.resize(std::max(out.size(),std::size_t(n+1)),B(0));
    auto scalar=B::from_strings((t.coefficient/a.denominator().rational()).str());
    if(ip>=2)scalar=-scalar;
    if(ip%2){B z;arb_set(acb_imagref(z.raw()),acb_realref(scalar.raw()));scalar=std::move(z);}
    out[n]+=scalar;
  }
  while(out.size()>1&&out.back().is_zero())out.pop_back();
  return out;
}
// Polynomial translation costs depend on degree, not requested Taylor order.
inline std::vector<B> shifted(std::vector<B> p,const B& center) {
  for(std::size_t i=p.size()-1;i-->0;)for(std::size_t j=i;j+1<p.size();++j)
    acb_addmul(p[j].raw(),center.raw(),p[j+1].raw(),B::precision());
  return p;
}
inline std::shared_ptr<const FiniteLagPlan> prepare(const Compiled& c) {
  if(!c.canonical)throw std::runtime_error("noncanonical connection");
  if(c.squares.size()>8 || c.canonical_entries.size()>100000)
    throw std::runtime_error("finite-lag preparation budget");
  const auto started=std::chrono::steady_clock::now();
  const auto check_time=[&]{if(std::chrono::steady_clock::now()-started>std::chrono::seconds(10))
    throw std::runtime_error("finite-lag preparation time budget");};
  std::vector<std::map<unsigned,Exact>> parts;
  for(unsigned letter=0;letter<c.letters.size();++letter) {
    check_time();auto p=decompose(c,c.letter_derivative(letter));
    for(auto it=p.begin();it!=p.end();)if(it->second.is_zero())it=p.erase(it);else ++it;
    if(p.size()>1)throw std::runtime_error("letter has multiple root parities");
    parts.push_back(std::move(p));
  }
  auto plan=std::make_shared<FiniteLagPlan>();
  std::vector<std::vector<std::pair<unsigned,unsigned>>> graph(c.dimension);
  for(const auto& e:c.canonical_entries)if(!e.coefficient.is_zero())for(const auto& [mask,a]:parts[e.letter]) {
    graph[e.row].push_back({e.column,mask});graph[e.column].push_back({e.row,mask});
  }
  const unsigned unset=~0u;plan->gauge.assign(c.dimension,unset);
  for(unsigned i=0;i<c.dimension;++i)if(plan->gauge[i]==unset) {
    plan->gauge[i]=0;std::vector<unsigned> queue{i};
    for(std::size_t at=0;at<queue.size();++at)for(auto [col,mask]:graph[queue[at]]) {
      unsigned g=plan->gauge[queue[at]]^mask;
      if(plan->gauge[col]==unset){plan->gauge[col]=g;queue.push_back(col);}
      else if(plan->gauge[col]!=g)throw std::runtime_error("inconsistent diagonal root gauge");
    }
  }
  std::vector<Exact> squares;
  for(unsigned mask=0;mask<(1u<<c.squares.size());++mask) {
    Exact p(c.field,1);for(unsigned r=0;r<c.squares.size();++r)if(mask&(1u<<r))p=p*c.squares[r];
    squares.push_back(std::move(p));
  }
  std::map<std::tuple<unsigned,unsigned,std::string>,unsigned> ids;
  std::map<std::string,unsigned> coefficient_ids;
  const auto add=[&](unsigned row,unsigned col,unsigned eps,const Exact& coefficient,const Rational& weight) {
    if(coefficient.is_zero() || weight.is_zero())return;
    check_time();auto key=std::make_tuple(col,eps,coefficient.str());
    auto [it,inserted]=ids.try_emplace(key,plan->products.size());
    if(inserted) {
      auto [ci,new_coefficient]=coefficient_ids.try_emplace(std::get<2>(key),plan->coefficients.size());
      if(new_coefficient) {
        plan->coefficients.emplace_back(polynomial(coefficient.numerator()),polynomial(coefficient.denominator()));
        const auto& [p,q]=plan->coefficients.back();
        plan->max_degree=std::max(plan->max_degree,unsigned(std::max(p.size(),q.size())-1));
      }
      plan->products.push_back({col,eps,ci->second,{}});
    }
    plan->products[it->second].targets.emplace_back(row,B::from_strings(weight.str()));
  };
  for(const auto& e:c.canonical_entries)for(const auto& [mask,a]:parts[e.letter])
    add(e.row,e.column,1,c.reduced(a*squares[plan->gauge[e.row]&mask]),e.coefficient);
  for(unsigned i=0;i<c.dimension;++i)if(plan->gauge[i]) {
    const auto& p=squares[plan->gauge[i]];
    add(i,i,0,c.reduced(p.derivative(0)/(p*Exact(c.field,2))),Rational(1));
  }
  return plan;
}
} // namespace finite_lag_detail
inline const FiniteLagPlan* finite_lag_plan(const Compiled& c) {
  if(!c.finite_lag_enabled)return nullptr;
  if(!c.finite_lag_attempted || c.finite_lag_bits!=B::precision()) {
    c.finite_lag_attempted=true;c.finite_lag_bits=B::precision();c.finite_lag.reset();
    try {c.finite_lag=finite_lag_detail::prepare(c);c.finite_lag_reason.clear();}
    catch(const std::runtime_error& e){c.finite_lag_reason=e.what();}
  }
  return c.finite_lag.get();
}
inline Chart finite_lag_chart(const Compiled& c,const FiniteLagPlan& plan,const Boundary& initial,
    const std::vector<B>& roots,const B& center,const B& step,unsigned order) {
  const unsigned d=c.dimension,w=c.epsilon_order+1;const auto bits=B::precision();
  if(roots.size()!=c.squares.size())throw std::runtime_error("finite-lag root dimensions");
  std::vector<B> factors(d,B(1)),endpoint_factors(d,B(1)),end_roots;
  for(unsigned r=0;r<roots.size();++r)
    end_roots.push_back(continue_polynomial_sqrt(c.square_polynomials[r],center,center+step,roots[r]));
  for(unsigned i=0;i<d;++i)for(unsigned r=0;r<roots.size();++r)if(plan.gauge[i]&(1u<<r)) {
    factors[i]=factors[i]*roots[r];endpoint_factors[i]=endpoint_factors[i]*end_roots[r];
  }
  for(unsigned i=0;i<d;++i)if(factors[i].contains_zero() || endpoint_factors[i].contains_zero())
    throw std::runtime_error("noninvertible finite-lag gauge");
  const auto& groups=plan.products;
  auto coefficients=plan.coefficients;
  // Translate each distinct rational coefficient once, sharing it over sources.
  for(auto& [p,q]:coefficients) {
    p=finite_lag_detail::shifted(std::move(p),center);
    q=finite_lag_detail::shifted(std::move(q),center);
    if(q[0].contains_zero())throw std::runtime_error("singular finite-lag denominator");
    // A finite-lag inverse can amplify independent interval errors even when
    // the physical connection is regular. Bound its scaled feedback norm
    // before allocating/solving: sum |q_j h^j| / |q_0| <= 3/4.
    B feedback,power(1);
    for(unsigned j=1;j<q.size();++j){power=power*step;feedback+=magnitude(q[j]*power);}
    if(!arb_le(acb_realref(feedback.raw()),acb_realref((magnitude(q[0])*B::from_strings("3/4")).raw())))
      throw std::runtime_error("finite-lag denominator feedback exceeds stability budget");
  }
  std::uint64_t cells=static_cast<std::uint64_t>(order+1)*d*w;
  for(const auto& g:groups)cells+=coefficients[g.coefficient].second.size()*w;
  if(cells>10000000)throw std::runtime_error("finite-lag workspace budget");
  AcbArray y(static_cast<slong>(order+1)*d*w);
  std::vector<std::unique_ptr<AcbArray>> buffers;
  for(const auto& g:groups)buffers.push_back(std::make_unique<AcbArray>(coefficients[g.coefficient].second.size()*w));
  const auto at=[&](unsigned n,unsigned i,unsigned k){return y.p+(static_cast<slong>(n)*d+i)*w+k;};
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<w;++k)acb_mul(at(0,i,k),initial[i][k].raw(),factors[i].raw(),bits);
  B value;
  for(unsigned n=0;n<order;++n) {
    for(unsigned gi=0;gi<groups.size();++gi) {
      const auto& g=groups[gi];const auto& [p,q]=coefficients[g.coefficient];
      auto f=buffers[gi]->p;const unsigned depth=q.size();
      for(unsigned k=0;k+g.epsilon<w;++k) {
        acb_zero(value.raw());
        for(unsigned m=0;m<p.size()&&m<=n;++m)if(!p[m].is_zero())
          acb_addmul(value.raw(),p[m].raw(),at(n-m,g.column,k),bits);
        for(unsigned m=1;m<depth&&m<=n;++m)if(!q[m].is_zero())
          acb_submul(value.raw(),q[m].raw(),f+((n-m)%depth)*w+k,bits);
        acb_div(value.raw(),value.raw(),q[0].raw(),bits);
        acb_set(f+(n%depth)*w+k,value.raw());
        for(const auto& [row,weight]:g.targets)acb_addmul(at(n+1,row,k+g.epsilon),value.raw(),weight.raw(),bits);
      }
    }
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<w;++k)acb_div_ui(at(n+1,i,k),at(n+1,i,k),n+1,bits);
  }
  Chart result{Boundary(d,std::vector<B>(w)),Boundary(d,std::vector<B>(w)),Boundary(d,std::vector<B>(w)),{},false,true};
  std::vector<B> powers(order+1,B(1));for(unsigned n=1;n<=order;++n)powers[n]=powers[n-1]*step;
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<w;++k) {
    auto& v=result.values[i][k];
    for(unsigned n=order+1;n-->0;){acb_mul(v.raw(),v.raw(),step.raw(),bits);acb_add(v.raw(),v.raw(),at(n,i,k),bits);}
    v=v/endpoint_factors[i];
    if(!v.is_finite())throw std::runtime_error("nonfinite finite-lag output");
    // Estimate the omitted transformed series in the physical endpoint basis.
    // Input radii remain present in every term; this is not a tail certificate.
    for(unsigned n=order-3;n<=order;++n) {
      B term;acb_mul(term.raw(),at(n,i,k),powers[n].raw(),bits);
      result.truncation_errors[i][k]+=magnitude(term/endpoint_factors[i]);
    }
    result.errors[i][k]=result.truncation_errors[i][k]+arithmetic_error(v);
    if(!result.errors[i][k].is_finite())throw std::runtime_error("nonfinite finite-lag error estimate");
  }
  return result;
}
} // namespace diffexp::transport

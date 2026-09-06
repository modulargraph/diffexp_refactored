#pragma once
struct WeightedTerm {unsigned row,column,epsilon;Exact coefficient;Rational weight;};
struct ProductGroup{unsigned column,epsilon;std::vector<B> p,q;std::vector<std::pair<unsigned,B>> targets;};
std::vector<B> complex_polynomial(const Exact&a){
 if(!a.denominator().is_rational())throw std::runtime_error("polynomial required");
 std::vector<B> out(1,B(0));for(auto&t:a.numerator_terms()){
  for(unsigned v=2;v<t.powers.size();++v)if(t.powers[v])throw std::runtime_error("remaining algebraic symbol");
  if(t.powers[0]>1024)throw std::runtime_error("degree budget");
  unsigned n=t.powers[0],ip=t.powers[1]%4;out.resize(std::max(out.size(),std::size_t(n+1)),B(0));
  auto scalar=B::from_strings((t.coefficient/a.denominator().rational()).str());if(ip>=2)scalar=-scalar;if(ip%2){B z;arb_set(acb_imagref(z.raw()),acb_realref(scalar.raw()));scalar=std::move(z);}out[n]=out[n]+scalar;
 }while(out.size()>1&&out.back().is_zero())out.pop_back();return out;
}
std::vector<ProductGroup> prepare_products(const std::vector<WeightedTerm>& terms){
 std::map<std::tuple<unsigned,unsigned,std::string>,unsigned> ids;std::vector<ProductGroup> out;
 for(auto&e:terms){auto key=std::make_tuple(e.column,e.epsilon,e.coefficient.str());auto[it,added]=ids.try_emplace(key,out.size());if(added){auto p=complex_polynomial(e.coefficient.numerator()),q=complex_polynomial(e.coefficient.denominator());if(q[0].contains_zero())throw std::runtime_error("singular product denominator");out.push_back({e.column,e.epsilon,std::move(p),std::move(q),{}});}out[it->second].targets.push_back({e.row,B::from_strings(e.weight.str())});}return out;
}
Boundary product_chart(const std::vector<ProductGroup>& groups,const Boundary& initial,const B& step,unsigned N,std::vector<Boundary>* retained=nullptr){
 unsigned d=initial.size(),w=initial[0].size();slong bits=B::precision();
 if(std::uint64_t(N+1)*(d+groups.size())*w>10000000)throw std::runtime_error("product workspace budget");
 transport::AcbArray y((N+1)*d*w);std::vector<std::unique_ptr<transport::AcbArray>> buffers;
 for(auto&g:groups)buffers.push_back(std::make_unique<transport::AcbArray>(g.q.size()*w));
 auto at=[&](unsigned n,unsigned i,unsigned k){return y.p+(n*d+i)*w+k;};
 for(unsigned i=0;i<d;++i)for(unsigned k=0;k<w;++k)acb_set(at(0,i,k),initial[i][k].raw());
 B val;
 for(unsigned n=0;n<N;++n){
  for(unsigned gi=0;gi<groups.size();++gi){const auto&g=groups[gi];auto f=buffers[gi]->p;unsigned depth=g.q.size();
   for(unsigned k=0;k+g.epsilon<w;++k){acb_zero(val.raw());
    for(unsigned m=0;m<g.p.size()&&m<=n;++m)if(!g.p[m].is_zero())acb_addmul(val.raw(),g.p[m].raw(),at(n-m,g.column,k),bits);
    for(unsigned m=1;m<depth&&m<=n;++m)if(!g.q[m].is_zero())acb_submul(val.raw(),g.q[m].raw(),f+((n-m)%depth)*w+k,bits);
    acb_div(val.raw(),val.raw(),g.q[0].raw(),bits);acb_set(f+(n%depth)*w+k,val.raw());
    for(auto&[row,weight]:g.targets)acb_addmul(at(n+1,row,k+g.epsilon),val.raw(),weight.raw(),bits);
   }
  }
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<w;++k)acb_div_ui(at(n+1,i,k),at(n+1,i,k),n+1,bits);
 }
 if(retained){retained->assign(N+1,Boundary(d,std::vector<B>(w)));for(unsigned n=0;n<=N;++n)for(unsigned i=0;i<d;++i)for(unsigned k=0;k<w;++k)acb_set((*retained)[n][i][k].raw(),at(n,i,k));}
 Boundary out(d,std::vector<B>(w));for(unsigned i=0;i<d;++i)for(unsigned k=0;k<w;++k)for(unsigned n=N+1;n-->0;){acb_mul(out[i][k].raw(),out[i][k].raw(),step.raw(),bits);acb_add(out[i][k].raw(),out[i][k].raw(),at(n,i,k),bits);}return out;
}

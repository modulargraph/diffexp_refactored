#pragma once
#include "diffexp/jet.hpp"
#include "diffexp/exact.hpp"
#include "diffexp/reference_comparison.hpp"
#include <chrono>
#include <iostream>
#include <tuple>

namespace diffexp {
struct CanonicalEntry {unsigned row,column,letter;Rational coefficient;};
using Boundary=std::vector<std::vector<kernel::ComplexBall>>;

inline std::map<int,Rational> linear_logs(const data::Expr& e) {
  using Form=std::map<int,Rational>;
  if(data::number(e))return {{-1,Rational(e.head)}};
  if(e.head=="Log" && e.args.size()==1 && e.args[0].head=="W" && e.args[0].args.size()==1)
    return {{static_cast<int>(data::integer(e.args[0].args[0]))-1,Rational(1)}};
  if(e.head=="neg") {auto a=linear_logs(e.args.at(0));for(auto& [i,c]:a)c=-c;return a;}
  if(e.args.size()!=2)throw std::invalid_argument("ancillary matrix is not linear in dlogs");
  auto a=linear_logs(e.args[0]),b=linear_logs(e.args[1]);
  if(e.head=="+" || e.head=="-") {
    for(const auto& [i,c]:b) {auto p=a.try_emplace(i,Rational(0)).first;p->second+=e.head=="+"?c:-c;}
    return a;
  }
  const auto scalar=[](const Form& f){return f.size()==1 && f.begin()->first==-1;};
  if(e.head=="*" && (scalar(a)||scalar(b))) {
    if(!scalar(b))std::swap(a,b);
    for(auto& [i,c]:a)c*=b.begin()->second;return a;
  }
  if(e.head=="/" && scalar(b)) {for(auto& [i,c]:a)c=c/b.begin()->second;return a;}
  throw std::invalid_argument("nonlinear log expression in canonical matrix");
}

inline std::vector<CanonicalEntry> read_dlog_matrix(const data::Expr& matrix,unsigned dimension,unsigned letters) {
  if(matrix.head!="List" || matrix.args.size()!=dimension)throw std::invalid_argument("canonical matrix shape");
  std::vector<CanonicalEntry> out;
  for(unsigned i=0;i<dimension;++i) {
    const auto& row=matrix.args[i];
    if(row.head!="List" || row.args.size()!=dimension)throw std::invalid_argument("canonical matrix row shape");
    for(unsigned j=0;j<dimension;++j)for(const auto& [letter,c]:linear_logs(row.args[j])) {
      if(c.is_zero())continue;
      if(letter<0 || static_cast<unsigned>(letter)>=letters)throw std::invalid_argument("invalid dlog letter index");
      out.push_back({i,j,static_cast<unsigned>(letter),c});
    }
  }
  return out;
}

inline Boundary read_boundary(const data::Expr& input,unsigned dimension,unsigned epsilon_order,slong bits) {
  Jet context(0,1,bits);Boundary out;
  if(input.head!="List" || input.args.size()!=dimension)throw std::invalid_argument("boundary shape");
  for(const auto& row:input.args) {
    if(row.head!="List" || row.args.size()!=epsilon_order+1)throw std::invalid_argument("boundary epsilon shape");
    std::vector<Jet::Ball> values;
    for(const auto& e:row.args)values.push_back(evaluate(e,context,{}).at(0));
    out.push_back(std::move(values));
  }
  return out;
}

// Truncated epsilon-triangular recurrence. A finite regular boundary at a
// singular center is supported; general logarithmic starts use the Frobenius
// kernel. This routine reports numerical values, never a truncation certificate.
inline Boundary canonical_chart(const std::vector<CanonicalEntry>& entries,
    const std::vector<Jet>& letters,const Boundary& incoming,
    const Jet::Ball& input_offset,const Jet::Ball& output_offset,unsigned order) {
  using B=Jet::Ball;
  if(letters.empty() || incoming.empty() || !order)throw std::invalid_argument("empty canonical chart");
  if(incoming[0].empty())throw std::invalid_argument("empty canonical epsilon window");
  for(const auto& row:incoming)if(row.size()!=incoming[0].size())throw std::invalid_argument("canonical boundary row shape");
  unsigned d=incoming.size(),eps=incoming[0].size()-1;
  const auto bits=letters[0].bits();
  if(letters[0].length()<order+3)throw std::invalid_argument("insufficient letter jet order");
  struct Entry {unsigned row,col;std::vector<B> regular;B residue;};
  std::map<std::pair<unsigned,unsigned>,Entry> combined;
  std::vector<Jet> derivatives;std::vector<unsigned> valuations;
  for(const auto& w:letters) {
    unsigned v=0;while(v<w.length() && w.at(v).is_zero())++v;
    if(v+order+2>=w.length())throw std::domain_error("insufficient nonzero letter prefix");
    if(v && !input_offset.is_zero())throw std::domain_error("singular chart needs a center boundary");
    auto u=w.shifted_down(v);derivatives.push_back(u.derivative()/u);valuations.push_back(v);
  }
  for(const auto& e:entries) {
    if(e.row>=d || e.column>=d || e.letter>=letters.size())throw std::invalid_argument("canonical entry shape");
    auto [p,inserted]=combined.try_emplace({e.row,e.column},Entry{e.row,e.column,std::vector<B>(order,B(0)),B(0)});
    const B c=B::from_strings(e.coefficient.str());
    for(unsigned n=0;n<order;++n) {
      auto term=derivatives[e.letter].at(n);
      if(!term.is_zero())acb_addmul(p->second.regular[n].raw(),c.raw(),term.raw(),bits);
    }
    if(valuations[e.letter])p->second.residue+=c*B(valuations[e.letter]);
  }
  std::vector<Entry> a;for(auto& [key,value]:combined)a.push_back(std::move(value));
  std::vector<B> coefficients(static_cast<std::size_t>(eps+1)*(order+1)*d,B(0));
  const auto at=[&](unsigned k,unsigned n,unsigned i)->B& {return coefficients[(static_cast<std::size_t>(k)*(order+1)+n)*d+i];};
  Boundary out(d,std::vector<B>(eps+1,B(0)));
  for(unsigned k=0;k<=eps;++k) {
    if(k) {
      for(unsigned n=1;n<=order;++n) {
        for(const auto& e:a) {
          auto& target=at(k,n,e.row);
          if(!e.residue.is_zero() && !at(k-1,n,e.col).is_zero())
            acb_addmul(target.raw(),e.residue.raw(),at(k-1,n,e.col).raw(),bits);
          for(unsigned m=0;m<n;++m)
            if(!e.regular[m].is_zero() && !at(k-1,n-1-m,e.col).is_zero())
              acb_addmul(target.raw(),e.regular[m].raw(),at(k-1,n-1-m,e.col).raw(),bits);
        }
        for(unsigned i=0;i<d;++i)acb_div_ui(at(k,n,i).raw(),at(k,n,i).raw(),n,bits);
      }
    }
    for(unsigned i=0;i<d;++i) {
      B initial(0),value(0);
      for(unsigned n=order;n>0;--n)initial=(initial+at(k,n,i))*input_offset;
      at(k,0,i)=incoming.at(i).at(k)-initial;
      for(unsigned n=order+1;n-->0;)value=value*output_offset+at(k,n,i);
      if(!value.is_finite())throw std::runtime_error("non-finite canonical coefficient");
      out[i][k]=std::move(value);
    }
    if(k<eps) {
      std::vector<B> residue_action(d,B(0));
      for(const auto& e:a)if(!e.residue.is_zero())
        acb_addmul(residue_action[e.row].raw(),e.residue.raw(),at(k,0,e.col).raw(),bits);
      for(const auto& r:residue_action)
        if(!r.contains_zero())throw std::runtime_error("boundary incompatible with a finite regular singular start");
    }
  }
  return out;
}

inline std::vector<Jet> henn_letters(const Jet& t) {
  std::vector<Jet> v;
  const char* start[]={"3","-1","1","1","-1"};
  const char* finish[]={"4","-113/47","281/149","349/257","-863/541"};
  for(unsigned i=0;i<5;++i) {
    auto s=t.constant(0);s.set(0,Jet::Ball::from_strings(start[i]));
    auto f=t.constant(0);f.set(0,Jet::Ball::from_strings(finish[i]));
    v.push_back(s+(f-s)*t);
  }
  auto [v1,v2,v3,v4,v5]=std::tuple{v[0],v[1],v[2],v[3],v[4]};
  auto delta=v1*v1*(v2-v5).pow(2)+(v2*v3+v4*(-v3+v5)).pow(2)+t.constant(2)*v1*(-v2*v2*v3+v4*(v3-v5)*v5+v2*(v3*v4+(v3+v4)*v5));
  auto root=delta.sqrt();
  std::vector<Jet> letters{
    v1,v2,v3,v4,v5,v3+v4,v4+v5,v1+v5,v1+v2,v2+v3,
    v1-v4,v2-v5,-v1+v3,-v2+v4,-v3+v5,
    v1+v2-v4,v2+v3-v5,-v1+v3+v4,-v2+v4+v5,v1-v3+v5,
    -v1-v2+v3+v4,-v2-v3+v4+v5,v1-v3-v4+v5,v1+v2-v4-v5,-v1+v2+v3-v5};
  for(const auto& base:std::vector<Jet>{
      v1*v2-v2*v3+v3*v4-v1*v5-v4*v5,-v1*v2+v2*v3-v3*v4-v1*v5+v4*v5,
      -v1*v2-v2*v3+v3*v4+v1*v5-v4*v5,v1*v2-v2*v3-v3*v4-v1*v5+v4*v5,
      -v1*v2+v2*v3-v3*v4+v1*v5-v4*v5})letters.push_back((base-root)/(base+root));
  letters.push_back(root);return letters;
}

inline int run_henn_nonplanar(const std::string& data_dir,unsigned order=140) {
  using B=Jet::Ball;
  if(order<10 || order>500)throw std::invalid_argument("Henn Taylor order must be between 10 and 500");
  B::set_precision(256);
  auto start=std::chrono::steady_clock::now();
  auto matrix=read_dlog_matrix(data::read_file(data_dir+"/XB_Atilde.txt"),108,31);
  auto value=read_boundary(data::read_file(data_dir+"/XB_Boundary_values_X0.txt"),108,4,256);
  auto reference=read_boundary(data::read_file(data_dir+"/XB_Boundary_values_X1.txt"),108,4,256);
  const char* centers[]={"0","1/3","7/10","19/20"};
  const char* bounds[]={"0","0.19698903994874387","0.5085056419319689","0.8409884239904127","1"};
  for(unsigned chart=0;chart<4;++chart) {
    Jet t(0,order+4,256);t.set(0,B::from_strings(centers[chart]));t.set(1,B(1));
    std::cout<<"Henn chart "<<chart+1<<"/4, Taylor order "<<order<<std::endl;
    value=canonical_chart(matrix,henn_letters(t),value,
      B::from_strings(bounds[chart])-B::from_strings(centers[chart]),
      B::from_strings(bounds[chart+1])-B::from_strings(centers[chart]),order);
  }
  double maximum=0;
  for(unsigned i=0;i<108;++i)for(unsigned k=0;k<=4;++k) {
    maximum=std::max(maximum,finite_reference_error(value[i][k],reference[i][k],256));
  }
  std::cout<<"108 masters, epsilon 0..4; maximum reference discrepancy = "<<maximum
           <<"; seconds = "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
  std::cout<<"Reference comparison only; finite Taylor truncation is not certified by this driver.\n";
  return std::isfinite(maximum) && maximum<1e-8 ? 0:1;
}
}  // namespace diffexp

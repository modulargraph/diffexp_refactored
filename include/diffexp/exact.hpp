#pragma once

#include "diffexp/kernel/scalar.hpp"
#include "diffexp/data.hpp"
#include <flint/fmpz_mpoly_q.h>
#include <flint/fmpz_poly.h>
#include <algorithm>
#include <cctype>
#include <memory>
#include <set>
#include <span>
#include <utility>

namespace diffexp {
using Rational = kernel::Rational;

// The field is owned by its values. Independent families can coexist in one
// process: neither a global variable list nor a kernel restart is necessary.
class ExactField {
 public:
  struct Context {
    std::vector<std::string> names;
    std::vector<const char*> symbols;
    fmpz_mpoly_ctx_t raw;
    explicit Context(std::vector<std::string> variables) : names(std::move(variables)) {
      if (names.empty()) throw std::invalid_argument("an exact field needs variables");
      std::set<std::string> seen;
      for (const auto& name : names) {
        if (name.empty() || !std::isalpha(static_cast<unsigned char>(name[0])) ||
            !std::all_of(name.begin(), name.end(), [](unsigned char c) {
              return std::isalnum(c) || c == '_';
            }) || !seen.insert(name).second)
          throw std::invalid_argument("invalid or duplicate field variable: " + name);
        symbols.push_back(name.c_str());
      }
      fmpz_mpoly_ctx_init(raw, static_cast<slong>(names.size()), ORD_DEGLEX);
    }
    Context(const Context&) = delete;
    ~Context() { fmpz_mpoly_ctx_clear(raw); }
  };
  explicit ExactField(std::vector<std::string> variables)
      : context_(std::make_shared<Context>(std::move(variables))) {}
  const std::vector<std::string>& variables() const { return context_->names; }
  std::size_t index(const std::string& name) const {
    auto i = std::find(variables().begin(), variables().end(), name);
    if (i == variables().end()) throw std::invalid_argument("undeclared variable: " + name);
    return static_cast<std::size_t>(i - variables().begin());
  }
  std::shared_ptr<Context> context() const { return context_; }
 private:
  std::shared_ptr<Context> context_;
};

class Exact {
 public:
  struct Term { Rational coefficient; std::vector<unsigned long> powers; };
  explicit Exact(const ExactField& field, long value = 0) : Exact(field.context()) {
    fmpz_mpoly_q_set_si(value_, value, ctx_->raw);
  }
  Exact(const ExactField& field, const std::string& expression) : Exact(field.context()) {
    // Check division and exponent domains before reaching FLINT routines that
    // assume valid field operands. Data syntax is never executed as code.
    auto parsed=parse(data::Reader(expression).read());
    swap(parsed);
  }
  Exact(const Exact& other) : Exact(other.ctx_) {
    fmpz_mpoly_q_set(value_, other.value_, ctx_->raw);
  }
  Exact(Exact&& other) noexcept : Exact(other.ctx_) {
    fmpz_mpoly_q_swap(value_, other.value_, ctx_->raw);
  }
  ~Exact() { fmpz_mpoly_q_clear(value_, ctx_->raw); }
  Exact& operator=(Exact other) noexcept { swap(other); return *this; }
  void swap(Exact& other) noexcept {
    // Each raw value must travel with its own context, even across fields.
    std::swap(value_[0], other.value_[0]);
    ctx_.swap(other.ctx_);
  }
  bool is_zero() const { return fmpz_mpoly_q_is_zero(value_, ctx_->raw); }
  bool is_rational() const { return fmpz_mpoly_q_is_fmpq(value_, ctx_->raw); }
  std::string str() const {
    char* p = fmpz_mpoly_q_get_str_pretty(value_, ctx_->symbols.data(), ctx_->raw);
    if (!p) throw std::bad_alloc();
    std::string result(p); flint_free(p); return result;
  }
  Rational rational() const {
    if (!is_rational()) throw std::domain_error("expression is not constant");
    auto n = numerator_terms(), d = denominator_terms();
    return n.empty() ? Rational(0) : n.at(0).coefficient / d.at(0).coefficient;
  }
  Exact constant(long n) const {
    Exact out(ctx_); fmpz_mpoly_q_set_si(out.value_, n, ctx_->raw); return out;
  }
  Exact constant(const Rational& n) const {
    Exact out(ctx_);
    fmpq_t q; fmpq_init(q);
    fmpq_set_str(q,n.str().c_str(),10); fmpq_canonicalise(q);
    fmpz_mpoly_q_set_fmpq(out.value_,q,ctx_->raw); fmpq_clear(q);
    return out;
  }
  // Parse another value in this value's existing field context.
  Exact parse(const std::string& expression) const {
    return parse(data::Reader(expression).read());
  }
  std::size_t variable_count() const { return ctx_->names.size(); }
  const std::vector<std::string>& variables() const {return ctx_->names;}
  Exact polynomial_lcm(const Exact& b) const {
    same_field(b);
    if (!denominator().is_one() || !b.denominator().is_one())
      throw std::invalid_argument("polynomial lcm requires integer polynomials");
    if (is_zero() || b.is_zero()) return constant(0);
    Exact g(ctx_);
    if (!fmpz_mpoly_gcd(fmpz_mpoly_q_numref(g.value_), fmpz_mpoly_q_numref(value_),
                        fmpz_mpoly_q_numref(b.value_),ctx_->raw))
      throw std::runtime_error("FLINT polynomial gcd failed");
    return (*this/g)*b;
  }
  bool is_one() const { return fmpz_mpoly_q_is_one(value_,ctx_->raw); }
  Exact variable(std::size_t i) const {
    if (i >= ctx_->names.size()) throw std::out_of_range("exact variable index");
    Exact out(ctx_); fmpz_mpoly_q_gen(out.value_, i, ctx_->raw); return out;
  }
  Exact derivative(std::size_t i) const {
    if (i >= ctx_->names.size()) throw std::out_of_range("derivative variable index");
    Exact n = numerator(), d = denominator(), dn(ctx_), dd(ctx_);
    fmpz_mpoly_derivative(fmpz_mpoly_q_numref(dn.value_), fmpz_mpoly_q_numref(value_), i, ctx_->raw);
    fmpz_mpoly_derivative(fmpz_mpoly_q_numref(dd.value_), fmpz_mpoly_q_denref(value_), i, ctx_->raw);
    return (dn*d-n*dd)/(d*d);
  }
  Exact numerator() const {
    Exact out(ctx_);
    fmpz_mpoly_set(fmpz_mpoly_q_numref(out.value_), fmpz_mpoly_q_numref(value_), ctx_->raw);
    return out;
  }
  Exact denominator() const {
    Exact out(ctx_);
    fmpz_mpoly_set(fmpz_mpoly_q_numref(out.value_), fmpz_mpoly_q_denref(value_), ctx_->raw);
    return out;
  }
  std::vector<Term> numerator_terms() const { return terms(fmpz_mpoly_q_numref(value_)); }
  std::vector<Term> denominator_terms() const { return terms(fmpz_mpoly_q_denref(value_)); }
  void univariate_polynomials(fmpz_poly_t numerator,fmpz_poly_t denominator,std::size_t variable)const {
    if(variable>=variable_count())throw std::out_of_range("univariate coefficient variable");
    if(!fmpz_mpoly_is_fmpz_poly(fmpz_mpoly_q_numref(value_),variable,ctx_->raw)||
       !fmpz_mpoly_is_fmpz_poly(fmpz_mpoly_q_denref(value_),variable,ctx_->raw))
      throw std::domain_error("coefficient depends on more than the selected variable");
    if(!fmpz_mpoly_get_fmpz_poly(numerator,fmpz_mpoly_q_numref(value_),variable,ctx_->raw)||
       !fmpz_mpoly_get_fmpz_poly(denominator,fmpz_mpoly_q_denref(value_),variable,ctx_->raw))
      throw std::domain_error("coefficient depends on more than the selected variable");
  }
  Exact from_univariate_polynomials(const fmpz_poly_t numerator,const fmpz_poly_t denominator,std::size_t variable)const {
    if(variable>=variable_count())throw std::out_of_range("univariate coefficient variable");
    if(fmpz_poly_is_zero(denominator))throw std::domain_error("zero univariate denominator");
    Exact out(ctx_);
    fmpz_mpoly_set_fmpz_poly(fmpz_mpoly_q_numref(out.value_),numerator,variable,ctx_->raw);
    fmpz_mpoly_set_fmpz_poly(fmpz_mpoly_q_denref(out.value_),denominator,variable,ctx_->raw);
    fmpz_mpoly_q_canonicalise(out.value_,ctx_->raw);return out;
  }
  Exact pow(unsigned long n) const {
    Exact a(*this), out = constant(1);
    while (n) { if (n & 1) out = out*a; n >>= 1; if (n) a = a*a; }
    return out;
  }
  // Simultaneous, exact substitution, including rational replacements. It
  // intentionally does not reinterpret a symbol in a different field.
  Exact substitute(std::span<const Exact> replacements) const {
    if (replacements.size() != ctx_->names.size())
      throw std::invalid_argument("substitution has wrong number of variables");
    for (const auto& x : replacements) same_field(x);
    const auto evaluate = [&](const std::vector<Term>& ts) {
      Exact sum = constant(0);
      for (const auto& term : ts) {
        Exact t(ctx_);
        // Polynomial coefficients are integers, not machine-sized integers.
        fmpz_t c; fmpz_init(c);
        fmpz_set_str(c, term.coefficient.str().c_str(), 10);
        fmpz_mpoly_q_set_fmpz(t.value_, c, ctx_->raw); fmpz_clear(c);
        for (std::size_t i=0; i<replacements.size(); ++i)
          if (term.powers[i]) t = t*replacements[i].pow(term.powers[i]);
        sum = sum+t;
      }
      return sum;
    };
    return evaluate(numerator_terms())/evaluate(denominator_terms());
  }
  friend Exact operator+(const Exact& a, const Exact& b) {
    a.same_field(b); Exact c(a.ctx_); fmpz_mpoly_q_add(c.value_,a.value_,b.value_,a.ctx_->raw); return c;
  }
  friend Exact operator-(const Exact& a, const Exact& b) {
    a.same_field(b); Exact c(a.ctx_); fmpz_mpoly_q_sub(c.value_,a.value_,b.value_,a.ctx_->raw); return c;
  }
  friend Exact operator-(const Exact& a) {
    Exact c(a.ctx_); fmpz_mpoly_q_neg(c.value_,a.value_,a.ctx_->raw); return c;
  }
  friend Exact operator*(const Exact& a, const Exact& b) {
    a.same_field(b); Exact c(a.ctx_); fmpz_mpoly_q_mul(c.value_,a.value_,b.value_,a.ctx_->raw); return c;
  }
  friend Exact operator/(const Exact& a, const Exact& b) {
    a.same_field(b);
    if (b.is_zero()) throw std::domain_error("exact division by zero");
    Exact c(a.ctx_); fmpz_mpoly_q_div(c.value_,a.value_,b.value_,a.ctx_->raw); return c;
  }
  friend bool operator==(const Exact& a, const Exact& b) {
    a.same_field(b); return fmpz_mpoly_q_equal(a.value_,b.value_,a.ctx_->raw);
  }
 private:
  explicit Exact(std::shared_ptr<ExactField::Context> context) : ctx_(std::move(context)) {
    fmpz_mpoly_q_init(value_, ctx_->raw);
    fmpz_mpoly_q_zero(value_, ctx_->raw);
  }
  void same_field(const Exact& b) const {
    if (ctx_ != b.ctx_) throw std::invalid_argument("operation mixes distinct exact fields");
  }
  Exact parse(const data::Expr& e) const {
    if(data::number(e)) {
      if(!std::all_of(e.head.begin(),e.head.end(),[](unsigned char c){return std::isdigit(c);}))
        throw std::invalid_argument("exact input requires integers and rational fractions");
      Exact out(ctx_); fmpz_t n;fmpz_init(n);
      const int status=fmpz_set_str(n,e.head.c_str(),10);
      if(status==0)fmpz_mpoly_q_set_fmpz(out.value_,n,ctx_->raw);
      fmpz_clear(n);
      if(status)throw std::invalid_argument("invalid integer");
      return out;
    }
    if(e.atom()) {
      auto i=std::find(ctx_->names.begin(),ctx_->names.end(),e.head);
      if(i==ctx_->names.end())throw std::invalid_argument("undeclared exact variable: "+e.head);
      return variable(i-ctx_->names.begin());
    }
    if(e.head=="neg" && e.args.size()==1)return -parse(e.args[0]);
    if(e.args.size()!=2)throw std::invalid_argument("unsupported exact expression: "+e.head);
    if(e.head=="^") {
      long n=data::integer(e.args[1]);
      if(n>1000000 || n< -1000000)throw std::invalid_argument("exact exponent exceeds resource limit");
      auto a=parse(e.args[0]);return n<0?constant(1)/a.pow(-n):a.pow(n);
    }
    auto a=parse(e.args[0]),b=parse(e.args[1]);
    if(e.head=="+")return a+b;
    if(e.head=="-")return a-b;
    if(e.head=="*")return a*b;
    if(e.head=="/")return a/b;
    throw std::invalid_argument("unsupported exact operator: "+e.head);
  }
  std::vector<Term> terms(const fmpz_mpoly_t p) const {
    std::vector<Term> out;
    if (!fmpz_mpoly_total_degree_fits_si(p, ctx_->raw))
      throw std::overflow_error("polynomial exponents exceed machine range");
    for (slong i=0; i<fmpz_mpoly_length(p,ctx_->raw); ++i) {
      fmpz_t c; fmpz_init(c); fmpz_mpoly_get_term_coeff_fmpz(c,p,i,ctx_->raw);
      char* s=fmpz_get_str(nullptr,10,c); fmpz_clear(c);
      if (!s) throw std::bad_alloc();
      std::string coefficient(s); flint_free(s);
      Term t{Rational(coefficient), std::vector<unsigned long>(ctx_->names.size())};
      fmpz_mpoly_get_term_exp_ui(t.powers.data(),p,i,ctx_->raw);
      out.push_back(std::move(t));
    }
    return out;
  }
  std::shared_ptr<ExactField::Context> ctx_;
  fmpz_mpoly_q_t value_;
};
}  // namespace diffexp

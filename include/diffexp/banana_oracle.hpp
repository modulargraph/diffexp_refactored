#pragma once

#include "diffexp/kernel/scalar.hpp"
#include <flint/acb_calc.h>
#include <flint/acb_hypgeom.h>
#include <flint/arb_calc.h>
#include <array>
#include <algorithm>

namespace diffexp::oracle {
using Ball = kernel::ComplexBall;
using Rational = kernel::Rational;

// Independent coordinate-space oracle, not a FeynmanTrick solution. Fixed
// D=2, epsilon^0, p^2=-1, two through five unit propagator powers. Positive square roots
// of squared masses and principal K0 on the positive real axis are used.
// There is no physical cut here; upper/lower Feynman rims agree.
struct BananaOptions {
  slong target_bits = 80; // Required absolute radius <= 2^-target_bits.
  slong working_bits = 192;
  slong max_evaluations = 200000; // Total callback budget, across all pieces.
  slong max_depth = 40;
};
struct BananaResult {
  Ball value;
  Ball lower_tail_bound;
  Ball upper_tail_bound;
  slong evaluations = 0;
  slong analytic_queries = 0;
  slong target_bits = 0;
  slong working_bits = 0;
  long logarithmic_lower_cutoff = 0;
  long radial_upper_cutoff = 0;
  unsigned propagators = 0;
  // Cutoffs refer to t=m_min*r. For the original Banana4 fixtures m_min=1.
  Rational minimum_squared_mass{1};
};
using Banana4Options=BananaOptions;
using Banana4Result=BananaResult;
namespace detail {
struct BananaContext {
  std::vector<Ball> masses=std::vector<Ball>(5);
  Ball momentum_scale{1}, jacobian{1};
  slong evaluations = 0, analytic_queries = 0, max_evaluations = 0;
};
inline int banana_density(acb_ptr out, const acb_t u, void* opaque,
                          slong order, slong bits) {
  auto& context = *static_cast<BananaContext*>(opaque);
  ++context.evaluations;
  if (order == 1) ++context.analytic_queries;
  if (order > 1 || context.evaluations > context.max_evaluations) {
    for (slong i=0; i < (order > 1 ? order : 1); ++i) acb_indeterminate(out+i);
    return 0;
  }
  Ball r, zero, j, k, argument;
  acb_exp(r.raw(),u,bits);
  // An order=1 query asks for a holomorphic enclosure on the entire complex
  // input rectangle, NOT a derivative or midpoint evaluation. Re(exp u)>0
  // guarantees every mass*exp(u) avoids K0's cut and zero. Acb evaluates the
  // whole ball, including its imaginary radius, for both kinds of query.
  if (!arb_is_positive(acb_realref(r.raw()))) {
    acb_indeterminate(out); return 0;
  }
  acb_mul(argument.raw(),r.raw(),context.momentum_scale.raw(),bits);
  acb_hypgeom_bessel_j(j.raw(),zero.raw(),argument.raw(),bits);
  acb_mul(out,r.raw(),r.raw(),bits);
  acb_mul(out,out,j.raw(),bits);
  acb_mul_2exp_si(out,out,static_cast<slong>(context.masses.size()-1));
  acb_mul(out,out,context.jacobian.raw(),bits);
  for (const auto& mass:context.masses) {
    acb_mul(argument.raw(),mass.raw(),r.raw(),bits);
    acb_hypgeom_bessel_k(k.raw(),zero.raw(),argument.raw(),bits);
    acb_mul(out,out,k.raw(),bits);
  }
  return 0;
}
} // namespace detail

inline BananaResult banana_bessel(
    const std::vector<Rational>& squared_masses,
    const BananaOptions& options = {}) {
  if (options.target_bits < 24 || options.target_bits > 512 ||
      options.working_bits < options.target_bits+48 || options.working_bits > 4096 ||
      options.max_evaluations < 1 || options.max_evaluations > 10000000 ||
      options.max_depth < 1 || options.max_depth > 100)
    throw std::invalid_argument("invalid banana oracle accuracy/resource limits");
  const ulong N=squared_masses.size();
  if(N<2 || N>5)throw std::invalid_argument("banana oracle supports one through four loops");
  Rational minimum=squared_masses.front();
  for(const auto& m:squared_masses) {
    if(m.sign()<=0)throw std::invalid_argument("banana oracle requires positive squared masses");
    if(m<minimum)minimum=m;
  }
  const auto bits=options.working_bits;
  detail::BananaContext context;
  context.max_evaluations=options.max_evaluations;context.masses.resize(N);
  // Rescale t=m_min*r so the tail proof sees mass ratios >=1, even when
  // input masses are below one. J0(t/m_min) remains bounded by one for real t.
  fmpq_t minimum_q;fmpq_init(minimum_q);fmpq_set_str(minimum_q,minimum.str().c_str(),10);
  const slong scale_bits=std::max<slong>(0,static_cast<slong>(fmpz_bits(fmpq_denref(minimum_q)))-
      static_cast<slong>(fmpz_bits(fmpq_numref(minimum_q)))+1);
  if(scale_bits>2048) {fmpq_clear(minimum_q);throw std::invalid_argument("banana oracle mass rescaling exceeds resource limit");}
  arb_set_fmpq(acb_realref(context.jacobian.raw()),minimum_q,bits);fmpq_clear(minimum_q);
  acb_inv(context.jacobian.raw(),context.jacobian.raw(),bits);
  acb_sqrt(context.momentum_scale.raw(),context.jacobian.raw(),bits);
  for (std::size_t i=0;i<N;++i) {
    fmpq_t q; fmpq_init(q);fmpq_set_str(q,(squared_masses[i]/minimum).str().c_str(),10);
    arb_set_fmpq(acb_realref(context.masses[i].raw()),q,bits);fmpq_clear(q);
    acb_sqrt(context.masses[i].raw(),context.masses[i].raw(),bits);
  }
  BananaResult result;
  result.propagators=N;result.minimum_squared_mass=minimum;
  result.target_bits=options.target_bits;result.working_bits=bits;
  const long A=options.target_bits+scale_bits+32, R=(options.target_bits+scale_bits)/N+16;
  result.logarithmic_lower_cutoff=-A; result.radial_upper_cutoff=R;
  Ball a(-A), L(2), term, polynomial;
  acb_exp(a.raw(),a.raw(),bits);
  acb_log(L.raw(),L.raw(),bits);
  acb_add_ui(L.raw(),L.raw(),A+1,bits);
  // K0(x)=integral exp(-x cosh t) dt <= E1(x/2)
  // <= log(2/x)+1 for 0<x<=1 (split E1 at 1).
  // K0(mr)<=K0(r) for m>=1 and |J0(r)|<=1 on the real axis.
  // Therefore low tail <= 2^(N-1)/m_min^2 * a^2 sum_{k=0}^N
  //   [N!/(N-k)!] (log(2/a)+1)^(N-k)/2^(k+1), a=exp(-A).
  ulong falling=1;
  for (ulong k=0;k<=N;++k) {
    acb_pow_ui(term.raw(),L.raw(),N-k,bits);
    acb_mul_ui(term.raw(),term.raw(),falling,bits);
    acb_mul_2exp_si(term.raw(),term.raw(),-static_cast<slong>(k+1));
    acb_add(polynomial.raw(),polynomial.raw(),term.raw(),bits);
    falling*=N-k;
  }
  acb_mul(result.lower_tail_bound.raw(),a.raw(),a.raw(),bits);
  acb_mul(result.lower_tail_bound.raw(),result.lower_tail_bound.raw(),polynomial.raw(),bits);
  acb_mul_2exp_si(result.lower_tail_bound.raw(),result.lower_tail_bound.raw(),N-1);
  acb_mul(result.lower_tail_bound.raw(),result.lower_tail_bound.raw(),context.jacobian.raw(),bits);
  // cosh(t)>=1+t^2/2 gives K0(x)<=sqrt(pi/(2x)) exp(-x).
  // For R>=1 and all masses>=1, the high tail is bounded by
  // 2^(N-1)(pi/2)^(N/2)/m_min^2 * integral_R^infty t^(1-N/2)exp(-Nt)dt
  // <= 2^(N-1)*pi^N/m_min^2 * exp(-N R), for N>=2.
  Ball pi, decay(-static_cast<long>(N)*R);
  arb_const_pi(acb_realref(pi.raw()),bits);
  acb_pow_ui(pi.raw(),pi.raw(),N,bits);
  acb_exp(decay.raw(),decay.raw(),bits);
  acb_mul(result.upper_tail_bound.raw(),pi.raw(),decay.raw(),bits);
  acb_mul_2exp_si(result.upper_tail_bound.raw(),result.upper_tail_bound.raw(),N-1);
  acb_mul(result.upper_tail_bound.raw(),result.upper_tail_bound.raw(),context.jacobian.raw(),bits);

  acb_calc_integrate_opt_t integration; acb_calc_integrate_opt_init(integration);
  integration->depth_limit=options.max_depth;
  integration->deg_limit=120;
  mag_t tolerance; mag_init(tolerance);
  mag_set_ui_2exp_si(tolerance,1,-options.target_bits-20);
  Ball endpoint(R), left(-A), right, piece;
  acb_log(endpoint.raw(),endpoint.raw(),bits);
  // Unit intervals keep the analytic test rectangles inside a useful strip.
  for (long j=-A;j<=R;++j) {
    const bool last=arb_le(acb_realref(endpoint.raw()),acb_realref(Ball(j+1).raw()));
    if(last) acb_set(right.raw(),endpoint.raw()); else acb_set_si(right.raw(),j+1);
    integration->eval_limit=options.max_evaluations-context.evaluations;
    if(integration->eval_limit<1) {
      mag_clear(tolerance);throw std::runtime_error("banana oracle exhausted total evaluation budget");
    }
    const int status=acb_calc_integrate(piece.raw(),detail::banana_density,&context,
        left.raw(),right.raw(),options.target_bits+20,tolerance,integration,bits);
    if(status!=ARB_CALC_SUCCESS || !piece.is_finite() || context.evaluations>options.max_evaluations) {
      mag_clear(tolerance);throw std::runtime_error("banana oracle finite-interval integration did not converge within limits");
    }
    acb_add(result.value.raw(),result.value.raw(),piece.raw(),bits);
    if(last) break;
    acb_set(left.raw(),right.raw());
  }
  mag_t tail;mag_init(tail);
  arb_get_mag(tail,acb_realref(result.lower_tail_bound.raw()));
  arb_add_error_mag(acb_realref(result.value.raw()),tail);
  arb_get_mag(tail,acb_realref(result.upper_tail_bound.raw()));
  arb_add_error_mag(acb_realref(result.value.raw()),tail);
  // The exact real-axis integral is real; intersecting its enclosure with the
  // real axis is valid only after checking that zero is in its imaginary ball.
  const bool real=arb_contains_zero(acb_imagref(result.value.raw()));
  if(real) arb_zero(acb_imagref(result.value.raw()));
  mag_set_ui_2exp_si(tolerance,1,-options.target_bits);
  const bool accurate=result.value.is_finite() && real &&
      mag_cmp(arb_radref(acb_realref(result.value.raw())),tolerance)<=0;
  mag_clear(tail);mag_clear(tolerance);
  if(!accurate)throw std::runtime_error("banana oracle failed requested absolute radius");
  result.evaluations=context.evaluations;result.analytic_queries=context.analytic_queries;
  return result;
}
inline BananaResult equal_banana_bessel(unsigned loops,const BananaOptions& options = {}) {
  if(loops<1 || loops>4)throw std::invalid_argument("banana oracle supports one through four loops");
  return banana_bessel(std::vector<Rational>(loops+1,Rational(1)),options);
}
// Preserve the original Banana4 API and its >=1 mass-domain contract.
inline Banana4Result banana4_bessel(const std::array<Rational,5>& squared_masses,
                                    const Banana4Options& options = {}) {
  for(const auto& m:squared_masses)if(m<Rational(1))
    throw std::invalid_argument("banana4 oracle wrapper requires squared masses >= 1; use banana_bessel for rescaling");
  return banana_bessel(std::vector<Rational>(squared_masses.begin(),squared_masses.end()),options);
}
inline Banana4Result banana4_equal_bessel(const Banana4Options& options = {}) {
  return banana4_bessel({Rational(1),Rational(1),Rational(1),Rational(1),Rational(1)},options);
}
inline Banana4Result banana4_unequal_bessel(const Banana4Options& options = {}) {
  return banana4_bessel({Rational(2),Rational("3/2"),Rational("4/3"),Rational("5/4"),Rational(1)},options);
}
} // namespace diffexp::oracle

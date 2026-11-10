#pragma once
#include "diffexp/adjoint_transport.hpp"
#include "diffexp/diagonal_log_gauge.hpp"
#include "diffexp/transport.hpp"
#include "diffexp/ft_ultraspherical_detail.hpp"
#include "diffexp/ft_spectral.hpp"


namespace diffexp::ft_ultraspherical {
using B = Jet::Ball;
namespace ultra = ft_ultraspherical_detail;
namespace sp = transport::spectral_detail;
// Explicit opt-in: raw affine paths and internal precision protection.
struct Options : ft_spectral::Options {
  Options(){diagonal_gauge=false;endpoint_clustering=false;}
  slong working_bits=512;
  // Original leg index and completed fraction of that leg, after acceptance.
  std::function<void(unsigned,double)> progress;
};
using Diagnostics=ft_spectral::Diagnostics;
// Solve g'=forcing-g*A in Chebyshev coefficient space. The leg contract,
// epsilon-zero block order, singularity geometry, refinement test and
// subdivision are those of ft_spectral; only the solve of one level differs.
inline std::optional<LaurentRows>
try_transport_single(const ExactEpsilonMatrix &matrix, const LaurentRows &initial,
              const ExactEpsilonMatrix &forcing,
              const std::vector<Exact> &vertices, const Options &options,
              Diagnostics &diagnostics, unsigned arithmetic_accuracy_goal = 0) {
  using namespace ft_spectral::detail;
  diagnostics = {};
  auto started = std::chrono::steady_clock::now();
  double preparing = 0;
  auto elapsed = [&] {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         started)
        .count();
  };
  auto check = [&] {
    if (elapsed() > options.seconds_budget)
      throw sp::Rejected("FT spectral time budget exhausted");
  };
  kernel::PrecisionScope precision_restore;
  try {
    if(options.working_bits<128||options.working_bits>16384)throw sp::Rejected("FT ultraspherical internal precision budget");
    B::set_precision(std::max<slong>({512,options.working_bits,precision_restore.previous}));
    const auto bits = B::precision();
    check_request(matrix, initial, forcing, vertices, options, bits);
    const auto d = initial.columns(), r = initial.coefficients.size();
    // A conservative routing policy, not a proof that the dynamics cannot
    // contract uncertainty. Decline this backend before expensive preparation;
    // callers may use a different transport with the original input intact.
    const auto input_tolerance=B::from_strings("1e-"+std::to_string(
        (arithmetic_accuracy_goal?arithmetic_accuracy_goal:options.accuracy_goal)+2));
    for(const auto& row:initial.coefficients)for(const auto& series:row)for(const auto& value:series){
      if(!value.is_finite())throw sp::Rejected("FT spectral nonfinite boundary");
      if(!sp::le(transport::arithmetic_error(value),input_tolerance*sp::midpoint_scale(value)))
        throw sp::Rejected("FT ultraspherical inherited input exceeds arm budget; use fallback");
    }
    auto window = open_window(matrix, initial, forcing, vertices, options,
                              arithmetic_accuracy_goal);
    const auto xi = window.xi, ei = window.ei;
    const int low = window.low;
    const unsigned w = window.w;
    auto &current = window.current;
    for(const auto& row:matrix)for(const auto& entry:row)
      if(!entry.is_zero()&&*exact_epsilon_valuation(entry,ei)<0)
        throw sp::Rejected("FT spectral connection needs epsilon normalization");
    for (unsigned leg = 0; leg + 1 < vertices.size(); ++leg) {
      const auto prep = std::chrono::steady_clock::now();
      check();
      const auto &from = vertices[leg];
      const auto width = vertices[leg + 1] - from;
      if (width.is_zero())
        continue;
      // Coefficient space works on the raw affine path: there is no
      // parameter clustering and no diagonal gauge to apply.
      if(options.endpoint_clustering||options.diagonal_gauge)throw sp::Rejected("coefficient prototype requires affine ungauged path");
      auto point = exact_point(from, xi, from + width * from.variable(xi));
      const auto coupling = couple(
          matrix, forcing, point, width, window, options.max_nodes,
          std::vector<const Exact *>(d, nullptr),
          [](const Exact &, unsigned, unsigned) {
            return std::pair<unsigned, bool>{0, true};
          });
      unsigned minimum_degree = coupling.minimum_degree;
      std::vector<ultra::EP> full_q;
      std::vector<std::vector<ultra::EP>> full_p(d),full_r(d);
      auto chebpoint=exact_point(from,xi,from+width*(from.variable(xi)+from.constant(1))/from.constant(2));
      unsigned polynomial_degree=0;
      for(unsigned i=0;i<d;++i){auto den=from.constant(1);for(unsigned j=0;j<d;++j)den=den.polynomial_lcm(matrix[j][i].denominator());
        std::vector<Exact> normalized;for(unsigned a=0;a<r;++a){auto v=forcing[a][i];v=low<0?v*v.variable(ei).pow(-low):v/v.variable(ei).pow(low);normalized.push_back(v);den=den.polynomial_lcm(v.denominator());}
        unsigned epsdegree=0;for(auto&t:den.numerator_terms())epsdegree=std::max(epsdegree,unsigned(t.powers[ei]));
        auto denseries=feynman::scalar_functional_detail::epsilon_series(den,ei,epsdegree);bool separable=!denseries[0].is_zero();for(auto&v:denseries)if(separable&&!(v/denseries[0]).is_rational())separable=false;
        auto compile_series=[&](const Exact&v){ultra::EP result;auto series=feynman::scalar_functional_detail::epsilon_series(v,ei,w-1);for(auto&c:series)result.push_back(ultra::coeff(c.substitute(chebpoint)));while(result.size()>1&&std::all_of(result.back().begin(),result.back().end(),[](const B&v){return v.is_zero();}))result.pop_back();return result;};
        if(separable){const auto&den0=denseries[0];full_q.push_back({ultra::coeff(den0.substitute(chebpoint))});
          for(unsigned j=0;j<d;++j)full_p[i].push_back(compile_series(-width/from.constant(2)*den0*matrix[j][i]));
          for(unsigned a=0;a<r;++a)full_r[i].push_back(compile_series(width/from.constant(2)*den0*normalized[a]));
        }else{full_q.push_back(ultra::bivariate(den.substitute(chebpoint)));
          for(unsigned j=0;j<d;++j)full_p[i].push_back(ultra::bivariate((-width/from.constant(2)*den*matrix[j][i]).substitute(chebpoint)));
          for(unsigned a=0;a<r;++a)full_r[i].push_back(ultra::bivariate((width/from.constant(2)*den*normalized[a]).substitute(chebpoint)));
        }
        for(auto&p:full_q.back())polynomial_degree=std::max(polynomial_degree,unsigned(p.size()-1));
        for(auto&entry:full_p[i])for(auto&p:entry)polynomial_degree=std::max(polynomial_degree,unsigned(p.size()-1));
        for(auto&entry:full_r[i])for(auto&p:entry)polynomial_degree=std::max(polynomial_degree,unsigned(p.size()-1));
      }
      if(polynomial_degree>32)throw sp::Rejected("FT ultraspherical polynomial degree budget");
      minimum_degree=std::max(minimum_degree,polynomial_degree);
      const auto partitioned = partition(d, coupling.edges, options, diagnostics);
      const auto &blocks = partitioned.members;
      const auto &block_id = partitioned.block_id;
      const auto geometry = path_geometry(coupling.denominators, xi, bits, check);
      if (sp::node_forecast(geometry.singularities, options.accuracy_goal) >
          (options.conservative ? 2 * options.max_nodes : 512))
        throw sp::Rejected("FT spectral interval too close to singularities");
      preparing +=
          std::chrono::duration<double>(std::chrono::steady_clock::now() - prep)
              .count();
      auto midpoint=current, delta=current;
      std::vector<unsigned> active;
      for(unsigned j=0;j<d;++j){bool uncertain=false;
        for(unsigned a=0;a<r;++a)for(unsigned e=0;e<w;++e){auto&v=current.coefficients[a][j][e];acb_get_mid(midpoint.coefficients[a][j][e].raw(),v.raw());delta.coefficients[a][j][e]=v-midpoint.coefficients[a][j][e];uncertain|=!delta.coefficients[a][j][e].is_zero();}
        if(uncertain&&!std::all_of(matrix[j].begin(),matrix[j].end(),[](const Exact&v){return v.is_zero();}))active.push_back(j);
      }
      const unsigned nr=r+active.size();
      Refinement refinement;
      std::vector<B> previous_affine;
      bool accepted = false;
      for (unsigned n : {8u,12u,16u,24u,32u,48u,64u,80u,96u,128u}) {
        if (n > options.max_nodes)
          break;
        check();
        // Same interval-count budget and retry contract as the collocation
        // backend: beyond the first level this only ends refinement.
        if (partitioned.largest * n > options.max_block_nodes) {
          if (!refinement.previous_n)
            throw sp::Rejected("FT spectral coupled-node budget");
          break;
        }
        const unsigned m = n + 1, lim = n + polynomial_degree + 2;
        // Two coefficient arrays (solution and derivative) of lim rows each.
        if(static_cast<std::uint64_t>(2)*lim*w*d*nr>options.max_cells)throw sp::Rejected("FT spectral solution storage budget");
        diagnostics.nodes.push_back(n);
        ultra::Multipliers multipliers(options.max_cells);
        // Only endpoint values are compared between levels.
        sp::Matrix state(w*d,nr);
        std::vector<std::shared_ptr<ultra::Factor>> factors;
        std::size_t factor_cells=0;
        for(unsigned b=0;b<blocks.size();++b){check();unsigned size=blocks[b].size();std::vector<ultra::V> q;std::vector<std::vector<ultra::V>> p;
          for(auto i:blocks[b]){q.push_back(full_q[i][0]);p.emplace_back();for(auto j:blocks[b])p.back().push_back(full_p[i][j][0]);}
          factor_cells+=static_cast<std::size_t>(size)*size*m*m;if(factor_cells>options.max_cells)throw sp::Rejected("FT ultraspherical factor storage budget");
          factors.push_back(std::make_shared<ultra::Factor>(ultra::build(q,p,n)));++diagnostics.factorizations;
        }
        std::vector<ultra::M> ustate(w*d),dstate(w*d);
        for(unsigned e=0;e<w;++e)for(unsigned b=0;b<blocks.size();++b){check();unsigned size=blocks[b].size();ultra::M rhs(size*m,ultra::V(nr));
          for(unsigned i=0;i<size;++i){unsigned row=blocks[b][i];
            for(unsigned a=0;a<r;++a){rhs[n*size+i][a]=midpoint.coefficients[a][row][e];if(e<full_r[row][a].size()){ultra::V unit(lim);unit[0]=B(1);auto forcingpoly=ultra::polyapply(full_r[row][a][e],unit);for(unsigned k=0;k<n;++k)rhs[k*size+i][a]+=forcingpoly[k];}}
            if(e==0)for(unsigned a=0;a<active.size();++a)if(active[a]==row)rhs[n*size+i][r+a]=B(1);
            for(unsigned j=0;j<d;++j)for(unsigned lag=0;lag<=e&&lag<full_p[row][j].size();++lag){
              if(lag==0&&block_id[j]==b)continue;
              const auto&p=full_p[row][j][lag];if(std::all_of(p.begin(),p.end(),[](const B&v){return v.is_zero();}))continue;
              if(ustate[(e-lag)*d+j].empty())throw std::logic_error("ultra epsilon/SCC dependency order");multipliers.add_product(rhs,i,size,p,ustate[(e-lag)*d+j],false);
            }
            for(unsigned lag=1;lag<=e&&lag<full_q[row].size();++lag)multipliers.add_product(rhs,i,size,full_q[row][lag],dstate[(e-lag)*d+row],true);
          }
          auto solved=factors[b]->solve(std::move(rhs));
          for(unsigned i=0;i<size;++i){unsigned row=blocks[b][i];auto&u=ustate[e*d+row];auto&deriv=dstate[e*d+row];u=ultra::M(lim,ultra::V(nr));deriv=ultra::M(lim,ultra::V(nr));
            for(unsigned a=0;a<nr;++a){B endpoint;for(unsigned k=0;k<=n;++k){const auto&v=solved[k?(k-1)*size+i:n*size+i][a];endpoint+=v;if(k==0)u[0][a]+=v;else{u[k][a]+=v/B(2);if(k>=2)u[k-2][a]-=v/B(2);deriv[k-1][a]=v*B(k);}}acb_set(acb_mat_entry(state.value,e*d+row,a),endpoint.raw());}
          }
        }
        // Keep finite affine operators separate between resolutions. The
        // midpoint of U*delta is zero for symmetric input noise, so comparing
        // only completed output midpoints would hide transfer discretization.
        std::vector<B> affine(nr*d*w);
        for(unsigned a=0;a<nr;++a)for(unsigned i=0;i<d;++i)for(unsigned e=0;e<w;++e)
          acb_set(affine[(a*d+i)*w+e].raw(),acb_mat_entry(state.value,e*d+i,a));
        std::vector<B> values(r * d * w), difference(r * d * w);
        for (unsigned a = 0; a < r; ++a)
          for (unsigned i = 0; i < d; ++i)
            for (unsigned e = 0; e < w; ++e) {
              auto id = (a * d + i) * w + e;
              acb_set(values[id].raw(),
                      acb_mat_entry(state.value, e * d + i, a));
              for(unsigned j=0;j<d;++j){auto it=std::find(active.begin(),active.end(),j);
                if(it==active.end()){if(i==j)values[id]+=delta.coefficients[a][j][e];continue;}
                const unsigned basis=r+std::distance(active.begin(),it);
                for(unsigned lag=0;lag<=e;++lag){B u;acb_set(u.raw(),acb_mat_entry(state.value,lag*d+i,basis));values[id]+=u*delta.coefficients[a][j][e-lag];}
              }
              if (!values[id].is_finite())
                throw sp::Rejected("FT spectral nonfinite endpoint");
              if (!previous_affine.empty()) {
                difference[id]=sp::difference(affine[id],previous_affine[id]);
                for(unsigned j=0;j<active.size();++j)for(unsigned lag=0;lag<=e;++lag){
                  const auto map_id=((r+j)*d+i)*w+lag;
                  // Full ball difference encloses the change of the finite
                  // transfer maps, including their arithmetic uncertainty.
                  const auto map_change=sp::upper(transport::magnitude(affine[map_id]-previous_affine[map_id]));
                  const auto noise=sp::upper(transport::magnitude(delta.coefficients[a][active[j]][e-lag]));
                  difference[id]+=map_change*noise;
                }
                difference[id]=sp::upper(difference[id]);
              }
            }
        if (!refinement.previous_difference.empty() && n > minimum_degree) {
          unsigned absolute = 0;
          if (auto tails = refinement.converged(values, difference, n, window,
                                                absolute)) {
            commit(current, values, *tails);
            accepted = true;
            diagnostics.absolute_stability_components += absolute;
            ++diagnostics.legs;
            break;
          }
        }
        refinement.advance(std::move(values), std::move(difference), n);
        previous_affine = std::move(affine);
      }
      if (!accepted)
        throw sp::Rejected(
            "FT spectral resolution or arithmetic accuracy exhausted");
      if(options.progress)options.progress(leg,1.);
    }
    diagnostics.preparation_seconds = preparing;
    diagnostics.numerical_seconds = elapsed() - preparing;
    return current;
  } catch (const std::exception &error) {
    diagnostics.reason = error.what();
    diagnostics.preparation_seconds = preparing;
    diagnostics.numerical_seconds = elapsed() - preparing;
    return std::nullopt;
  }
}
// Completed-arm transport with optional bisection of legs; see
// ft_spectral::detail::subdivide for the retry contract. Segment attempts do
// not report progress themselves: the driver reports each accepted fraction.
inline std::optional<LaurentRows>
try_transport(const ExactEpsilonMatrix& matrix,const LaurentRows& initial,
              const ExactEpsilonMatrix& forcing,const std::vector<Exact>& vertices,
              const Options& options,Diagnostics& diagnostics) {
  return ft_spectral::detail::subdivide(initial,vertices,options,diagnostics,
    [&](const LaurentRows& state,const std::vector<Exact>& path,const Options& local,Diagnostics& attempt,unsigned goal) {
      if(!goal)return try_transport_single(matrix,state,forcing,path,local,attempt);
      auto quiet=local;quiet.progress={};
      return try_transport_single(matrix,state,forcing,path,quiet,attempt,goal);
    },[&](unsigned leg,double fraction){if(options.progress)options.progress(leg,fraction);});
}
} // namespace diffexp::ft_ultraspherical

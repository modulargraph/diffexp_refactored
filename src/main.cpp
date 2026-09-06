#include "diffexp/kernel/json_codec.hpp"
#include "diffexp/canonical.hpp"
#include "diffexp/cli.hpp"
#include "diffexp/transport.hpp"
#include "diffexp/planar.hpp"
#include "diffexp/mpl.hpp"
#include "diffexp/rational_transport.hpp"
#include "diffexp/singular_transport.hpp"
#include "diffexp/original_banana_unequal.hpp"
#include "diffexp/bubble.hpp"
#include "diffexp/certified_bubble.hpp"
#include "diffexp/sunrise.hpp"
#include "diffexp/banana_oracle.hpp"
#include "diffexp/direct_examples.hpp"
#if defined(__unix__) || defined(__APPLE__)
#include "diffexp/recursion_graph.hpp"
#include "diffexp/family_config.hpp"
#include "diffexp/recursion_pipeline.hpp"
#include "diffexp/level_cache.hpp"
#include "diffexp/fire_modular.hpp"
#include "diffexp/henn_boundary.hpp"
#endif
#include <iostream>
#include <string>
#include <charconv>
#include <optional>
#include <chrono>

int main(int argc, char** argv) {
  bool json_mode=false;
  try {
    const std::string command = argc > 1 ? argv[1] : "--help";
    if (command == "--version") {
      std::cout << "DiffExp 2.1.1 (standalone C++20 / FLINT)\n";
    } else if(command == "transport") {
      json_mode=true;if(argc!=3)throw std::invalid_argument("transport requires a JSON file or -");
      return diffexp::transport::run_file(argv[2]);
    } else if(command == "family-template") {
#if defined(__unix__) || defined(__APPLE__)
      if(argc!=3)throw std::invalid_argument("family-template requires one example name");
      std::cout<<boost::json::serialize(diffexp::family_config::describe(diffexp::feynman::example_family(argv[2])))<<'\n';
#else
      throw std::runtime_error("Feynman preparation requires a POSIX platform");
#endif
    } else if(command == "prepare" || command == "ft") {
#if defined(__unix__) || defined(__APPLE__)
      for(int i=3;i<argc;++i)if(std::string(argv[i])=="--json")json_mode=true;
      auto& progress=json_mode?std::cerr:std::cout;
      if(argc<3)throw std::invalid_argument(command+" requires a family configuration JSON file");
      auto configuration=diffexp::family_config::load(argv[2]);
      auto options=configuration.preparation;
      auto numerical=configuration.numerical;
      const auto preparation_start=std::chrono::steady_clock::now();
      diffexp::AdjointConditioningStats conditioning;
      numerical.adjoint.conditioning_stats=&conditioning;
      unsigned epsilon_order=configuration.epsilon_order;
      bool numerical_cache=true;
      bool epsilon_order_given=std::string(argv[2])=="-" || std::filesystem::exists(argv[2]);
      auto cache=std::filesystem::temp_directory_path()/"diffexp-exact-cache";
      std::string henn_basis;
      std::filesystem::path fire_prime;
      for(int i=3;i<argc;++i) {
        const std::string option=argv[i];
        if(option=="--json")continue;
        if(command=="ft" && option=="--no-numerical-cache"){numerical_cache=false;continue;}
        if(i+1>=argc)throw std::invalid_argument("missing "+command+" option value");
        const char* option_value=argv[++i];
        if(option=="--fire")options.reduction.provider.executable=option_value;
        else if(option=="--fire-prime")fire_prime=option_value;
        else if(option=="--cache")cache=option_value;
        else if(option=="--henn-basis")henn_basis=option_value;
        else if(option=="--fire-seconds" || option=="--level-seconds" || option=="--total-seconds" ||
            option=="--fire-threads" || option=="--fire-simplifier-threads" || option=="--fire-memory-mib") {
          const std::string value=option_value;unsigned number=0;
          const auto parsed=std::from_chars(value.data(),value.data()+value.size(),number);
          if(parsed.ec!=std::errc{} || parsed.ptr!=value.data()+value.size() || !number || number>86400)
            throw std::invalid_argument("preparation resource budgets must be integers in 1..86400");
          if(option=="--fire-seconds")options.reduction.provider.timeout_seconds=number;
          else if(option=="--level-seconds")options.reduction.total_timeout_seconds=number;
          else if(option=="--total-seconds")options.total_timeout_seconds=number;
          else if(option=="--fire-threads")options.reduction.provider.threads=number;
          else if(option=="--fire-simplifier-threads")options.reduction.provider.simplifier_threads=number;
          else options.reduction.provider.memory_bytes=static_cast<std::size_t>(number)*1024*1024;
        }
        else if(command=="ft" && option=="--ft-transport") {
          const std::string method=option_value;
          if(method!="auto" && method!="spectral" && method!="taylor")throw std::invalid_argument("FT transport must be auto, spectral or taylor");
          numerical.ordinary_method=method=="auto"?diffexp::recursion::OrdinaryMethod::automatic:method=="spectral"?diffexp::recursion::OrdinaryMethod::spectral:diffexp::recursion::OrdinaryMethod::taylor;
        }
        else if(command=="ft" && option=="--method") {
          const std::string method=option_value;
          if(method!="adjoint" && method!="values" && method!="factored" && method!="auto")
            throw std::invalid_argument("FT method must be adjoint, factored, auto or values");
          numerical.observable_adjoint=method!="values";
          numerical.linear_method=method=="factored"?diffexp::recursion::LinearMethod::factored:
            method=="auto"?diffexp::recursion::LinearMethod::automatic:diffexp::recursion::LinearMethod::adjoint;
        }
        else if(command=="ft") {
          const std::string value=option_value;unsigned number=0;
          const auto parsed=std::from_chars(value.data(),value.data()+value.size(),number);
          if(parsed.ec!=std::errc{} || parsed.ptr!=value.data()+value.size())throw std::invalid_argument("FT numerical option must be a nonnegative integer");
          if(option=="--epsilon-order") {
            if(number>100)throw std::invalid_argument("epsilon order must be in 0..100");
            epsilon_order=number;epsilon_order_given=true;
          }
          else if(option=="--endpoint-order")numerical.endpoint_order=number;
          else if(option=="--ordinary-order")numerical.ordinary_order=number;
          else if(option=="--precision-bits")numerical.working_bits=number;
          else if(option=="--leaf-digits")numerical.leaf_digits=number;
          else if(option=="--ft-transport-digits"){if(!number || number>100000)throw std::invalid_argument("FT transport digits must be in 1..100000");numerical.spectral.accuracy_goal=number;}
          else throw std::invalid_argument("unknown FT option: "+option);
        } else throw std::invalid_argument("unknown prepare option: "+option);
      }
      options.reduction.batch_cache_directory=cache/"fire-batches";
      options.reduction.pending_cache_directory=cache/"fire-pending";
      if(numerical_cache){numerical.endpoint_cache_directory=cache/"affine-series";numerical.ordinary_cache_directory=cache/"ordinary-transport";}
      auto requested=configuration.integrals;
      std::optional<diffexp::henn::CanonicalBasis> canonical_basis;
      if(!henn_basis.empty()) {
        if(!diffexp::family_config::same_geometry(configuration.family,diffexp::feynman::example_family("henn_double_pentagon_x0")))
          throw std::invalid_argument("--henn-basis requires the published Henn X0 geometry and denominator order; a family label alone is insufficient");
        canonical_basis=diffexp::henn::read_x0(henn_basis);
        requested=canonical_basis->scalar_targets;
        if(!epsilon_order_given)epsilon_order=4;
      }
      diffexp::artifacts::Store store(cache);unsigned reused=0,built=0;
      auto provider=[&](const auto& basis,const auto& dimension,const auto& field,auto parameter,const auto& sources,const auto& budget) {
        diffexp::level::Provider fallback;
        if(!fire_prime.empty()) {
          diffexp::fire_modular::Options modular;modular.executable=fire_prime;modular.cache_directory=cache/"fire-modular";
          modular.progress=[](const std::string& message){std::cerr<<message<<'\n';};
          auto session=std::make_shared<diffexp::fire_modular::Session>(basis,dimension,field,modular,budget.batch_cache_directory);
          fallback=[session](const auto& batch,const auto& limits){return (*session)(batch,limits);};
        }
        auto result=diffexp::cached_level::prepare(store,basis,dimension,field,parameter,sources,budget,fallback);
        if(result.cache_hit)++reused;else if(result.result.success)++built;
        return std::move(result.result);
      };
      auto graph=diffexp::recursion::prepare(configuration.family,std::move(requested),options,provider,
        [&](const auto& event) {
          progress<<"Level "<<event.depth<<": "<<event.physical_propagators<<" propagators; ";
          if(event.scalar_leaf)progress<<"analytic scalar boundary";
          else progress<<event.masters<<" masters; "<<event.sources<<" shared source integrals";
          progress<<"; elapsed "<<event.elapsed_seconds<<" s\n"<<std::flush;
        });
      progress<<"Exact recursion prepared: "<<graph.nodes.size()<<" levels, "<<built<<" systems built, "<<reused<<" systems reused.\n"
        <<"Verified cache: "<<std::filesystem::absolute(cache)<<"\n";
      boost::json::object report{{"schema","DiffExp.FeynmanTrick/v1"},{"family",configuration.family.name},
        {"status",command=="prepare"?"prepared":"completed"},{"levels",graph.nodes.size()},
        {"systems_built",built},{"systems_reused",reused}};
      report["timings"]=boost::json::object{{"preparation_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-preparation_start).count()}};
      if(command=="prepare")progress<<"This prepares equations and recursive operations; it does not report a numerical integral.\n";
      else {
        numerical.progress=[&progress,start=std::chrono::steady_clock::now()](std::size_t depth,const std::string& phase,int high) {
          const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
          progress<<"Numerical level "<<depth+1<<": "<<phase<<"; epsilon through "<<high<<"; elapsed "<<elapsed<<" s\n"<<std::flush;
        };
        const auto numerical_start=std::chrono::steady_clock::now();
        const auto raw_high=canonical_basis?std::max(0,diffexp::henn::needed_scalar_high(*canonical_basis,epsilon_order)):static_cast<int>(epsilon_order);
        diffexp::recursion::Evaluator evaluator(graph,numerical);auto result=evaluator.evaluate(raw_high);
        report["timings"].as_object()["numerical_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-numerical_start).count();
        if(canonical_basis)result=diffexp::henn::project_boundary(*canonical_basis,graph.nodes.front().requested,result,epsilon_order,numerical.working_bits);
        if(json_mode) {
          // Display enough decimal digits for the working binary precision;
          // Arb adds outward rounding to the printed interval radius.
          const long output_digits=std::max<long>(45,(static_cast<long>(numerical.working_bits)*30103+99999)/100000+3);
          boost::json::array rows;
          for(const auto& row:result.values) {
            boost::json::array powers;
            for(int k=result.low;k<=static_cast<int>(epsilon_order);++k)
              powers.push_back(diffexp::ball_json(row.at(k-result.low),output_digits));
            rows.push_back(std::move(powers));
          }
          report["epsilon_low"]=result.low;report["epsilon_high"]=epsilon_order;
          report["coefficient_order"]="component,epsilon";
          report["coefficients"]=std::move(rows);
          report["component_basis"]=canonical_basis?"canonical":"requested_scalar";
          report["ball_encoding"]="Arb decimal interval strings for real and imaginary parts";
          report["omitted_tails_certified"]=result.taylor_tail_certified;
          report["working_bits"]=numerical.working_bits;
          report["persistent_numerical_cache"]=numerical_cache;
        } else {
          if(!result.taylor_tail_certified)progress<<"Retained series values (printed radii exclude omitted endpoint and transport tails):\n";
          for(unsigned row=0;row<result.values.size();++row)for(int k=result.low;k<=static_cast<int>(epsilon_order);++k) {
            progress<<(canonical_basis?"Canonical component ":"Scalar ")<<row+1<<" epsilon^"<<k<<" = ";acb_printn(result.values[row][k-result.low].raw(),32,0);progress<<'\n';
          }
        }
        const auto& stats=evaluator.statistics();
        boost::json::array rejected;for(const auto& reason:stats.spectral_rejections)rejected.emplace_back(reason);
        report["ft_transport"]=boost::json::object{{"spectral_attempts",stats.spectral_attempts},{"spectral_accepted",stats.spectral_accepted},{"spectral_legs",stats.spectral_legs},{"spectral_reused",stats.spectral_reused},{"rejections",rejected}};
        report["timings"].as_object()["ordinary_seconds"]=stats.ordinary_seconds;
        report["timings"].as_object()["spectral_seconds"]=stats.spectral_seconds;
        progress<<"Exact endpoint plans "<<stats.exact_plans<<"; series built "<<stats.endpoint_series_built<<"; verified series reused "<<stats.endpoint_series_reused<<"; numerical refinements "<<stats.refinements<<"; numerical cache hits "<<stats.cache_hits<<".\n"
          <<"Linear transport selections: adjoint "<<stats.adjoint_selections<<", factored "<<stats.factored_selections<<"; conditioning method fallbacks "<<stats.conditioning_method_fallbacks<<".\n"
          <<"Ordinary continuation: "<<stats.ordinary_checkpoints.loaded<<" saved arms resumed, "<<stats.ordinary_checkpoints.completed_reused
          <<" already complete; "<<stats.ordinary_checkpoints.saved<<" accepted charts saved.\n"
          <<"Child-window preflights "<<stats.demand_preflights<<"; local operator reuses "<<stats.local_operator_reuses<<".\n"
          <<"Centered chart actions "<<conditioning.centered_charts<<"; conditioning subdivisions "<<conditioning.conditioning_subdivisions<<".\n"
          <<(result.taylor_tail_certified?"Certified scalar-leaf coefficient enclosures.\n":"Numerical recursive series result; omitted endpoint and transport tails are not certified.\n");
        if(canonical_basis) {
          auto poles=diffexp::henn::audit_negative_poles(result,diffexp::Rational("1/100000000000000000000"),numerical.working_bits);
          progress<<"Canonical negative-pole audit at absolute tolerance 1e-20: "<<(poles.pass?"passed":"FAILED")
            <<" ("<<poles.checked_coefficients<<" coefficients checked).\n";
          report["negative_poles_pass"]=poles.pass;
          if(!poles.pass)throw std::runtime_error("canonical negative-pole audit failed");
        }
      }
      report["timings"].as_object()["total_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-preparation_start).count();
      if(json_mode)std::cout<<boost::json::serialize(report)<<'\n';
#else
      throw std::invalid_argument("native FIRE preparation currently requires a POSIX platform");
#endif
    } else if(command == "singular-endpoint") {
      if(argc!=2)throw std::invalid_argument("singular-endpoint takes no arguments; its JSON includes chart metadata and plot samples");
      std::cout<<boost::json::serialize(diffexp::singular_endpoint_example())<<'\n';
    } else if (command == "sunrise") {
      if(argc>3)throw std::invalid_argument("sunrise accepts an optional epsilon order");
      unsigned order=0;
      if(argc==3) {
        const std::string value=argv[2];auto parsed=std::from_chars(value.data(),value.data()+value.size(),order);
        if(parsed.ec!=std::errc{} || parsed.ptr!=value.data()+value.size())throw std::invalid_argument("epsilon order must be nonnegative integer");
      }
      diffexp::Jet::Ball::set_precision(384);auto result=diffexp::feynman::sunrise(order);
      for(unsigned k=0;k<=order;++k)std::cout<<"Sunrise epsilon^"<<k<<" = "<<result.at(k).real_midpoint(30)<<'\n';
      std::cout<<"Native two-level Feynman recursion; numerical series result with omitted Taylor tails not yet certified.\n";
    } else if(command == "banana4-oracle") {
      if(argc>3 || (argc==3 && std::string(argv[2])!="equal" && std::string(argv[2])!="unequal"))
        throw std::invalid_argument("banana4-oracle accepts equal or unequal");
      auto result=argc==3 && std::string(argv[2])=="unequal"?diffexp::oracle::banana4_unequal_bessel():diffexp::oracle::banana4_equal_bessel();
      std::cout<<"Independent certified Banana4 finite-part reference: ";acb_printn(result.value.raw(),40,0);
      std::cout<<"\nBessel integral with rigorous endpoint tails; this is not the Feynman-trick solver.\n";
    } else if (command == "bubble") {
      if(argc!=2)throw std::invalid_argument("bubble takes no arguments");
      auto result=diffexp::certified_feynman_bubble();
      for(unsigned k=0;k<result.coefficients.size();++k) {
        std::cout<<"Bubble epsilon^"<<k<<" = ";acb_printn(result.coefficients[k].raw(),30,0);std::cout<<'\n';
      }
      std::cout<<"Certified absolute radius <= "<<result.absolute_tolerance<<"; "<<result.charts<<" series charts.\n";
    } else if (command == "banana-unequal") {
      if((argc!=3 && argc!=5) || (argc==5 && std::string(argv[3])!="--route"))
        throw std::invalid_argument("banana-unequal requires ORIGINAL_DIR [--route momentum-first|mass-first]");
      const auto route=diffexp::original_banana_unequal_route(argc==5?argv[4]:"momentum-first");
      return diffexp::run_original_banana_unequal_route(argv[2],route);
    } else if (command == "banana-equal") {
      if(argc<3 || argc>7 || argc%2==0)throw std::invalid_argument("banana-equal requires ORIGINAL_DIR [--route real|contour] [--to ENDPOINT]");
      std::string route="real";diffexp::Rational endpoint(20);
      for(int i=3;i<argc;i+=2) {
        const std::string option=argv[i];
        if(option=="--route")route=argv[i+1];
        else if(option=="--to")endpoint=diffexp::Rational(argv[i+1]);
        else throw std::invalid_argument("unknown equal banana option");
      }
      if(route!="real" && route!="contour")
        throw std::invalid_argument("equal banana route must be real or contour");
      if(route=="contour" && endpoint!=diffexp::Rational(20))
        throw std::invalid_argument("the saved contour comparison ends at20; later endpoints use the real route");
      return route=="real"?diffexp::run_original_banana_equal_real(argv[2],endpoint):diffexp::run_original_banana_equal(argv[2]);
    } else if (command == "mpl") {
      if(argc>3 || (argc==3 && std::string(argv[2])!="--short"))throw std::invalid_argument("mpl accepts only --short");
      return diffexp::run_mpl(argc==2);
    } else if (command == "planar") {
      if(argc!=4)throw std::invalid_argument("planar requires ancillary directory and family");
      return diffexp::run_planar(argv[2],argv[3]);
    } else if (command == "series") {
      if(argc!=3)throw std::invalid_argument("series requires one JSON request file");
      return diffexp::run_series_file(argv[2]);
    } else if (command == "henn-nonplanar") {
      if(argc<3 || argc>4)throw std::invalid_argument("henn-nonplanar requires ancillary directory and optional Taylor order");
      unsigned order=140;
      if(argc==4) {
        const std::string value=argv[3];
        auto parsed=std::from_chars(value.data(),value.data()+value.size(),order);
        if(parsed.ec!=std::errc{} || parsed.ptr!=value.data()+value.size())
          throw std::invalid_argument("Henn Taylor order must be a positive integer");
      }
      return diffexp::run_henn_nonplanar(argv[2],order);
    } else if (command == "backend-info") {
#ifdef DIFFEXP_WITH_KERNEL_RUNTIME
      auto information=boost::json::parse(diffexp::kernel::backend_info_json()).as_object();
      information["compatibility_runtime"]=true;
#else
      boost::json::object information{{"standalone_cpp",true},{"librarylink",false},{"compatibility_runtime",false}};
#endif
      information["package_version"]="2.1.1";information["native_transport"]=true;
      information["feynman_configuration_schema"]="DiffExp.FeynmanFamily/v1";
      information["transport_schema"]="DiffExp.Transport/v1";
      std::cout<<boost::json::serialize(information)<<'\n';
    } else if (command == "kernel") {
#ifdef DIFFEXP_WITH_KERNEL_RUNTIME
      // One request/response per line; persistent handles live until EOF.
      std::string request;
      while (std::getline(std::cin, request)) {
        if (!request.empty())
          std::cout << diffexp::kernel::run_recurrence_json(request) << std::endl;
      }
      diffexp::kernel::reset_solver_sessions();
#else
      throw std::invalid_argument("prepared-input compatibility runtime is not enabled in this build");
#endif
    } else if (command == "--help") {
      std::cout << "Usage: diffexp --version | backend-info | series FILE_OR_MINUS | transport FILE_OR_MINUS | singular-endpoint | kernel | henn-nonplanar DATA [ORDER] | planar DATA FAMILY | mpl [--short] | banana-equal ORIGINAL_DIR | banana-unequal ORIGINAL_DIR | bubble | sunrise [EPSILON_ORDER] | banana4-oracle [equal|unequal] | family-template NAME | prepare CONFIGURATION [--fire PATH] [--cache DIRECTORY] [--henn-basis FILE] | ft CONFIGURATION [--fire PATH] [--cache DIRECTORY] [--henn-basis FILE] [--epsilon-order N] [--endpoint-order N] [--ordinary-order N] [--precision-bits N] [--leaf-digits N]\n";
      std::cout << "Use --json with prepare/ft for one JSON response; progress is written to stderr.\nPreparation budgets for prepare/ft: --fire-seconds N --level-seconds N --total-seconds N.\n";
      std::cout << "Modular preparation: --fire-prime /path/to/FIRE7p; durable samples under --cache, no symbolic fallback.\n";
      std::cout << "Provider resources: --fire-threads N --fire-simplifier-threads N --fire-memory-mib N. FT ordinary solver: --ft-transport auto|spectral|taylor (default auto); --ft-transport-digits N overrides the local spectral accuracy target (default leaf digits + 12). FT execution: --method adjoint|factored|auto|values (default adjoint); --no-numerical-cache disables numerical checkpoint files while retaining exact reductions.\n";
      std::cout << "Direct singular example: singular-endpoint (JSON chart, sectors, endpoint and plot samples).\n";
      std::cout << "Original equal banana: --route real|contour (default real; local +i0 singular matching); --to ENDPOINT (default 20; real route supports endpoints >=20, including historical 32).\n";
      std::cout << "Original unequal banana: --route momentum-first|mass-first (default momentum-first).\n";
    } else {
      std::cerr << "Unknown command: " << command << '\n';
      return 2;
    }
  } catch (const std::exception& e) {
    if(json_mode)std::cout<<boost::json::serialize(boost::json::object{{"schema","DiffExp.Error/v1"},{"status","error"},{"message",e.what()}})<<'\n';
    std::cerr << "diffexp: " << e.what() << '\n';
    return 1;
  }
}

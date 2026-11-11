#include "diffexp/fire_modular.hpp"
#include "diffexp/families.hpp"
#include <iostream>
using namespace diffexp;
int main(int argc,char** argv){try{
  ExactField field({"d"});Exact d(field,"d");
  auto text=fire_modular::detail::table({{{1,0},{{{1,0},d.constant(1)}}},{{2,0},{{{1,0},(d-d.constant(2))/d.constant(3)}}}},fire::SymbolMap(field.variables()));
  fire::SymbolMap symbols(field.variables());auto parsed=fire::import_table_text(text,field,2,&symbols);
  if(parsed.at({2,0}).at({1,0})!=(d-d.constant(2))/d.constant(3))throw std::runtime_error("modular reconstructed table roundtrip");
  // The active variables are those of the dimension and the contraction coefficients.
  {ExactField wide({"x","unused","dim"});Exact x(wide,"x"),dim(wide,"dim");
    auto raw=ibp::quadratic_family(feynman::banana(1,{Rational(1),Rational(2)}),dim);raw.physical[0].constant=x/(x+dim.constant(1));
    if(fire_modular::detail::active_variables(ibp::PropagatorBasis(raw),dim)!=std::vector<std::size_t>{0,2})throw std::runtime_error("active variable detection");
    // A text table of residues is what a table of the same constants was.
    fire::Residues words{{{1,0},{{{1,0},1}}},{{2,0},{{{1,0},12345}}}};
    if(fire::import_table_text(fire_modular::detail::table(words),field,2)!=fire::import_table_text(fire_modular::detail::table({{{1,0},{{{1,0},d.constant(1)}}},{{2,0},{{{1,0},d.constant(12345)}}}},fire::SymbolMap(field.variables())),field,2))
      throw std::runtime_error("residue table codec");
    fire_modular::detail::validate(words,ibp::PropagatorBasis(ibp::quadratic_family(feynman::banana(1,{Rational(1),Rational(2)}),d)),{{2,0}});
    bool rejected=false;words[{2,0}][{3,0}]=5;try{fire_modular::detail::validate(words,ibp::PropagatorBasis(ibp::quadratic_family(feynman::banana(1,{Rational(1),Rational(2)}),d)),{{2,0}});}catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("residues with a missing terminal row accepted");}
  // Residues are read once into machine words; absent terms are zero.
  {const auto p=modular::prime(1);const auto third=d.constant(7)/d.constant(3);
    auto words=fire_modular::detail::residues({{{1,0},{{{1,0},d.constant(1)}}},{{2,0},{{{1,0},third}}}},{5},p);fire_modular::detail::Probe sample{{5},words,7};
    if(fire_modular::detail::coefficient(sample,{{2,0},{1,0}})!=modular::rational_mod(Rational("7/3"),p)||fire_modular::detail::coefficient(sample,{{1,0},{1,0}})!=1||fire_modular::detail::coefficient(sample,{{1,0},{2,0}})!=0)
      throw std::runtime_error("finite-field residue rows");
    if(fire_modular::detail::masters(sample)!=std::set<ibp::Integral>{{1,0}})throw std::runtime_error("residue row terminal integrals");
    bool nonconstant=false;try{fire_modular::detail::residues({{{2,0},{{{1,0},d}}}},{5},p);}catch(const std::runtime_error&){nonconstant=true;}
    if(!nonconstant)throw std::runtime_error("symbolic finite-field output accepted");}
  // Ordinary, slice, planning and held-out points never share an ordinal.
  namespace ordinals=fire_modular::ordinals;
  if(ordinals::slice(4096)||!ordinals::slice(ordinals::slice_base)||!ordinals::slice(ordinals::slice_base+ordinals::slice_stride*ordinals::slice_axes-1)||ordinals::slice(ordinals::planning_base)||
     ordinals::validation(ordinals::planning_base+32)||!ordinals::validation(ordinals::validation_point(2,0))||ordinals::validation_point(14,2)!=101402)throw std::runtime_error("sample ordinal ranges overlap");
  if(argc==1){std::cout<<"modular provider codecs passed; supply FIRE7p for integration\n";return 0;}
  fire_modular::Options opts;opts.executable=argv[1];
  if(argc>2)opts.cache_directory=argv[2];else {auto name=(std::filesystem::temp_directory_path()/"de3-modular-validation-XXXXXX").string();std::vector<char> tmp(name.begin(),name.end());tmp.push_back(0);if(!mkdtemp(tmp.data()))throw std::runtime_error("cannot create private modular validation directory");opts.cache_directory=tmp.data();}
  opts.progress=[](const auto& s){std::cerr<<s<<'\n';};
  ibp::PropagatorBasis basis(ibp::quadratic_family(feynman::banana(1,{Rational(1),Rational(2)}),d));
  fire_modular::Session session(basis,d,field,opts);fire::Options limits;limits.timeout_seconds=120;
  auto reduced=session({{2,1},{1,2}},limits);if(!reduced.success)throw std::runtime_error(reduced.reason);
  auto semantic=fire_batch::identity(basis,d,field,{{2,1},{1,2}},limits);
  auto legacy_identity=artifacts::detail::sha256(artifacts::detail::canonical(boost::json::object{{"schema","DiffExp3.NativeModularProvider/v1"},{"science",semantic.json_value()},{"active",boost::json::array{0}},{"prime_table","FIRE7-primes-1-through-16"}}));
  if(reduced.directory.filename()!=legacy_identity)throw std::runtime_error("existing FIRE checkpoint identity changed");
  ibp::Generator generator(basis,d);ibp::ExactReducer reducer(d);
  for(int a=-1;a<=4;++a)for(int b=-1;b<=4;++b)for(auto row:generator.relations({a,b}))reducer.insert(std::move(row));
  for(const auto& [a,row]:reduced.reductions){ibp::Relation residual{{a,d.constant(1)}};ibp::add_scaled(residual,row,d.constant(-1));
    if(!reducer.reduce(residual).remainder.empty())throw std::runtime_error("modular FIRE differs from independent exact native IBPs");}
  fire_modular::Session resumed(basis,d,field,opts);auto again=resumed({{2,1},{1,2}},limits);
  if(!again.success||again.reductions!=reduced.reductions)throw std::runtime_error("modular completion recovery differs");
  // Losing the final record must not repeat completed finite-field work.
  const auto completed=reduced.directory/"completed.json";auto completion=fire_modular::detail::read(completed);
  std::filesystem::remove(completed);auto offline_options=opts;offline_options.executable="/nonexistent/FIRE7p";
  fire_modular::Session offline(basis,d,field,offline_options);auto rebuilt=offline({{2,1},{1,2}},limits);
  if(!rebuilt.success||rebuilt.reductions!=reduced.reductions)throw std::runtime_error("durable sample recovery required fresh FIRE work: "+rebuilt.reason);
  // A valid envelope checksum is not sufficient: independently stored probes
  // must still reject a changed analytic result on every cache load.
  auto changed=reduced.reductions;changed.at({2,1}).begin()->second=changed.at({2,1}).begin()->second+d.constant(1);
  auto altered=completion;altered["tables"]=fire_modular::detail::table(changed,symbols);fire_modular::detail::save(completed,altered);
  fire_modular::Session invalid(basis,d,field,offline_options);auto rejected=invalid({{2,1},{1,2}},limits);
  if(rejected.success||rejected.reason.find("independent sample replay")==std::string::npos)throw std::runtime_error("changed analytic reconstruction accepted");
  fire_modular::detail::save(completed,completion);
  const auto sample_file=reduced.directory/artifacts::detail::string(completion.at("validation_samples").as_array().front());
  auto sample_text=fire::read_text(sample_file);fire_modular::detail::publish(sample_file,sample_text+"x");
  fire_modular::Session damaged(basis,d,field,offline_options);auto damaged_result=damaged({{2,1},{1,2}},limits);
  fire_modular::detail::publish(sample_file,sample_text);
  if(damaged_result.success)throw std::runtime_error("damaged independent probe accepted");
  std::cout<<"modular FIRE bubble reductions, independent exact native witnesses, offline sample recovery and corrupted result/probe rejection passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

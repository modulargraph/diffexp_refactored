#pragma once
#include "diffexp/numerical_rows_io.hpp"
#include <charconv>

namespace diffexp::adjoint_checkpoint {
namespace json=boost::json;namespace fs=std::filesystem;
struct Statistics {std::size_t loaded=0,saved=0,completed_reused=0;};
struct Options {
  fs::path directory;
  std::size_t max_bytes=64*1024*1024;
  Statistics* statistics=nullptr;
};
// The algorithm revision must change when retained arithmetic or chart
// selection changes. Observers and counters are deliberately not scientific
// inputs; precision, budgets, exact equations and original balls are.
inline std::string identity(const ExactEpsilonMatrix& matrix,const LaurentRows& initial,
    const ExactEpsilonMatrix& forcing,const std::vector<Exact>& vertices,const AdjointOptions& options) {
  if(vertices.empty())throw std::invalid_argument("checkpoint path is empty");
  json::array symbols,path;
  for(const auto& name:vertices.front().variables())symbols.emplace_back(name);
  const auto check_field=[&](const Exact& entry) {
    if(entry.variables()!=vertices.front().variables())throw std::invalid_argument("checkpoint exact field mismatch");
  };
  for(const auto& vertex:vertices){check_field(vertex);path.emplace_back(vertex.str());}
  for(const auto* m:{&matrix,&forcing})for(const auto& row:*m)for(const auto& entry:row)check_field(entry);
  for(const auto& row:initial.coefficients)for(const auto& series:row)for(const auto& b:series)
    if(!b.is_finite())throw std::invalid_argument("checkpoint initial ball is nonfinite");
  const json::object payload{{"algorithm","DiffExp3.AdjointOriginalPath/v1"},
    {"flint_version",FLINT_VERSION},{"precision_bits",Jet::Ball::precision()},
    {"field_symbols",symbols},{"connection",numerical_rows_io::exact_matrix(matrix)},
    {"forcing",numerical_rows_io::exact_matrix(forcing)},{"initial",numerical_rows_io::exact_rows(initial)},
    {"path",path},{"taylor_order",options.taylor_order},{"max_charts_per_leg",options.max_charts_per_leg},
    {"max_taylor_cells",options.max_taylor_cells},{"polynomial_recurrence",options.polynomial_recurrence},
    {"centered_input",options.centered_input},{"max_centered_map_cells",options.max_centered_map_cells},
    {"max_conditioning_halvings",options.max_conditioning_halvings},{"max_rows_per_batch",options.max_rows_per_batch}};
  return artifacts::detail::sha256(artifacts::detail::canonical(payload));
}
namespace detail {
inline void limits(const Options& options) {
  if(!options.max_bytes || options.max_bytes>1024ULL*1024*1024)
    throw std::invalid_argument("ordinary checkpoint file budget");
}
inline void publish(const fs::path& target,const std::string& bytes) {
  fs::create_directories(target.parent_path());
  auto name=target.string()+".new.XXXXXX";std::vector<char> tmp(name.begin(),name.end());tmp.push_back(0);
  int fd=::mkstemp(tmp.data());if(fd<0)throw std::runtime_error("cannot create ordinary checkpoint staging file");
  try {
    std::size_t offset=0;
    while(offset<bytes.size()) {
      const auto n=::write(fd,bytes.data()+offset,bytes.size()-offset);
      if(n<0 && errno==EINTR)continue;
      if(n<=0)throw std::runtime_error("cannot write ordinary checkpoint");
      offset+=static_cast<std::size_t>(n);
    }
    if(::fsync(fd))throw std::runtime_error("cannot sync ordinary checkpoint");
    const auto status=::close(fd);fd=-1;if(status)throw std::runtime_error("cannot close ordinary checkpoint");
    fs::rename(tmp.data(),target);artifacts::detail::sync_directory(target.parent_path());
  } catch(...) {if(fd>=0)::close(fd);::unlink(tmp.data());throw;}
}
inline json::value read(const fs::path& path,std::size_t max_bytes) {
  const auto size=fs::file_size(path);
  if(size>max_bytes)throw std::length_error("ordinary checkpoint exceeds file budget");
  std::ifstream in(path,std::ios::binary);std::string bytes(static_cast<std::size_t>(size),'\0');
  if(!in.read(bytes.data(),static_cast<std::streamsize>(size)) || in.peek()!=std::char_traits<char>::eof())
    throw std::runtime_error("ordinary checkpoint read failed or size changed");
  return json::parse(bytes);
}
inline json::object encode(const AdjointContinuation& saved,const std::string& key) {
  return {{"schema","DiffExp3.AdjointContinuation/v1"},{"identity",key},
    {"certificate","uncertified retained numerical checkpoint"},{"leg",saved.leg},
    {"accepted_charts",saved.accepted_charts},
    {"parameter_ieee_bits",std::to_string(std::bit_cast<std::uint64_t>(saved.parameter))},
    {"rows",numerical_rows_io::exact_rows(saved.rows)}};
}
inline AdjointContinuation decode(const json::value& envelope,const std::string& key) {
  const auto& root=envelope.as_object();artifacts::detail::keys(root,{"payload","sha256"});
  const auto& payload=root.at("payload");
  if(artifacts::detail::string(root.at("sha256"))!=artifacts::detail::sha256(artifacts::detail::canonical(payload)))
    throw std::invalid_argument("ordinary checkpoint checksum mismatch");
  const auto& o=payload.as_object();
  artifacts::detail::keys(o,{"schema","identity","certificate","leg","accepted_charts","parameter_ieee_bits","rows"});
  if(o.at("schema")!="DiffExp3.AdjointContinuation/v1" || artifacts::detail::string(o.at("identity"))!=key ||
      o.at("certificate")!="uncertified retained numerical checkpoint")
    throw std::invalid_argument("ordinary checkpoint identity/schema mismatch");
  const auto leg=artifacts::detail::integer(o.at("leg")),charts=artifacts::detail::integer(o.at("accepted_charts"));
  const auto bits=artifacts::detail::string(o.at("parameter_ieee_bits"));std::uint64_t parsed=0;
  const auto converted=std::from_chars(bits.data(),bits.data()+bits.size(),parsed);
  if(leg<0 || charts<1 || converted.ec!=std::errc{} || converted.ptr!=bits.data()+bits.size() || std::to_string(parsed)!=bits)
    throw std::invalid_argument("ordinary checkpoint continuation metadata");
  return {static_cast<unsigned>(leg),static_cast<unsigned>(charts),std::bit_cast<double>(parsed),numerical_rows_io::read_rows(o.at("rows"))};
}
} // namespace detail

// Optional I/O adapter around the pure recurrence. A different problem gets
// a different key; a corrupt matching file is an error, never a silent restart.
// Full input rows bind factored maps to their child expression as well.
inline LaurentRows transport(const ExactEpsilonMatrix& matrix,LaurentRows initial,
    const ExactEpsilonMatrix& forcing,const std::vector<Exact>& vertices,
    const AdjointOptions& numerical,const Options& storage={}) {
  if(storage.directory.empty())return transport_adjoint_rows(matrix,std::move(initial),forcing,vertices,numerical);
  detail::limits(storage);
  if(numerical.continuation)throw std::invalid_argument("explicit continuation cannot be combined with ordinary checkpoint storage");
  const auto key=identity(matrix,initial,forcing,vertices,numerical);
  const auto file=storage.directory/(key+".json");auto options=numerical;
  if(fs::exists(file)) {
    auto saved=detail::decode(detail::read(file,storage.max_bytes),key);
    options.continuation=std::make_shared<const AdjointContinuation>(std::move(saved));
  }
  const auto observer=options.continuation_observer;
  options.continuation_observer=[&](const AdjointContinuation& saved) {
    const auto payload=detail::encode(saved,key);
    const auto bytes=artifacts::detail::canonical(json::object{{"payload",payload},
      {"sha256",artifacts::detail::sha256(artifacts::detail::canonical(payload))}});
    if(bytes.size()>storage.max_bytes)throw std::length_error("ordinary checkpoint exceeds file budget");
    detail::publish(file,bytes);
    if(storage.statistics)++storage.statistics->saved;
    if(observer)observer(saved);
  };
  auto result=transport_adjoint_rows(matrix,std::move(initial),forcing,vertices,options);
  if(storage.statistics && options.continuation) {
    ++storage.statistics->loaded;
    const auto& saved=*options.continuation;
    bool complete=saved.parameter==1;
    for(std::size_t leg=saved.leg+1;leg+1<vertices.size();++leg)complete&=vertices[leg]==vertices[leg+1];
    if(complete)++storage.statistics->completed_reused;
  }
  return result;
}
} // namespace diffexp::adjoint_checkpoint

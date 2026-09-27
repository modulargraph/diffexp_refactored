#pragma once
#include "diffexp/ft_spectral_checkpoint.hpp"
#include "diffexp/ft_ultraspherical.hpp"

namespace diffexp::ft_ultraspherical_checkpoint {
using Storage = ft_spectral_checkpoint::Storage;
namespace json = boost::json;
inline constexpr const char* algorithm = "DiffExp.FTUltrasphericalCompletedArm/v2";
inline std::string identity(const ExactEpsilonMatrix& matrix,
    const LaurentRows& initial, const ExactEpsilonMatrix& forcing,
    const std::vector<Exact>& vertices, const ft_ultraspherical::Options& options,
    const AdjointOptions& native = {}) {
  return artifacts::detail::sha256(artifacts::detail::canonical(json::object{
      {"algorithm",algorithm},
      {"spectral_identity",ft_spectral_checkpoint::identity(matrix,initial,forcing,vertices,options,native)},
      {"working_bits",options.working_bits}}));
}
inline std::optional<LaurentRows> try_transport(
    const ExactEpsilonMatrix& matrix, const LaurentRows& initial,
    const ExactEpsilonMatrix& forcing, const std::vector<Exact>& vertices,
    const ft_ultraspherical::Options& options, ft_spectral::Diagnostics& diagnostics,
    const AdjointOptions& native = {}, const Storage& storage = {}) {
  if(storage.directory.empty())
    return ft_ultraspherical::try_transport(matrix,initial,forcing,vertices,options,diagnostics);
  adjoint_checkpoint::detail::limits({storage.directory,storage.max_bytes,nullptr});
  if(native.continuation)throw std::invalid_argument("FT ultraspherical completed cache cannot use partial continuation");
  if(vertices.empty()||matrix.size()!=initial.columns()||forcing.size()!=initial.coefficients.size())
    throw std::invalid_argument("FT ultraspherical checkpoint input shape");
  for(const auto* m:{&matrix,&forcing})for(const auto&row:*m)
    if(row.size()!=initial.columns())throw std::invalid_argument("FT ultraspherical checkpoint matrix shape");
  const auto [xi,ei]=path_epsilon_variables(vertices.front());(void)xi;
  const auto low=ft_spectral_checkpoint::detail::expected_low(initial,forcing,ei);
  const auto legs=ft_spectral_checkpoint::detail::leg_count(vertices);
  const auto key=identity(matrix,initial,forcing,vertices,options,native);
  const auto file=storage.directory/(key+".json");
  if(std::filesystem::exists(file)) {
    auto envelope=adjoint_checkpoint::detail::read(file,storage.max_bytes);
    const auto& root=envelope.as_object();
    artifacts::detail::keys(root,{"payload","sha256"});
    auto payload=root.at("payload").as_object();
    if(payload.at("schema")!=algorithm || artifacts::detail::string(root.at("sha256"))!=artifacts::detail::sha256(artifacts::detail::canonical(payload)))
      throw std::invalid_argument("FT ultraspherical checkpoint schema/checksum mismatch");
    // Reuse the typed decoder after validating our distinct outer schema/hash.
    payload["schema"]=ft_spectral_checkpoint::algorithm;
    auto rows=ft_spectral_checkpoint::detail::decode(json::object{
      {"payload",payload},{"sha256",artifacts::detail::sha256(artifacts::detail::canonical(payload))}},
      key,initial,low,legs,options,diagnostics);
    if(storage.reuse_counter)++*storage.reuse_counter;
    if(options.progress)for(unsigned leg=0;leg+1<vertices.size();++leg)
      if(vertices[leg]!=vertices[leg+1])options.progress(leg,1.);
    return rows;
  }
  auto rows=ft_ultraspherical::try_transport(matrix,initial,forcing,vertices,options,diagnostics);
  if(!rows)return std::nullopt;
  if(diagnostics.legs!=legs)throw std::logic_error("FT ultraspherical checkpoint requires a complete arm");
  ft_spectral_checkpoint::detail::validate(*rows,initial,low);
  // Store only the last retained grid diagnostic: the shared decoder bounds
  // collocation histories to nine grids, while this backend also supports N80.
  // This changes no scientific result or acceptance claim.
  auto retained=diagnostics;
  if(!retained.nodes.empty())retained.nodes={retained.nodes.back()};
  auto payload=ft_spectral_checkpoint::detail::encode(*rows,key,retained);
  payload["schema"]=algorithm;
  const auto bytes=artifacts::detail::canonical(json::object{{"payload",payload},
    {"sha256",artifacts::detail::sha256(artifacts::detail::canonical(payload))}});
  if(bytes.size()>storage.max_bytes)throw std::length_error("FT ultraspherical checkpoint file budget");
  adjoint_checkpoint::detail::publish(file,bytes);
  return rows;
}
}

#pragma once
#include "diffexp/scheduler.hpp"
#include <boost/json.hpp>
#include <array>
#include <bit>
#include <cerrno>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string_view>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace diffexp::artifacts {
namespace json=boost::json;
namespace fs=std::filesystem;
namespace detail {
// SHA-256 (FIPS 180-4), independent of the optional compatibility runtime.
inline std::string sha256(std::string_view bytes) {
  constexpr std::array<std::uint32_t,64> k={
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
  std::array<std::uint32_t,8> h={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
  std::vector<unsigned char> data(bytes.begin(),bytes.end());
  const auto length=static_cast<std::uint64_t>(data.size())*8;data.push_back(0x80);
  while(data.size()%64!=56)data.push_back(0);
  for(int shift=56;shift>=0;shift-=8)data.push_back(static_cast<unsigned char>(length>>shift));
  for(std::size_t offset=0;offset<data.size();offset+=64) {
    std::array<std::uint32_t,64> w{};
    for(unsigned i=0;i<16;++i)for(unsigned j=0;j<4;++j)w[i]=(w[i]<<8)|data[offset+4*i+j];
    for(unsigned i=16;i<64;++i) {
      auto s0=std::rotr(w[i-15],7)^std::rotr(w[i-15],18)^(w[i-15]>>3);
      auto s1=std::rotr(w[i-2],17)^std::rotr(w[i-2],19)^(w[i-2]>>10);
      w[i]=w[i-16]+s0+w[i-7]+s1;
    }
    auto a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],z=h[7];
    for(unsigned i=0;i<64;++i) {
      auto t1=z+(std::rotr(e,6)^std::rotr(e,11)^std::rotr(e,25))+((e&f)^(~e&g))+k[i]+w[i];
      auto t2=(std::rotr(a,2)^std::rotr(a,13)^std::rotr(a,22))+((a&b)^(a&c)^(b&c));
      z=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
    }
    h[0]+=a;h[1]+=b;h[2]+=c;h[3]+=d;h[4]+=e;h[5]+=f;h[6]+=g;h[7]+=z;
  }
  constexpr char hex[]="0123456789abcdef";std::string result;
  for(auto word:h)for(int shift=28;shift>=0;shift-=4)result.push_back(hex[(word>>shift)&15]);
  return result;
}
inline bool hash_string(const std::string& s) {
  return s.size()==64 && s.find_first_not_of("0123456789abcdef")==std::string::npos;
}
inline json::value ordered(const json::value& value,unsigned depth=0) {
  if(depth>128)throw std::invalid_argument("artifact JSON nesting exceeds limit");
  // Exact scientific values/balls must use exact strings. Floating-point JSON
  // numbers introduce ambiguous precision and are excluded from this codec.
  if(value.is_double())throw std::invalid_argument("artifact JSON requires exact strings instead of floating-point numbers");
  if(value.is_object()) {
    std::map<std::string,json::value> sorted;
    for(const auto& item:value.as_object())sorted.emplace(std::string(item.key()),ordered(item.value(),depth+1));
    json::object result;for(auto& [key,v]:sorted)result.emplace(key,std::move(v));return result;
  }
  if(value.is_array()) {json::array result;for(const auto& v:value.as_array())result.push_back(ordered(v,depth+1));return result;}
  return value;
}
inline std::string canonical(const json::value& value) {return json::serialize(ordered(value));}
inline void keys(const json::object& object,std::initializer_list<const char*> required) {
  if(object.size()!=required.size())throw std::invalid_argument("artifact schema has missing or unknown fields");
  for(auto key:required)if(!object.contains(key))throw std::invalid_argument("artifact schema is missing a required field");
}
inline std::string string(const json::value& value) {return std::string(value.as_string());}
inline long integer(const json::value& value) {
  if(!value.is_int64() || value.as_int64() < -1000000 || value.as_int64()>1000000)
    throw std::invalid_argument("artifact integer field outside supported range");
  return static_cast<long>(value.as_int64());
}
inline json::object demand_json(const Demand& d) {
  d.validate();
  if(d.epsilon_low < -1000000 || d.epsilon_high>1000000 || d.taylor_order>1000000 ||
      d.precision_bits>1000000 || d.publication_digits>1000000)
    throw std::invalid_argument("artifact resource guarantee exceeds limit");
  return {{"epsilon_low",d.epsilon_low},{"epsilon_high",d.epsilon_high},
      {"taylor_order",static_cast<std::int64_t>(d.taylor_order)},
      {"precision_bits",static_cast<std::int64_t>(d.precision_bits)},
      {"publication_digits",static_cast<std::int64_t>(d.publication_digits)}};
}
inline Demand demand_from(const json::value& value) {
  const auto& o=value.as_object();keys(o,{"epsilon_low","epsilon_high","taylor_order","precision_bits","publication_digits"});
  auto low=integer(o.at("epsilon_low")),high=integer(o.at("epsilon_high")),t=integer(o.at("taylor_order")),
       b=integer(o.at("precision_bits")),p=integer(o.at("publication_digits"));
  if(t<0 || b<0 || p<0)throw std::invalid_argument("negative artifact resource guarantee");
  Demand result{static_cast<int>(low),static_cast<int>(high),static_cast<unsigned>(t),static_cast<unsigned>(b),static_cast<unsigned>(p)};
  result.validate();return result;
}
inline void sync_directory(const fs::path& path) {
  int fd=::open(path.c_str(),O_RDONLY|O_DIRECTORY);
  if(fd<0)throw std::runtime_error("cannot open artifact directory for sync");
  const int status=::fsync(fd);::close(fd);
  if(status)throw std::runtime_error("cannot sync artifact directory");
}
} // namespace detail

struct ParentIdentity {std::string role,semantic_id;};
struct Identity {
  std::string kind; // exact_equation, receiving_basis, boundary, transport, observable
  std::string algorithm_version;
  json::object family,normalization,branch,geometry,boundary,scientific_inputs;
  json::array ordered_basis;
  std::vector<ParentIdentity> parents;
  json::object json_value() const {
    const std::set<std::string> kinds={"exact_equation","receiving_basis","boundary","transport","observable"};
    if(!kinds.contains(kind) || algorithm_version.empty() || family.empty() || ordered_basis.empty() ||
        normalization.empty() || branch.empty() || geometry.empty() || boundary.empty())
      throw std::invalid_argument("incomplete scientific artifact identity");
    // Explicit {status:not-applicable} is permitted for genuinely irrelevant
    // axes; absence is never silently interpreted as mathematical equality.
    if(kind!="exact_equation" && parents.empty())
      throw std::invalid_argument("derived artifact identity requires ordered semantic parents");
    json::array p;std::set<std::string> roles;
    for(const auto& parent:parents) {
      if(parent.role.empty() || !roles.insert(parent.role).second || !detail::hash_string(parent.semantic_id))
        throw std::invalid_argument("invalid artifact parent identity");
      p.push_back(json::object{{"role",parent.role},{"semantic_id",parent.semantic_id}});
    }
    return {{"schema","DiffExp3.ArtifactIdentity/v1"},{"kind",kind},{"algorithm_version",algorithm_version},
      {"family",family},{"ordered_basis",ordered_basis},{"normalization",normalization},{"branch",branch},
      {"geometry",geometry},{"boundary",boundary},{"parents",p},{"scientific_inputs",scientific_inputs}};
  }
  static Identity from_json(const json::object& o) {
    detail::keys(o,{"schema","kind","algorithm_version","family","ordered_basis","normalization","branch","geometry","boundary","parents","scientific_inputs"});
    if(detail::string(o.at("schema"))!="DiffExp3.ArtifactIdentity/v1")throw std::invalid_argument("unknown scientific identity schema");
    Identity result;result.kind=detail::string(o.at("kind"));result.algorithm_version=detail::string(o.at("algorithm_version"));
    result.family=o.at("family").as_object();result.ordered_basis=o.at("ordered_basis").as_array();
    result.normalization=o.at("normalization").as_object();result.branch=o.at("branch").as_object();
    result.geometry=o.at("geometry").as_object();result.boundary=o.at("boundary").as_object();
    result.scientific_inputs=o.at("scientific_inputs").as_object();
    for(const auto& parent:o.at("parents").as_array()) {
      const auto& p=parent.as_object();detail::keys(p,{"role","semantic_id"});
      result.parents.push_back({detail::string(p.at("role")),detail::string(p.at("semantic_id"))});
    }
    result.json_value();return result;
  }
  std::string key() const {return detail::sha256(detail::canonical(json_value()));}
};

// These are producer-supplied claims, not proofs invented by the store. Demand
// precision_bits/publication_digits describe resources and NEVER imply a
// certified numerical error. An independent verifier owns the evidence.
struct Certificate {
  std::string type="uncertified_numerical"; // or exact, certified_enclosure
  std::string verifier,scope;
  json::object evidence;
  void validate() const {
    if(type!="uncertified_numerical" && type!="exact" && type!="certified_enclosure")
      throw std::invalid_argument("unknown artifact certificate type");
    if(type!="uncertified_numerical" && (verifier.empty() || scope.empty() || evidence.empty()))
      throw std::invalid_argument("artifact certificate claim lacks verifier, scope or evidence");
  }
};
struct CertificateRequirement {
  std::string type,verifier,scope;
  bool accepts(const Certificate& c) const {
    return (type.empty() || type==c.type) && (verifier.empty() || verifier==c.verifier) &&
        (scope.empty() || scope==c.scope);
  }
};
struct Record {
  std::string semantic_id,content_id,payload_sha256;
  Demand guarantee;
  Certificate certificate;
  json::value payload;
  std::vector<std::string> parent_content_ids;
};

class Store {
 public:
  explicit Store(fs::path root,std::size_t maximum_bytes=64*1024*1024,
                 std::size_t maximum_candidates=1024)
      :root_(fs::absolute(std::move(root))),maximum_bytes_(maximum_bytes),maximum_candidates_(maximum_candidates) {
    if(!maximum_bytes || maximum_bytes>1024ULL*1024*1024 || !maximum_candidates || maximum_candidates>100000)
      throw std::invalid_argument("invalid artifact store limits");
    fs::create_directories(root_);detail::sync_directory(root_);detail::sync_directory(root_.parent_path());
  }
  fs::path path(const std::string& semantic_id,const std::string& content_id) const {
    if(!detail::hash_string(semantic_id) || !detail::hash_string(content_id))throw std::invalid_argument("invalid artifact content address");
    return root_/semantic_id/(content_id+".json");
  }
  Record put(const Identity& identity,const Demand& guarantee,const json::value& payload,
             const Certificate& certificate={},const std::vector<std::string>& parents={}) {
    auto scientific=identity.json_value();auto resources=detail::demand_json(guarantee);certificate.validate();
    if(parents.size()!=identity.parents.size())throw std::invalid_argument("artifact parent provenance count mismatch");
    for(std::size_t i=0;i<parents.size();++i) {
      if(!detail::hash_string(parents[i]))throw std::invalid_argument("invalid parent content identity");
      // Bind each concrete parent to the matching ordered scientific parent.
      read_address(identity.parents[i].semantic_id,parents[i]);
    }
    const auto semantic=identity.key(),payload_hash=detail::sha256(detail::canonical(payload));
    json::array parent_json;for(const auto& p:parents)parent_json.emplace_back(p);
    json::object claim{{"type",certificate.type},{"verifier",certificate.verifier},{"scope",certificate.scope},
      {"evidence",certificate.evidence},{"payload_sha256",payload_hash}};
    json::object envelope{{"schema","DiffExp3.ArtifactRecord/v1"},{"identity",scientific},
      {"semantic_id",semantic},{"guarantee",resources},{"payload",payload},{"payload_sha256",payload_hash},
      {"certificate",claim},{"parent_content_ids",parent_json}};
    const auto content=detail::sha256(detail::canonical(envelope));envelope.emplace("content_id",content);
    const auto bytes=detail::canonical(envelope);
    if(bytes.size()>maximum_bytes_)throw std::invalid_argument("artifact payload exceeds store size limit");
    const auto target=path(semantic,content);fs::create_directories(target.parent_path());detail::sync_directory(root_);
    atomic_publish(target,bytes);
    return read(identity,content);
  }
  Record read(const Identity& identity,const std::string& content_id) const {
    return read_address(identity.key(),content_id,&identity);
  }
  std::optional<Record> lookup(const Identity& identity,const Demand& demand,
                              const CertificateRequirement& certificate={}) const {
    detail::demand_json(demand);const auto directory=root_/identity.key();
    if(!fs::exists(directory))return std::nullopt;
    std::vector<std::string> candidates;
    for(const auto& entry:fs::directory_iterator(directory)) {
      if(entry.path().extension()!=".json")continue;
      if(!entry.is_regular_file() || fs::is_symlink(entry.symlink_status()))throw std::runtime_error("invalid artifact directory entry");
      candidates.push_back(entry.path().stem().string());
      if(candidates.size()>maximum_candidates_)throw std::runtime_error("artifact lookup candidate budget exhausted");
    }
    std::sort(candidates.begin(),candidates.end());std::optional<Record> best;
    for(const auto& candidate:candidates) {
      auto record=read(identity,candidate); // Corruption fails closed, even for a weaker resource record.
      if(record.guarantee.dominates(demand) && certificate.accepts(record.certificate) &&
          (!best || best->guarantee.dominates(record.guarantee)))best=std::move(record);
    }
    return best;
  }
 private:
  fs::path root_;std::size_t maximum_bytes_,maximum_candidates_;
  json::object read_json(const fs::path& file) const {
    if(fs::is_symlink(fs::symlink_status(file)) || !fs::is_regular_file(file))throw std::runtime_error("artifact is missing or is not a regular file");
    const auto size=fs::file_size(file);
    if(size>maximum_bytes_)throw std::runtime_error("artifact exceeds configured size limit");
    std::ifstream input(file,std::ios::binary);std::string bytes(size,'\0');input.read(bytes.data(),size);
    if(!input)throw std::runtime_error("cannot read complete artifact");
    json::parse_options options;options.max_depth=128;
    auto value=json::parse(bytes,{},options);
    if(detail::canonical(value)!=bytes)throw std::runtime_error("artifact is not in the canonical storage encoding");
    return value.as_object();
  }
  Record read_address(const std::string& semantic,const std::string& content,const Identity* expected=nullptr) const {
    auto envelope=read_json(path(semantic,content));
    detail::keys(envelope,{"schema","identity","semantic_id","guarantee","payload","payload_sha256","certificate","parent_content_ids","content_id"});
    if(detail::string(envelope.at("schema"))!="DiffExp3.ArtifactRecord/v1" ||
        detail::string(envelope.at("semantic_id"))!=semantic || detail::string(envelope.at("content_id"))!=content)
      throw std::runtime_error("artifact schema or content address mismatch");
    auto scientific=envelope.at("identity").as_object();
    detail::keys(scientific,{"schema","kind","algorithm_version","family","ordered_basis","normalization","branch","geometry","boundary","parents","scientific_inputs"});
    if(detail::string(scientific.at("schema"))!="DiffExp3.ArtifactIdentity/v1" || detail::sha256(detail::canonical(scientific))!=semantic)
      throw std::runtime_error("artifact scientific identity checksum mismatch");
    auto parsed_identity=Identity::from_json(scientific);
    if(expected && detail::canonical(scientific)!=detail::canonical(expected->json_value()))
      throw std::runtime_error("artifact scientific identity mismatch");
    envelope.erase("content_id");
    if(detail::sha256(detail::canonical(envelope))!=content)throw std::runtime_error("artifact envelope checksum mismatch");
    Record result;result.semantic_id=semantic;result.content_id=content;
    result.guarantee=detail::demand_from(envelope.at("guarantee"));result.payload=envelope.at("payload");
    result.payload_sha256=detail::string(envelope.at("payload_sha256"));
    if(detail::sha256(detail::canonical(result.payload))!=result.payload_sha256)throw std::runtime_error("artifact payload checksum mismatch");
    const auto& claim=envelope.at("certificate").as_object();
    detail::keys(claim,{"type","verifier","scope","evidence","payload_sha256"});
    if(detail::string(claim.at("payload_sha256"))!=result.payload_sha256)throw std::runtime_error("artifact certificate payload binding mismatch");
    result.certificate={detail::string(claim.at("type")),detail::string(claim.at("verifier")),detail::string(claim.at("scope")),claim.at("evidence").as_object()};
    result.certificate.validate();
    for(const auto& p:envelope.at("parent_content_ids").as_array()) {
      auto id=detail::string(p);if(!detail::hash_string(id))throw std::runtime_error("invalid stored parent content identity");
      result.parent_content_ids.push_back(std::move(id));
    }
    if(result.parent_content_ids.size()!=scientific.at("parents").as_array().size())throw std::runtime_error("artifact stored parent provenance count mismatch");
    return result;
  }
  void atomic_publish(const fs::path& target,const std::string& bytes) {
    auto name=(target.parent_path()/".publish-XXXXXX").string();std::vector<char> writable(name.begin(),name.end());writable.push_back('\0');
    int fd=::mkstemp(writable.data());if(fd<0)throw std::runtime_error("cannot create artifact staging file");
    const fs::path temporary(writable.data());
    try {
      std::size_t offset=0;
      while(offset<bytes.size()) {
        const auto n=::write(fd,bytes.data()+offset,std::min<std::size_t>(bytes.size()-offset,1024*1024));
        if(n<0 && errno==EINTR)continue;
        if(n<=0)throw std::runtime_error("cannot write artifact staging file");offset+=n;
      }
      if(::fchmod(fd,0444) || ::fsync(fd))throw std::runtime_error("cannot sync immutable artifact");
      if(::close(fd)) {fd=-1;throw std::runtime_error("cannot close artifact staging file");}fd=-1;
      // link is atomic and never replaces an existing immutable content file.
      if(::link(temporary.c_str(),target.c_str())!=0) {
        if(errno!=EEXIST)throw std::runtime_error("cannot atomically publish artifact");
        if(detail::canonical(read_json(target))!=bytes)throw std::runtime_error("immutable artifact content collision or corruption");
      }
      fs::remove(temporary);detail::sync_directory(target.parent_path());
    } catch(...) {if(fd>=0)::close(fd);std::error_code ignored;fs::remove(temporary,ignored);throw;}
  }
};
} // namespace diffexp::artifacts

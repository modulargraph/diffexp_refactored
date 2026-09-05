#pragma once

#include "diffexp/exact.hpp"
#include <map>

namespace diffexp {
// Exact affine maps of shared ball sources. Coordinates are materialized only
// at an explicit enclosure boundary; exact cancellation happens before Arb.
class LinearUncertainty {
 public:
  using Ball = kernel::ComplexBall;
  struct Source { std::string identity; std::vector<Ball> values; };
  using SourceRef = std::shared_ptr<const Source>;
  using Matrix = std::vector<std::vector<Rational>>;
  explicit LinearUncertainty(SourceRef source)
      : offset_(source ? source->values.size() : 0,Ball(0)) {
    if (!source || source->identity.empty() || source->values.empty())
      throw std::invalid_argument("invalid uncertainty source");
    Matrix m(offset_.size(),std::vector<Rational>(offset_.size(),Rational(0)));
    for (std::size_t i=0;i<m.size();++i) m[i][i]=Rational(1);
    sources_.emplace(source->identity,Term{source,std::move(m)});
  }
  LinearUncertainty transformed(const Matrix& map) const {
    if (map.empty()) throw std::invalid_argument("empty linear map");
    for (const auto& row:map)
      if (row.size()!=offset_.size()) throw std::invalid_argument("linear map shape mismatch");
    LinearUncertainty out(*this);
    out.offset_.assign(map.size(),Ball(0));
    for (std::size_t i=0;i<map.size();++i)
      for (std::size_t j=0;j<offset_.size();++j)
        if (!map[i][j].is_zero()) out.offset_[i]+=Ball::from_strings(map[i][j].str())*offset_[j];
    for (auto& [id,term]:out.sources_) {
      const auto& old=sources_.at(id).map;
      term.map.assign(map.size(),std::vector<Rational>(old[0].size(),Rational(0)));
      for (std::size_t i=0;i<map.size();++i)
        for (std::size_t j=0;j<map[i].size();++j)
          if (!map[i][j].is_zero())
            for (std::size_t k=0;k<old[j].size();++k) term.map[i][k]+=map[i][j]*old[j][k];
    }
    return out;
  }
  std::vector<Ball> independent_hull() const {
    auto out=offset_;
    for (const auto& [id,term]:sources_)
      for (std::size_t i=0;i<out.size();++i)
        for (std::size_t j=0;j<term.source->values.size();++j)
          if (!term.map[i][j].is_zero()) {
            const auto value=term.map[i][j]==Rational(1) ? term.source->values[j] :
                Ball::from_strings(term.map[i][j].str())*term.source->values[j];
            if (out[i].is_zero()) out[i]=value; else out[i]+=value;
          }
    return out;
  }
 private:
  struct Term { SourceRef source; Matrix map; };
  std::vector<Ball> offset_;
  std::map<std::string,Term> sources_;
};
}  // namespace diffexp

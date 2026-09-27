#pragma once
#include "diffexp/epsilon_gauge.hpp"
#include "diffexp/fuchsify.hpp"

namespace diffexp::recursion {
// J = to_basis * Y. This changes only the internal physical basis; public
// requested integrals and the immutable recursion graph retain their meaning.
struct PhysicalBasisTransform {
  ExactEpsilonMatrix original_connection,to_basis,from_basis,connection;
};
inline void verify_physical_basis(const PhysicalBasisTransform& basis,
    const ExactEpsilonMatrix& original,std::size_t xi) {
  const auto n=original.size();
  if(!n || n>256)throw std::invalid_argument("physical basis dimension budget");
  for(const auto* matrix:{&original,&basis.original_connection,&basis.to_basis,
      &basis.from_basis,&basis.connection}) {
    if(matrix->size()!=n)throw std::invalid_argument("physical basis row count");
    for(const auto& row:*matrix) {
      if(row.size()!=n)throw std::invalid_argument("physical basis column count");
      for(const auto& value:row)if(value.variables()!=original[0][0].variables())
        throw std::invalid_argument("physical basis field ordering");
    }
  }
  if(basis.original_connection!=original)
    throw std::invalid_argument("physical basis belongs to a different connection");
  const auto check_identity=[&](const ExactEpsilonMatrix& product) {
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)
      if(product[i][j]!=original[0][0].constant(i==j))
        throw std::invalid_argument("physical basis inverse identity failed");
  };
  check_identity(fuchsify::detail::multiply(basis.to_basis,basis.from_basis));
  check_identity(fuchsify::detail::multiply(basis.from_basis,basis.to_basis));
  const auto left=fuchsify::detail::multiply(basis.connection,basis.to_basis);
  const auto right=fuchsify::detail::multiply(basis.to_basis,original);
  for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)
    if(left[i][j]!=right[i][j]+basis.to_basis[i][j].derivative(xi))
      throw std::invalid_argument("physical basis differential identity failed");
}
} // namespace diffexp::recursion

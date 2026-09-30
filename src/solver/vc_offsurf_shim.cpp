// C entry point for off-surface virtual-casing evaluation.
// Upstream virtual-casing exposes ComputeBextOffSurf only as a C++ method;
// this mirrors the style of its own C wrappers in src/virtual-casing.cpp.
#include <virtual-casing.hpp>

#include <algorithm>
#include <vector>

extern "C" void VirtualCasingComputeBextOffSurfD(double* Bext, const double* B, long Nt, long Np,
                                                 const double* Xt, long n_points, const void* ctx) {
  const std::vector<double> B_(B, B + 3 * Nt * Np);
  const std::vector<double> Xt_(Xt, Xt + 3 * n_points);
  const auto& virtual_casing = *static_cast<const VirtualCasing<double>*>(ctx);
  const std::vector<double> Bext_ = virtual_casing.ComputeBextOffSurf(B_, Xt_);
  std::copy(Bext_.begin(), Bext_.end(), Bext);
}

#include <vector>

#include <virtual-casing.hpp>

extern "C" void VirtualCasingComputeBextOffSurfD(
    double* Bext,
    const double* B,
    const double* Xt,
    long src_Nt,
    long src_Np,
    long trg_Nt,
    long trg_Np,
    const void* ctx) {
  const auto* virtual_casing = reinterpret_cast<const VirtualCasing<double>*>(ctx);
  if (!virtual_casing) {
    return;
  }
  const std::vector<double> B_vec(B, B + 3 * src_Nt * src_Np);
  const std::vector<double> Xt_vec(Xt, Xt + 3 * trg_Nt * trg_Np);
  const std::vector<double> Bext_vec = virtual_casing->ComputeBextOffSurf(B_vec, Xt_vec);
  const auto count = static_cast<long>(Bext_vec.size());
  for (long i = 0; i < count; ++i) {
    Bext[i] = Bext_vec[static_cast<std::size_t>(i)];
  }
}

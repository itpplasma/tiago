// ISO C adapter between Tiago and VMEC++; see vmecpp_adapter.h.
//
// Adapted in part from VMEC++ (https://github.com/proximafusion/vmecpp,
// SPDX-License-Identifier: MIT, Copyright Proxima Fusion GmbH): the model
// bindings of src/vmecpp/cpp/vmecpp/vmec/pybind11/pybind_vmec.cc (VmecModel)
// and the implicit VJP of src/vmecpp/autodiff.py, ported to C++ with a block-tridiagonal
// LU instead of per-cotangent GMRES.
#include "vmecpp_adapter.h"

#include <Eigen/Dense>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <memory>
#include <numbers>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "vmecpp/common/vmec_indata/vmec_indata.h"
#include "vmecpp/vmec/geometry/geometry.h"
#include "vmecpp/vmec/geometry/vmec_geometry.h"
#include "vmecpp/vmec/output_quantities/output_quantities.h"
#include "vmecpp/vmec/vmec/vmec.h"

namespace {

std::string& ErrorMessage() {
  thread_local std::string message;
  return message;
}

int Fail(const std::string& message) {
  ErrorMessage() = message;
  return 1;
}

template <typename Function>
int Guard(Function&& function) {
  try {
    function();
    ErrorMessage().clear();
    return 0;
  } catch (const std::exception& error) {
    return Fail(error.what());
  } catch (...) {
    return Fail("unknown VMEC++ error");
  }
}

// ---- active state layout (pybind_vmec.cc) --------------------------------

std::vector<std::span<double>> ActiveSpans(vmecpp::FourierGeometry& x,
                                           const vmecpp::Sizes& s) {
  std::vector<std::span<double>> out = {x.rmncc};
  if (s.lthreed) out.push_back(x.rmnss);
  if (s.lasym) out.push_back(x.rmnsc);
  if (s.lasym && s.lthreed) out.push_back(x.rmncs);
  out.push_back(x.zmnsc);
  if (s.lthreed) out.push_back(x.zmncs);
  if (s.lasym) out.push_back(x.zmncc);
  if (s.lasym && s.lthreed) out.push_back(x.zmnss);
  out.push_back(x.lmnsc);
  if (s.lthreed) out.push_back(x.lmncs);
  if (s.lasym) out.push_back(x.lmncc);
  if (s.lasym && s.lthreed) out.push_back(x.lmnss);
  return out;
}

std::vector<std::span<double>> ActiveSpans(vmecpp::FourierForces& x,
                                           const vmecpp::Sizes& s) {
  std::vector<std::span<double>> out = {x.frcc};
  if (s.lthreed) out.push_back(x.frss);
  if (s.lasym) out.push_back(x.frsc);
  if (s.lasym && s.lthreed) out.push_back(x.frcs);
  out.push_back(x.fzsc);
  if (s.lthreed) out.push_back(x.fzcs);
  if (s.lasym) out.push_back(x.fzcc);
  if (s.lasym && s.lthreed) out.push_back(x.fzss);
  out.push_back(x.flsc);
  if (s.lthreed) out.push_back(x.flcs);
  if (s.lasym) out.push_back(x.flcc);
  if (s.lasym && s.lthreed) out.push_back(x.flss);
  return out;
}

template <typename FourierObject>
Eigen::VectorXd FlattenActive(FourierObject& x, const vmecpp::Sizes& s) {
  const auto spans = ActiveSpans(x, s);
  Eigen::Index total = 0;
  for (const auto& sp : spans) total += static_cast<Eigen::Index>(sp.size());
  Eigen::VectorXd out(total);
  Eigen::Index offset = 0;
  for (const auto& sp : spans) {
    const auto n = static_cast<Eigen::Index>(sp.size());
    out.segment(offset, n) = Eigen::Map<const Eigen::VectorXd>(sp.data(), n);
    offset += n;
  }
  return out;
}

template <typename FourierObject>
void UnflattenActive(FourierObject& x, const vmecpp::Sizes& s,
                     const Eigen::VectorXd& flat) {
  auto spans = ActiveSpans(x, s);
  Eigen::Index offset = 0;
  for (auto& sp : spans) {
    const auto n = static_cast<Eigen::Index>(sp.size());
    Eigen::Map<Eigen::VectorXd>(sp.data(), n) = flat.segment(offset, n);
    offset += n;
  }
  if (offset != flat.size()) {
    throw std::runtime_error("state vector has the wrong length");
  }
}

// ---- the model at a converged state (VmecModel in pybind_vmec.cc) --------

class Model {
 public:
  Model(const vmecpp::VmecINDATA& indata, int ns,
        const vmecpp::HotRestartState& initial_state) {
    auto vmec_or = vmecpp::Vmec::FromIndata(indata, nullptr, 1,
                                            vmecpp::OutputMode::kSilent);
    if (!vmec_or.ok()) {
      throw std::runtime_error(std::string(vmec_or.status().message()));
    }
    vmec_ = std::move(vmec_or.value());
    vmecpp::Vmec& v = *vmec_;
    v.always_fix_m1_gauge_ = true;
    v.fc_.ns_old = 0;
    v.fc_.delt0r = v.indata_.delt;
    v.fc_.ns_min = 3;
    v.fc_.nsval = ns;
    const Eigen::VectorXi& ns_array = v.indata_.ns_array;
    int index = static_cast<int>(ns_array.size()) - 1;
    for (int i = 0; i < ns_array.size(); ++i) {
      if (ns_array[i] == ns) {
        index = i;
        break;
      }
    }
    if (index >= 0 && index < v.indata_.ftol_array.size()) {
      v.fc_.ftolv = v.indata_.ftol_array[index];
    }
    if (index >= 0 && index < v.indata_.niter_array.size()) {
      v.fc_.niterv = v.indata_.niter_array[index];
    }
    double delt0 = v.indata_.delt;
    auto initialized =
        v.InitializeRadial(vmecpp::VmecCheckpoint::NONE, INT_MAX, ns, 0, delt0,
                           initial_state);
    if (!initialized.ok()) {
      throw std::runtime_error(std::string(initialized.status().message()));
    }
  }

  vmecpp::Vmec& vmec() { return *vmec_; }
  const vmecpp::Sizes& sizes() const { return vmec_->s_; }
  int ns() const { return vmec_->fc_.ns; }

  void Evaluate() {
    bool need_restart = false;
    std::string error_message;
    vmec_->fc_.restart_reason = vmecpp::RestartReason::NO_RESTART;
#ifdef _OPENMP
#pragma omp parallel num_threads(1)
#endif
    {
      auto s = vmec_->m_[0]->update(
          *vmec_->decomposed_x_[0], *vmec_->physical_x_[0],
          *vmec_->decomposed_f_[0], *vmec_->physical_f_[0], need_restart,
          last_preconditioner_update_, last_full_update_nestor_, vmec_->fc_, 2,
          2, vmecpp::VmecCheckpoint::NONE, INT_MAX, false, true);
      if (!s.ok()) error_message = std::string(s.status().message());
    }
    if (!error_message.empty()) throw std::runtime_error(error_message);
  }

  Eigen::VectorXd GetState() {
    return FlattenActive(*vmec_->decomposed_x_[0], vmec_->s_);
  }
  void SetState(const Eigen::VectorXd& flat) {
    UnflattenActive(*vmec_->decomposed_x_[0], vmec_->s_, flat);
    exact_primal_valid_ = false;
  }

#ifdef VMECPP_ENABLE_ENZYME
  void EnsurePrimal() {
    vmecpp::IdealMhdModel& model = *vmec_->m_[0];
    const int gs = static_cast<int>(model.r1_e.size());
    if (!exact_primal_valid_ ||
        exact_primal_.size() != static_cast<Eigen::Index>(20 * gs)) {
      exact_primal_.setZero(20 * gs);
      model.packGeometry(*vmec_->decomposed_x_[0], *vmec_->physical_x_[0],
                         exact_primal_.data(), gs, true);
      exact_primal_valid_ = true;
    }
  }

  Eigen::VectorXd HessianVectorProduct(const Eigen::VectorXd& v) {
    RequireNoLforbal();
    vmecpp::IdealMhdModel& model = *vmec_->m_[0];
    const int gs = static_cast<int>(model.r1_e.size());
    EnsurePrimal();
    Eigen::VectorXd dgeom = Eigen::VectorXd::Zero(20 * gs);
    vmec_->physical_x_backup_[0]->setZero();
    UnflattenActive(*vmec_->physical_x_backup_[0], vmec_->s_, v);
    model.packGeometry(*vmec_->physical_x_backup_[0], *vmec_->physical_x_[0],
                       dgeom.data(), gs, false);
    model.applyExactForceJacobian(exact_primal_.data(), dgeom.data(), gs,
                                  *vmec_->physical_f_[0],
                                  *vmec_->decomposed_f_[0], true);
    return FlattenActive(*vmec_->decomposed_f_[0], vmec_->s_);
  }

  Eigen::VectorXd HessianVectorProductTranspose(const Eigen::VectorXd& w) {
    RequireNoLforbal();
    vmecpp::IdealMhdModel& model = *vmec_->m_[0];
    const int gs = static_cast<int>(model.r1_e.size());
    EnsurePrimal();
    vmec_->decomposed_f_[0]->setZero();
    UnflattenActive(*vmec_->decomposed_f_[0], vmec_->s_, w);
    vmec_->physical_x_backup_[0]->setZero();
    model.applyExactForceJacobianTranspose(
        exact_primal_.data(), gs, *vmec_->decomposed_f_[0],
        *vmec_->physical_f_[0], *vmec_->physical_x_[0],
        *vmec_->physical_x_backup_[0], true);
    return FlattenActive(*vmec_->physical_x_backup_[0], vmec_->s_);
  }

  // Cotangents of the half-grid (mu0 p, iota, current) from a force cotangent.
  void ProfileVjp(const Eigen::VectorXd& force_bar, double* out) {
    RequireNoLforbal();
    if (vmec_->indata_.gamma != 0.0 || vmec_->indata_.lasym) {
      throw std::runtime_error("profile derivatives need gamma = 0 and lasym = false");
    }
    vmecpp::IdealMhdModel& model = *vmec_->m_[0];
    const int gs = static_cast<int>(model.r1_e.size());
    EnsurePrimal();
    const int n_half = vmec_->fc_.ns - 1;
    Eigen::VectorXd chip_bar = Eigen::VectorXd::Zero(n_half);
    Eigen::VectorXd pres_bar = Eigen::VectorXd::Zero(n_half);
    Eigen::VectorXd chip_profile_bar = Eigen::VectorXd::Zero(n_half);
    Eigen::VectorXd curr_bar = Eigen::VectorXd::Zero(n_half);
    vmec_->decomposed_f_[0]->setZero();
    UnflattenActive(*vmec_->decomposed_f_[0], vmec_->s_, force_bar);
    model.profileVjp(exact_primal_.data(), gs, *vmec_->decomposed_f_[0],
                     *vmec_->physical_f_[0], chip_bar.data(), pres_bar.data(),
                     chip_profile_bar.data(), curr_bar.data());
    const Eigen::VectorXd& phip_half = vmec_->p_[0]->phipH;
    const double flip = vmec_->fc_.haveToFlipTheta ? -1.0 : 1.0;
    for (int k = 0; k < n_half; ++k) {
      out[k] = pres_bar[k];
      out[n_half + k] = flip * chip_profile_bar[k] * phip_half[k];
      out[2 * n_half + k] = curr_bar[k];
    }
  }
#endif  // VMECPP_ENABLE_ENZYME

  // Transpose of the state-to-geometry map of MakeGeometry (stellarator-
  // symmetric, fixed boundary); the poloidal-flux cotangent is zero.
  Eigen::VectorXd GeometryStateVjp(const double* coefficient_bar) {
    const int ns = vmec_->fc_.ns;
    const int mpol = vmec_->s_.mpol;
    const int ntor = vmec_->s_.ntor;
    const int size = ns * mpol * (ntor + 1);
    vmecpp::FourierGeometry state_bar(&vmec_->s_, vmec_->r_[0].get(), ns);
    state_bar.setZero();
    const double lambda_scale = vmec_->constants_.lamscale;
    const Eigen::VectorXd& phip_full = vmec_->p_[0]->phipF;
    auto scale = [&](int j, int m, int n, bool lambda) {
      double value = (m == 0 ? 1.0 : std::numbers::sqrt2) *
                     (n == 0 ? 1.0 : std::numbers::sqrt2);
      if (lambda) {
        if (j >= phip_full.size() || phip_full[j] == 0.0) {
          throw std::runtime_error("invalid lambda scale");
        }
        value *= lambda_scale / phip_full[j];
      }
      return value;
    };
    auto add_block = [&](std::span<double> destination, int block, bool lambda,
                         int skipped_mode) {
      const int offset = block * size;
      for (int j = 0; j < ns; ++j) {
        for (int m = 0; m < mpol; ++m) {
          if (m == skipped_mode) continue;
          for (int n = 0; n <= ntor; ++n) {
            const int index = (j * mpol + m) * (ntor + 1) + n;
            destination[index] +=
                scale(j, m, n, lambda) * coefficient_bar[offset + index];
          }
        }
      }
    };
    const double sigma = -vmec_->indata_.signgs;
    auto add_m1_pair = [&](std::span<double> first, std::span<double> second,
                           int first_block, int second_block) {
      for (int j = 0; j < ns; ++j) {
        for (int n = 0; n <= ntor; ++n) {
          const int index = (j * mpol + 1) * (ntor + 1) + n;
          const double c = scale(j, 1, n, false);
          const double a = coefficient_bar[first_block * size + index];
          const double b = coefficient_bar[second_block * size + index];
          first[index] += c * (a + sigma * b);
          second[index] += c * (sigma * a - b);
        }
      }
    };
    add_block(state_bar.rmncc, 0, false, -1);
    if (vmec_->s_.lthreed) {
      if (mpol > 1) {
        add_block(state_bar.rmnss, 1, false, 1);
        add_block(state_bar.zmncs, 5, false, 1);
        add_m1_pair(state_bar.rmnss, state_bar.zmncs, 1, 5);
      } else {
        add_block(state_bar.rmnss, 1, false, -1);
      }
    }
    add_block(state_bar.zmnsc, 4, false, -1);
    if (vmec_->s_.lthreed && mpol == 1) add_block(state_bar.zmncs, 5, false, -1);
    add_block(state_bar.lmnsc, 8, true, -1);
    if (vmec_->s_.lthreed) add_block(state_bar.lmncs, 9, true, -1);
    return FlattenActive(state_bar, vmec_->s_);
  }

 private:
  void RequireNoLforbal() const {
    if (vmec_->m_[0]->lforbal) {
      throw std::runtime_error("exact derivatives need lforbal = false");
    }
  }

  std::unique_ptr<vmecpp::Vmec> vmec_;
  int last_preconditioner_update_ = 0;
  int last_full_update_nestor_ = 0;
  bool exact_primal_valid_ = false;
  Eigen::VectorXd exact_primal_;
};

// ---- index sets of the implicit adjoint (autodiff.py) --------------------

struct Layout {
  int ns, mpol, ntor, modes;  // modes = mpol * (ntor + 1) per surface
  bool lthreed;
  int span;                   // ns * modes
  // offsets of the active spans, -1 if absent
  int r_cc, r_ss, z_sc, z_cs, l_sc, l_cs;
};

Layout MakeLayout(const vmecpp::Sizes& s, int ns) {
  Layout l{};
  l.ns = ns;
  l.mpol = s.mpol;
  l.ntor = s.ntor;
  l.modes = s.mpol * (s.ntor + 1);
  l.lthreed = s.lthreed;
  l.span = ns * l.modes;
  int k = 0;
  l.r_cc = k++ * l.span;
  l.r_ss = s.lthreed ? k++ * l.span : -1;
  l.z_sc = k++ * l.span;
  l.z_cs = s.lthreed ? k++ * l.span : -1;
  l.l_sc = k++ * l.span;
  l.l_cs = s.lthreed ? k++ * l.span : -1;
  return l;
}

std::vector<int> BoundaryEntries(const Layout& l) {
  std::vector<int> out;
  const int edge = (l.ns - 1) * l.modes;
  for (int start : {l.r_cc, l.r_ss, l.z_sc, l.z_cs}) {
    if (start < 0) continue;
    for (int i = start + edge; i < start + l.span; ++i) out.push_back(i);
  }
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<int> GaugeEntries(const Layout& l, int mpol_geometry,
                              int ntor_geometry) {
  std::vector<int> out;
  if (!l.lthreed || mpol_geometry <= 1) return out;
  for (int j = 1; j < l.ns - 1; ++j) {
    for (int n = 1; n <= ntor_geometry; ++n) {
      out.push_back(l.z_cs + j * l.modes + 1 * (l.ntor + 1) + n);
    }
  }
  return out;
}

}  // namespace

// ---- handle ---------------------------------------------------------------

// Block-tridiagonal linear system A x = b, one block per radial surface
// (VMEC's forces on surface j depend on the geometry of surfaces j - 1, j and
// j + 1 only). Block LU without pivoting across blocks (block Thomas), partial
// pivoting within: O(ns m^3) time and O(ns m^2) memory for blocks of size m,
// against O((ns m)^3) and O((ns m)^2) for a dense LU.
class BlockTridiagonal {
 public:
  // diagonal[k] = A(k, k), lower[k] = A(k, k - 1) (k >= 1),
  // upper[k] = A(k, k + 1) (k + 1 < blocks); blocks may be empty.
  void Factor(std::vector<Eigen::MatrixXd> diagonal,
              std::vector<Eigen::MatrixXd> lower,
              std::vector<Eigen::MatrixXd> upper) {
    const int blocks = static_cast<int>(diagonal.size());
    lu_.assign(blocks, {});
    gain_.assign(blocks, {});
    lower_ = std::move(lower);
    for (int k = 0; k < blocks; ++k) {
      Eigen::MatrixXd schur = std::move(diagonal[k]);
      if (k > 0 && schur.size() > 0 && gain_[k - 1].size() > 0) {
        schur.noalias() -= lower_[k] * gain_[k - 1];
      }
      if (schur.rows() == 0) continue;
      lu_[k].compute(schur);
      if (k + 1 < blocks) gain_[k] = lu_[k].solve(upper[k]);  // S_k^-1 A(k, k+1)
    }
  }

  // b and x are concatenated block by block.
  Eigen::VectorXd Solve(const Eigen::VectorXd& b) const {
    const int blocks = static_cast<int>(lu_.size());
    std::vector<int> offset(blocks + 1, 0);
    for (int k = 0; k < blocks; ++k) offset[k + 1] = offset[k] + Size(k);
    Eigen::VectorXd x(b.size());
    for (int k = 0; k < blocks; ++k) {  // forward: y_k = S_k^-1 (b_k - L_k y_k-1)
      if (Size(k) == 0) continue;
      Eigen::VectorXd r = b.segment(offset[k], Size(k));
      if (k > 0 && Size(k - 1) > 0) r.noalias() -= lower_[k] * x.segment(offset[k - 1], Size(k - 1));
      x.segment(offset[k], Size(k)) = lu_[k].solve(r);
    }
    for (int k = blocks - 2; k >= 0; --k) {  // backward: x_k = y_k - G_k x_k+1
      if (Size(k) == 0 || Size(k + 1) == 0) continue;
      x.segment(offset[k], Size(k)).noalias() -= gain_[k] * x.segment(offset[k + 1], Size(k + 1));
    }
    return x;
  }

 private:
  int Size(int k) const { return static_cast<int>(lu_[k].rows()); }
  std::vector<Eigen::PartialPivLU<Eigen::MatrixXd>> lu_;
  std::vector<Eigen::MatrixXd> lower_, gain_;
};

struct tiago_vmecpp {
  vmecpp::VmecINDATA indata;
  std::optional<vmecpp::OutputQuantities> output;
  std::unique_ptr<Model> model;
  // last converged equilibrium, the start of the next solve (hot restart)
  std::optional<vmecpp::WOutFileContents> last_wout;
  bool hot_restart = true;
  int status = -1;
  // adjoint cache for the current solve
  bool factorized = false;
  std::vector<int> interior, prescribed, gauge;  // interior: grouped by surface
  BlockTridiagonal system;                        // transposed interior Hessian
};

namespace {

Eigen::VectorXd& Vector(vmecpp::VmecINDATA& in, const std::string& name) {
  if (name == "am") return in.am;
  if (name == "ai") return in.ai;
  if (name == "ac") return in.ac;
  if (name == "aphi") return in.aphi;
  throw std::runtime_error("unknown input array: " + name);
}

double* Scalar(vmecpp::VmecINDATA& in, const std::string& name) {
  if (name == "pres_scale") return &in.pres_scale;
  if (name == "curtor") return &in.curtor;
  if (name == "phiedge") return &in.phiedge;
  if (name == "bloat") return &in.bloat;
  if (name == "spres_ped") return &in.spres_ped;
  if (name == "gamma") return &in.gamma;
  return nullptr;
}

const vmecpp::WOutFileContents& Wout(const tiago_vmecpp* h) {
  if (!h->output) throw std::runtime_error("no equilibrium solved yet");
  return h->output->wout;
}

void CopyVector(const Eigen::VectorXd& v, double* out, int length) {
  if (v.size() != length) {
    throw std::runtime_error("wrong length " + std::to_string(length) +
                             " (have " + std::to_string(v.size()) + ")");
  }
  std::copy(v.data(), v.data() + length, out);
}

}  // namespace

extern "C" const char* tiago_vmecpp_error(void) {
  return ErrorMessage().c_str();
}

extern "C" int tiago_vmecpp_create(const char* json_path, tiago_vmecpp** output) {
  if (json_path == nullptr || output == nullptr) return Fail("null argument");
  *output = nullptr;
  return Guard([&] {
    auto handle = std::make_unique<tiago_vmecpp>();
    handle->indata = vmecpp::VmecINDATA::FromFile(json_path);
    if (handle->indata.lfreeb) {
      throw std::runtime_error("Tiago runs VMEC++ fixed-boundary (lfreeb = false)");
    }
    if (handle->indata.lasym) {
      throw std::runtime_error("the VMEC++ adjoint needs lasym = false");
    }
    *output = handle.release();
  });
}

extern "C" void tiago_vmecpp_destroy(tiago_vmecpp* handle) { delete handle; }

extern "C" int tiago_vmecpp_get_input(const tiago_vmecpp* handle, const char* name,
                                      int i, int j, double* value) {
  return Guard([&] {
    auto& in = const_cast<vmecpp::VmecINDATA&>(handle->indata);
    const std::string key(name);
    if (double* s = Scalar(in, key)) {
      *value = *s;
    } else if (key == "rbc" || key == "zbs") {
      const auto& m = key == "rbc" ? in.rbc : in.zbs;
      if (i < 1 || i > m.rows() || j < 1 || j > m.cols()) {
        throw std::runtime_error("boundary index out of range");
      }
      *value = m(i - 1, j - 1);
    } else {
      const Eigen::VectorXd& v = Vector(in, key);
      *value = (i >= 1 && i <= v.size()) ? v[i - 1] : 0.0;
    }
  });
}

extern "C" int tiago_vmecpp_set_input(tiago_vmecpp* handle, const char* name, int i,
                                      int j, double value) {
  return Guard([&] {
    auto& in = handle->indata;
    const std::string key(name);
    if (key == "ftol") {
      in.ftol_array.setConstant(value);  // every multi-grid step
    } else if (key == "niter") {
      in.niter_array.setConstant(static_cast<int>(value));
    } else if (key == "hot_restart") {
      handle->hot_restart = value != 0.0;
    } else if (double* s = Scalar(in, key)) {
      *s = value;
    } else if (key == "rbc" || key == "zbs") {
      auto& m = key == "rbc" ? in.rbc : in.zbs;
      if (i < 1 || i > m.rows() || j < 1 || j > m.cols()) {
        throw std::runtime_error("boundary index out of range");
      }
      m(i - 1, j - 1) = value;
    } else {
      Eigen::VectorXd& v = Vector(in, key);
      if (i < 1) throw std::runtime_error("array index out of range");
      if (i > v.size()) {
        const auto old = v.size();
        v.conservativeResize(i);
        for (Eigen::Index k = old; k < i; ++k) v[k] = 0.0;
      }
      v[i - 1] = value;
    }
  });
}

extern "C" int tiago_vmecpp_input_length(const tiago_vmecpp* handle,
                                         const char* name, int* length) {
  return Guard([&] {
    auto& in = const_cast<vmecpp::VmecINDATA&>(handle->indata);
    *length = static_cast<int>(Vector(in, name).size());
  });
}

extern "C" int tiago_vmecpp_get_int(const tiago_vmecpp* handle, const char* name,
                                    int* value) {
  return Guard([&] {
    const auto& in = handle->indata;
    const std::string key(name);
    if (key == "ns") {
      *value = in.ns_array[in.ns_array.size() - 1];
    } else if (key == "mpol") {
      *value = in.mpol;
    } else if (key == "ntor") {
      *value = in.ntor;
    } else if (key == "nfp") {
      *value = in.nfp;
    } else if (key == "ncurr") {
      *value = in.ncurr;
    } else if (key == "lasym") {
      *value = in.lasym ? 1 : 0;
    } else if (key == "ntheta") {
      *value = in.ntheta;
    } else if (key == "nzeta") {
      *value = in.nzeta;
    } else if (key == "lfreeb") {
      *value = in.lfreeb ? 1 : 0;
    } else if (key == "status") {
      *value = handle->status;
    } else if (key == "signgs") {
      *value = Wout(handle).signgs;
    } else if (key == "mnmax") {
      *value = Wout(handle).mnmax;
    } else if (key == "mnmax_nyq") {
      *value = Wout(handle).mnmax_nyq;
    } else if (key == "have_to_flip_theta") {
      if (!handle->model) throw std::runtime_error("no equilibrium solved yet");
      *value = handle->model->vmec().fc_.haveToFlipTheta ? 1 : 0;
    } else {
      throw std::runtime_error("unknown integer: " + key);
    }
  });
}

extern "C" int tiago_vmecpp_get_string(const tiago_vmecpp* handle, const char* name,
                                       char* buffer, int length) {
  return Guard([&] {
    const std::string key(name);
    std::string value;
    if (key == "pmass_type") {
      value = handle->indata.pmass_type;
    } else if (key == "piota_type") {
      value = handle->indata.piota_type;
    } else if (key == "pcurr_type") {
      value = handle->indata.pcurr_type;
    } else {
      throw std::runtime_error("unknown string: " + key);
    }
    std::memset(buffer, ' ', static_cast<size_t>(length));
    std::memcpy(buffer, value.data(),
                std::min(value.size(), static_cast<size_t>(length)));
  });
}

extern "C" int tiago_vmecpp_solve(tiago_vmecpp* handle, int* converged) {
  *converged = 0;
  return Guard([&] {
    handle->output.reset();
    handle->model.reset();
    handle->factorized = false;
    handle->status = 1;
    const int last = static_cast<int>(handle->indata.ns_array.size()) - 1;
    const int ns = handle->indata.ns_array[last];
    // Hot restart: a single grid at the final ns from the last converged
    // equilibrium (inner surfaces; the boundary comes from the input), which
    // is close for the small parameter steps of a fit. Cold multigrid start
    // otherwise, and as fallback.
    vmecpp::VmecINDATA indata = handle->indata;
    std::optional<vmecpp::HotRestartState> initial;
    if (handle->hot_restart && handle->last_wout && handle->last_wout->ns == ns) {
      indata.ns_array = Eigen::VectorXi::Constant(1, ns);
      indata.ftol_array = Eigen::VectorXd::Constant(1, handle->indata.ftol_array[last]);
      indata.niter_array = Eigen::VectorXi::Constant(1, handle->indata.niter_array[last]);
      initial.emplace(*handle->last_wout, indata);
    }
    auto result = vmecpp::run(indata, initial, 1, vmecpp::OutputMode::kSilent, nullptr, true);
    if (!result.ok() && initial) {
      indata = handle->indata;
      result = vmecpp::run(indata, std::nullopt, 1, vmecpp::OutputMode::kSilent, nullptr, true);
    }
    if (!result.ok()) {
      ErrorMessage() = std::string(result.status().message());
      return;
    }
    handle->output = std::move(*result);
    handle->last_wout = handle->output->wout;
    vmecpp::HotRestartState restart(handle->output->wout, indata);
    handle->model = std::make_unique<Model>(indata, ns, restart);
    handle->status = 0;
    *converged = 1;
  });
}

extern "C" int tiago_vmecpp_geometry(const tiago_vmecpp* handle, double* coefficients,
                                     double* toroidal_flux, double* poloidal_flux) {
  return Guard([&] {
    if (!handle->output) throw std::runtime_error("no equilibrium solved yet");
    const vmecpp::Geometry g = vmecpp::MakeGeometry(
        handle->output->indata, handle->output->vmec_internal_results,
        vmecpp::GeometryCoefficientState::kPhysical);
    const auto& d = g.dimensions;
    const size_t size = static_cast<size_t>(d.ns) * d.mpol * (d.ntor + 1);
    const std::vector<double>* blocks[12] = {
        &g.coefficients.r_cc,      &g.coefficients.r_ss,
        &g.coefficients.r_sc,      &g.coefficients.r_cs,
        &g.coefficients.z_sc,      &g.coefficients.z_cs,
        &g.coefficients.z_cc,      &g.coefficients.z_ss,
        &g.coefficients.lambda_sc, &g.coefficients.lambda_cs,
        &g.coefficients.lambda_cc, &g.coefficients.lambda_ss};
    for (int b = 0; b < 12; ++b) {
      double* out = coefficients + b * size;
      if (blocks[b]->size() == size) {
        std::copy(blocks[b]->begin(), blocks[b]->end(), out);
      } else {
        std::fill(out, out + size, 0.0);
      }
    }
    std::copy(g.toroidal_flux.begin(), g.toroidal_flux.end(), toroidal_flux);
    std::copy(g.poloidal_flux.begin(), g.poloidal_flux.end(), poloidal_flux);
  });
}

extern "C" int tiago_vmecpp_radial(const tiago_vmecpp* handle, double* phip_full,
                                   double* phip_half, double* iota_half,
                                   double* current_half, double* mass_half,
                                   double* lamscale) {
  return Guard([&] {
    if (!handle->model) throw std::runtime_error("no equilibrium solved yet");
    auto& v = handle->model->vmec();
    const int ns = v.fc_.ns;
    const auto& p = *v.p_[0];
    CopyVector(p.phipF.head(ns), phip_full, ns);
    CopyVector(p.phipH.head(ns - 1), phip_half, ns - 1);
    CopyVector(p.iotaH.head(ns - 1), iota_half, ns - 1);
    CopyVector(p.currH.head(ns - 1), current_half, ns - 1);
    CopyVector(p.massH.head(ns - 1), mass_half, ns - 1);
    *lamscale = v.constants_.lamscale;
  });
}

extern "C" int tiago_vmecpp_wout(const tiago_vmecpp* handle, const char* name,
                                 double* values, int length) {
  return Guard([&] {
    const auto& w = Wout(handle);
    const std::string key(name);
    auto ints = [&](const Eigen::VectorXi& v) {
      if (v.size() != length) throw std::runtime_error("wrong length for " + key);
      for (int k = 0; k < length; ++k) values[k] = v[k];
    };
    auto matrix = [&](const vmecpp::RowMatrixXd& m, int modes) {
      // stored with the mode index fastest, radius slowest
      if (static_cast<int>(m.size()) != length) {
        throw std::runtime_error("wrong length for " + key);
      }
      const bool radius_rows = m.rows() == w.ns && m.cols() == modes;
      for (int j = 0; j < w.ns; ++j) {
        for (int k = 0; k < modes; ++k) {
          values[j * modes + k] = radius_rows ? m(j, k) : m(k, j);
        }
      }
    };
    if (key == "phi") CopyVector(w.phi, values, length);
    else if (key == "presf") CopyVector(w.presf, values, length);
    else if (key == "iotaf") CopyVector(w.iotaf, values, length);
    else if (key == "iotas") CopyVector(w.iotas, values, length);
    else if (key == "jcurv") CopyVector(w.jcurv, values, length);
    else if (key == "jcuru") CopyVector(w.jcuru, values, length);
    else if (key == "buco") CopyVector(w.buco, values, length);
    else if (key == "chi") CopyVector(w.chi, values, length);
    else if (key == "mass") CopyVector(w.mass, values, length);
    else if (key == "xm") ints(w.xm);
    else if (key == "xn") ints(w.xn);
    else if (key == "xm_nyq") ints(w.xm_nyq);
    else if (key == "xn_nyq") ints(w.xn_nyq);
    else if (key == "rmnc") matrix(w.rmnc, w.mnmax);
    else if (key == "zmns") matrix(w.zmns, w.mnmax);
    else if (key == "bsupumnc") matrix(w.bsupumnc, w.mnmax_nyq);
    else if (key == "bsupvmnc") matrix(w.bsupvmnc, w.mnmax_nyq);
    else if (key == "bmnc") matrix(w.bmnc, w.mnmax_nyq);
    else if (key == "ctor") values[0] = w.ctor;
    else if (key == "b0") values[0] = w.b0;
    else if (key == "Aminor_p") values[0] = w.Aminor_p;
    else if (key == "Rmajor_p") values[0] = w.Rmajor_p;
    else if (key == "phiedge") values[0] = handle->indata.phiedge;
    else throw std::runtime_error("unknown wout field: " + key);
  });
}

extern "C" int tiago_vmecpp_adjoint(tiago_vmecpp* handle, int ncot,
                                    const double* geometry_bar, double* boundary_bar,
                                    double* profile_bar, int max_block) {
  return Guard([&] {
#ifndef VMECPP_ENABLE_ENZYME
    (void)ncot; (void)geometry_bar; (void)boundary_bar; (void)profile_bar;
    (void)max_block;
    throw std::runtime_error(
        "this VMEC++ build has no exact force Jacobian (VMECPP_ENABLE_ENZYME)");
#else
    if (!handle->model) throw std::runtime_error("no equilibrium solved yet");
    Model& model = *handle->model;
    vmecpp::Vmec& v = model.vmec();
    const int ns = v.fc_.ns;
    const int mpol = v.s_.mpol;
    const int ntor = v.s_.ntor;
    const Layout layout = MakeLayout(v.s_, ns);
    const Eigen::VectorXd state = model.GetState();
    const int state_size = static_cast<int>(state.size());

    if (!handle->factorized) {
      model.SetState(state);
      model.Evaluate();
      std::vector<int> boundary = BoundaryEntries(layout);
      handle->gauge = GaugeEntries(layout, v.s_.mpolGeometry, v.s_.ntorGeometry);
      std::vector<char> excluded(state_size, 0);
      for (int i : boundary) excluded[i] = 1;
      for (int i : handle->gauge) excluded[i] = 1;
      // structural null space: entries with zero Hessian row or column
      std::mt19937_64 generator(0);
      std::normal_distribution<double> normal;
      Eigen::VectorXd column = Eigen::VectorXd::Zero(state_size);
      Eigen::VectorXd row = Eigen::VectorXd::Zero(state_size);
      for (int probe = 0; probe < 6; ++probe) {
        Eigen::VectorXd p(state_size);
        for (int i = 0; i < state_size; ++i) p[i] = normal(generator);
        column = column.cwiseMax(model.HessianVectorProductTranspose(p).cwiseAbs());
        row = row.cwiseMax(model.HessianVectorProduct(p).cwiseAbs());
      }
      const double threshold =
          1.0e-9 * std::max({column.maxCoeff(), row.maxCoeff(), 1.0});
      // interior unknowns grouped by radial surface, in state order within one
      std::vector<std::vector<int>> by_surface(ns);
      auto surface_of = [&](int i) { return (i % layout.span) / layout.modes; };
      for (int i = 0; i < state_size; ++i) {
        if (!excluded[i] && column[i] > threshold && row[i] > threshold) {
          by_surface[surface_of(i)].push_back(i);
        }
      }
      handle->interior.clear();
      std::vector<int> position(state_size, -1);  // index within its surface block
      int largest = 0;
      for (int j = 0; j < ns; ++j) {
        for (int k = 0; k < static_cast<int>(by_surface[j].size()); ++k) {
          position[by_surface[j][k]] = k;
          handle->interior.push_back(by_surface[j][k]);
        }
        largest = std::max(largest, static_cast<int>(by_surface[j].size()));
      }
      handle->prescribed = boundary;
      handle->prescribed.insert(handle->prescribed.end(), handle->gauge.begin(),
                                handle->gauge.end());
      if (handle->interior.empty()) throw std::runtime_error("the interior force operator is null");
      if (largest > max_block) {
        throw std::runtime_error("a surface block of " + std::to_string(largest) +
                                 " unknowns exceeds max_block = " + std::to_string(max_block));
      }
      // Blocks of A = H^T by colored probing with forward-mode Hessian-vector
      // products (about twice as fast as the reverse ones): one per (state
      // slot, surface class j mod 3) sets that slot on every third surface;
      // each result row belongs to exactly one probed column because columns
      // of surface j only reach rows of surfaces j - 1 .. j + 1.
      std::vector<Eigen::MatrixXd> diagonal(ns), lower(ns), upper(ns);
      auto size_of = [&](int j) { return static_cast<int>(by_surface[j].size()); };
      for (int j = 0; j < ns; ++j) {
        diagonal[j] = Eigen::MatrixXd::Zero(size_of(j), size_of(j));
        if (j > 0) lower[j] = Eigen::MatrixXd::Zero(size_of(j), size_of(j - 1));
        if (j + 1 < ns) upper[j] = Eigen::MatrixXd::Zero(size_of(j), size_of(j + 1));
      }
      std::vector<int> starts;
      for (int start : {layout.r_cc, layout.r_ss, layout.z_sc, layout.z_cs, layout.l_sc,
                        layout.l_cs}) {
        if (start >= 0) starts.push_back(start);
      }
      Eigen::VectorXd probe = Eigen::VectorXd::Zero(state_size);
      for (int start : starts) {
        for (int mode = 0; mode < layout.modes; ++mode) {
          for (int color = 0; color < 3; ++color) {
            std::vector<int> columns;
            for (int j = color; j < ns; j += 3) {
              const int i = start + j * layout.modes + mode;
              if (position[i] >= 0) columns.push_back(i);
            }
            if (columns.empty()) continue;
            for (int i : columns) probe[i] = 1.0;
            const Eigen::VectorXd result = model.HessianVectorProduct(probe);
            for (int i : columns) probe[i] = 0.0;
            // column i of H (surface j) holds row i of A = H^T in the blocks
            // (j, j - 1), (j, j) and (j, j + 1)
            for (int i : columns) {
              const int j = surface_of(i), c = position[i];
              for (int r = 0; r < size_of(j); ++r) diagonal[j](c, r) = result[by_surface[j][r]];
              if (j > 0) {
                for (int r = 0; r < size_of(j - 1); ++r) lower[j](c, r) = result[by_surface[j - 1][r]];
              }
              if (j + 1 < ns) {
                for (int r = 0; r < size_of(j + 1); ++r) upper[j](c, r) = result[by_surface[j + 1][r]];
              }
            }
          }
        }
      }
      handle->system.Factor(std::move(diagonal), std::move(lower), std::move(upper));
      handle->factorized = true;
    }

    const int coefficient_size = 12 * ns * mpol * (ntor + 1);
    const int n = static_cast<int>(handle->interior.size());
    const int boundary_size = 2 * mpol * (2 * ntor + 1);
    const int modes = layout.modes;
    const int edge = ns - 1;
    const int mpol_geometry = v.s_.mpolGeometry;
    const int ntor_geometry = v.s_.ntorGeometry;
    for (int c = 0; c < ncot; ++c) {
      // The Enzyme passes use the model's work buffers (and the cached primal)
      // as scratch; re-evaluate at the converged state before every cotangent.
      model.SetState(state);
      model.Evaluate();
      const Eigen::VectorXd state_bar =
          model.GeometryStateVjp(geometry_bar + static_cast<size_t>(c) * coefficient_size);
      Eigen::VectorXd rhs(n);
      for (int r = 0; r < n; ++r) rhs[r] = state_bar[handle->interior[r]];
      const Eigen::VectorXd adjoint = handle->system.Solve(rhs);
      Eigen::VectorXd embedded = Eigen::VectorXd::Zero(state_size);
      for (int r = 0; r < n; ++r) embedded[handle->interior[r]] = adjoint[r];
      const Eigen::VectorXd coupling = model.HessianVectorProductTranspose(embedded);
      Eigen::VectorXd full = Eigen::VectorXd::Zero(state_size);
      for (int i : handle->prescribed) full[i] = state_bar[i] - coupling[i];
      // fold the pinned gauge back onto the boundary entries it derives from
      for (int i : handle->gauge) {
        const int surface = (i - layout.z_cs) / modes;
        const int mode = (i - layout.z_cs) % modes;
        full[layout.z_cs + edge * modes + mode] +=
            std::sqrt(surface / (ns - 1.0)) * full[i];
        full[i] = 0.0;
      }
      // transpose of the boundary parser (Boundaries::ensureM1Constrained etc.)
      auto value = [&](int start, int m, int nn) {
        return start < 0 ? 0.0 : full[start + edge * modes + m * (ntor + 1) + nn];
      };
      Eigen::MatrixXd rbcc = Eigen::MatrixXd::Zero(mpol, ntor + 1);
      Eigen::MatrixXd zbsc = rbcc, rbss = rbcc, zbcs = rbcc;
      for (int m = 0; m < mpol; ++m) {
        for (int nn = 0; nn <= ntor; ++nn) {
          if (m >= mpol_geometry || nn > ntor_geometry) continue;
          const double scale = (m == 0 ? 1.0 : std::numbers::sqrt2) *
                               (nn == 0 ? 1.0 : std::numbers::sqrt2);
          rbcc(m, nn) = value(layout.r_cc, m, nn) / scale;
          zbsc(m, nn) = value(layout.z_sc, m, nn) / scale;
          if (layout.lthreed) {
            rbss(m, nn) = value(layout.r_ss, m, nn) / scale;
            zbcs(m, nn) = value(layout.z_cs, m, nn) / scale;
          }
        }
      }
      if (layout.lthreed && mpol > 1) {
        for (int nn = 0; nn <= ntor; ++nn) {
          const double r_bar = rbss(1, nn), z_bar = zbcs(1, nn);
          rbss(1, nn) = 0.5 * (r_bar + z_bar);
          zbcs(1, nn) = 0.5 * (r_bar - z_bar);
        }
      }
      if (v.fc_.haveToFlipTheta) {
        for (int m = 1; m < mpol; ++m) {
          const double parity = m % 2 == 0 ? 1.0 : -1.0;
          rbcc.row(m) *= parity;
          zbsc.row(m) *= -parity;
          if (layout.lthreed) {
            rbss.row(m) *= -parity;
            zbcs.row(m) *= parity;
          }
        }
      }
      double* out = boundary_bar + static_cast<size_t>(c) * boundary_size;
      std::fill(out, out + boundary_size, 0.0);
      auto at = [&](int block, int m, int signed_n) -> double& {
        return out[(block * mpol + m) * (2 * ntor + 1) + signed_n + ntor];
      };
      for (int m = 0; m < mpol; ++m) {
        for (int sn = -ntor; sn <= ntor; ++sn) {
          const int target = std::abs(sn);
          const double sign = sn > 0 ? 1.0 : (sn < 0 ? -1.0 : 0.0);
          at(0, m, sn) += rbcc(m, target);
          if (layout.lthreed && m > 0) at(0, m, sn) += sign * rbss(m, target);
          if (m > 0) at(1, m, sn) += zbsc(m, target);
          if (layout.lthreed) at(1, m, sn) -= sign * zbcs(m, target);
        }
      }
      model.ProfileVjp(-embedded, profile_bar + static_cast<size_t>(c) * 3 * (ns - 1));
    }
#endif
  });
}

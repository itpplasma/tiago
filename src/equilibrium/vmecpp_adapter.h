/* ISO C interface between Tiago (Fortran) and VMEC++ (C++).
 *
 * Fixed-boundary equilibria with VMEC's m = 1 gauge pinned, and the implicit
 * adjoint of the converged state: for cotangents on the output geometry
 * (VMEC++'s product-basis coefficients) it returns the cotangents on the input
 * boundary (rbc, zbs) and on the half-grid (mu0 p, iota, current) profiles.
 *
 * The adjoint needs VMEC++'s exact force Jacobian (an Enzyme build). Parts of
 * this adapter are adapted from VMEC++ (proximafusion/vmecpp, MIT licence):
 * the VmecModel bindings in pybind_vmec.cc and the implicit VJP in
 * autodiff.py; upstream maintains only the C++ and Python APIs.
 */
#ifndef TIAGO_VMECPP_ADAPTER_H_
#define TIAGO_VMECPP_ADAPTER_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct tiago_vmecpp tiago_vmecpp;

/* All functions return 0 on success; tiago_vmecpp_error() describes a failure. */
const char* tiago_vmecpp_error(void);

/* Create from a VMEC++ JSON input file. */
int tiago_vmecpp_create(const char* json_path, tiago_vmecpp** output);
void tiago_vmecpp_destroy(tiago_vmecpp* handle);

/* Input values by name: scalars (pres_scale, curtor, phiedge, bloat,
 * spres_ped, gamma), 1-based index i into arrays (am, ai, ac, aphi: i = k + 1
 * for coefficient k), (i, j) = (m + 1, n + ntor + 1) for rbc and zbs. Setting
 * an array element beyond its length extends the array with zeros. Setting
 * "ftol" or "niter" sets every entry of ftol_array or niter_array. */
int tiago_vmecpp_get_input(const tiago_vmecpp* handle, const char* name, int i,
                           int j, double* value);
int tiago_vmecpp_set_input(tiago_vmecpp* handle, const char* name, int i, int j,
                           double value);
int tiago_vmecpp_input_length(const tiago_vmecpp* handle, const char* name,
                              int* length);

/* Integers of the input: ns (last ns_array entry), mpol, ntor, nfp, ncurr,
 * lasym, signgs, ntheta, nzeta, lfreeb, and after a solve mnmax, mnmax_nyq,
 * have_to_flip_theta (0/1), status of the last solve (0 = converged). */
int tiago_vmecpp_get_int(const tiago_vmecpp* handle, const char* name,
                         int* value);
/* Profile types (pmass_type, piota_type, pcurr_type) into buffer. */
int tiago_vmecpp_get_string(const tiago_vmecpp* handle, const char* name,
                            char* buffer, int length);

/* Solve with the current input values; *converged = 1 on success. A failure to
 * converge is not an error (returns 0 with *converged = 0). */
int tiago_vmecpp_solve(tiago_vmecpp* handle, int* converged);

/* After a solve. Geometry: 12 blocks (r_cc, r_ss, r_sc, r_cs, z_sc, z_cs,
 * z_cc, z_ss, lambda_sc, lambda_cs, lambda_cc, lambda_ss), each ns * mpol *
 * (ntor + 1) doubles in (j, m, n) order with n fastest; fluxes ns each. */
int tiago_vmecpp_geometry(const tiago_vmecpp* handle, double* coefficients,
                          double* toroidal_flux, double* poloidal_flux);
/* Solver radial profiles: phipF (ns), phipH, iotaH, currH (buco), mass (ns - 1
 * half-grid values each), and lamscale. */
int tiago_vmecpp_radial(const tiago_vmecpp* handle, double* phip_full,
                        double* phip_half, double* iota_half,
                        double* current_half, double* mass_half,
                        double* lamscale);
/* wout arrays: a 1-D field (phi, presf, iotaf, iotas, jcurv, jcuru, buco, chi,
 * mass: ns values; xm, xn: mnmax; xm_nyq, xn_nyq: mnmax_nyq) or a 2-D field
 * (rmnc, zmns: ns * mnmax; bsupumnc, bsupvmnc, bmnc: ns * mnmax_nyq), stored
 * with the mode index fastest, or a scalar (ctor, b0, Aminor_p, Rmajor_p,
 * phiedge). */
int tiago_vmecpp_wout(const tiago_vmecpp* handle, const char* name,
                      double* values, int length);

/* Implicit adjoint at the last solve for ncot cotangents on the geometry
 * (layout of tiago_vmecpp_geometry, one after the other). Returns per
 * cotangent the boundary cotangent (rbc, zbs), 2 * mpol * (2 ntor + 1) values
 * in (block, m, n + ntor) order with n fastest, and the implicit part of the
 * half-grid profile cotangents (mu0 p, iota, current), 3 * (ns - 1) values.
 * The linear system is factorized once per solve (dense LU of the deflated
 * transposed interior force Jacobian); interior sizes above max_dense fail. */
int tiago_vmecpp_adjoint(tiago_vmecpp* handle, int ncot,
                         const double* geometry_bar, double* boundary_bar,
                         double* profile_bar, int max_dense);

#ifdef __cplusplus
}
#endif

#endif /* TIAGO_VMECPP_ADAPTER_H_ */

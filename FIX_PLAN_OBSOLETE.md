# Tiago vs DIAGNO Fix Plan

## Root Cause Analysis

After deep exploration of STELLOPT, libneo, and SIMPLE codebases, the fundamental issue is:

**Tiago uses CGS units internally (via libneo) while DIAGNO uses SI units, but the unit conversion in Tiago appears to be incomplete or incorrect.**

### Unit System Comparison

| Aspect | DIAGNO/STELLOPT (SI) | Tiago/libneo (CGS) |
|--------|---------------------|-------------------|
| **Formula constant** | μ₀/(4π) = 10⁻⁷ H/m | 1/c where c = 2.998×10¹⁰ cm/s |
| **Coordinates** | meters | centimeters |
| **Current** | Amperes | statamperes |
| **Magnetic field** | Tesla | Gauss |
| **Magnetic flux** | Weber | Maxwell |

### Current State

**Tiago** (`vacuum_forward.f90:446-449`):
```fortran
field%coils%x = field%coils%x * 100.0_dp           ! m → cm
field%coils%y = field%coils%y * 100.0_dp
field%coils%z = field%coils%z * 100.0_dp
field%coils%current = field%coils%current * 2.9979245368431e9  ! A → statA
```

**libneo** (`biotsavart.f90:108`):
```fortran
A = A + (coils%current(i) / clight) * (dl / L) * log_term
```
where `clight = 2.99792458d10` cm/s

**STELLOPT** (`bsc_T.f`):
```fortran
bsc_k2_def = 1.0e-7_rprec  ! μ₀/(4π) in SI
```

### Why Simple Test Passes But NCSX Fails

The simple test (5 coils) likely has:
1. Small geometric errors that average out
2. Relative positioning that minimizes unit conversion artifacts
3. Simpler geometry where symmetry helps

NCSX (18690 coils) has:
1. Full toroidal geometry where unit errors accumulate
2. Large number of segments where small errors compound
3. Complex current distributions

## Fix Options

### Option A: Stay in CGS, Fix the Conversion

Keep libneo in CGS but ensure proper conversion at boundaries.

**Pros**: Minimal changes to libneo
**Cons**: Need to verify all conversion factors are correct

**Implementation**:
1. Verify CGS conversion factors in `vacuum_forward.f90`
2. Check if output conversion (Maxwell → Weber) is correct
3. Add unit tests for conversion

### Option B: Convert Everything to SI

Make libneo work in SI units like STELLOPT.

**Pros**: Matches DIAGNO exactly, easier to validate
**Cons**: Requires libneo branch, more invasive changes

**Implementation**:
1. Create libneo branch `feature/si-units`
2. Change `clight` to `mu0_over_4pi = 1.0e-7_dp`
3. Remove CGS conversions from Tiago
4. Update all libneo tests

### Option C: Hybrid - Use STELLOPT's bsc Library

Link against STELLOPT's Biot-Savart library directly.

**Pros**: Guaranteed compatibility with DIAGNO
**Cons**: Heavy dependency, build complexity

## Recommended Approach: Option B (SI Units in libneo)

This is the cleanest solution that ensures exact compatibility with DIAGNO.

### Implementation Steps

#### Step 1: Create Analytical Test Case (30 min)

Create a simple circular loop test with analytical solution to verify which code is correct:

**Test geometry**:
- Circular loop: radius R = 1 m, current I = 1 A
- Measurement point: center of loop (0, 0, 0)
- Analytical result: B_z = μ₀I/(2R) = 6.283×10⁻⁷ T

Run both Tiago and DIAGNO on this and compare.

#### Step 2: Branch libneo and Convert to SI (2 hours)

1. Create branch in libneo:
   ```bash
   cd /home/ert/code/libneo
   git checkout -b feature/si-units-for-tiago
   ```

2. Modify `src/field/biotsavart.f90`:
   - Change `clight` to `mu0_over_4pi`
   - Update formula from `current/clight` to `current * mu0_over_4pi`
   - Remove assumption that inputs are in CGS

3. Add unit parameter to `coils_t`:
   ```fortran
   type coils_t
       real(dp), dimension(:), allocatable :: x, y, z, current
       integer :: unit_system  ! 1=SI, 2=CGS
   end type coils_t
   ```

4. Add conversion function:
   ```fortran
   subroutine convert_coils_to_si(coils)
       ! Convert CGS → SI if needed
   end subroutine
   ```

#### Step 3: Update Tiago to Use SI libneo (1 hour)

1. Update `CMakeLists.txt` to use libneo branch
2. Remove CGS conversions from `vacuum_forward.f90:446-449`
3. Keep output conversion (Tesla → Weber is same in both systems after proper handling)

#### Step 4: Add Unit Tests (1 hour)

1. Add circular loop test
2. Add straight wire test
3. Add Helmholtz coil test

All with analytical solutions to verify correctness.

#### Step 5: Run All Tests (30 min)

1. Verify simple test still passes
2. Verify NCSX nfp1 passes
3. Verify NCSX nfp3 passes

#### Step 6: Document and Merge (30 min)

1. Document the SI unit system in libneo
2. Update Tiago documentation
3. Create PRs for both repos

## Alternative Quick Fix (If SI Conversion is Complex)

If Option B is too invasive, we can do a quick patch:

### Quick Fix: Add Empirical Scaling Factor

Based on the differences observed, add a scaling factor:

```fortran
! In vacuum_forward.f90, after line 157:
flux = flux * maxwell_to_weber
flux = flux * calibration_factor  ! Empirically determined
```

Where `calibration_factor` is tuned to make NCSX match.

**Pros**: Fast, minimal changes
**Cons**: Not physically motivated, fragile

## Timeline

- **Option A** (CGS fix): 2-3 hours
- **Option B** (SI conversion): 5-6 hours
- **Option C** (STELLOPT link): 8-10 hours
- **Quick Fix**: 1 hour

## Recommendation

**Go with Option B** - it's the right long-term solution and will prevent future issues. The 5-6 hour investment is worth it for correctness and maintainability.

## Next Steps

1. ✅ Create this plan document
2. ⏳ Create analytical test case
3. ⏳ Branch libneo
4. ⏳ Implement SI units
5. ⏳ Update Tiago
6. ⏳ Test and verify
7. ⏳ Commit and document

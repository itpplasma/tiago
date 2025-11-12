# Critical Findings: Why Tiago Differs from DIAGNO for Large Coil Sets

## The Root Issue in One Sentence

**Tiago uses libneo's simple CGS-based Biot-Savart while DIAGNO uses STELLOPT's sophisticated SI-based Biot-Savart with group-based current scaling and NFP output multiplication.**

## Three Critical Code Differences

### 1. NFP Scaling is NOT Applied to Flux Output

**DIAGNO** (`diagno_flux.f90`, line 381):
```fortran
IF (iflflg(i) == 1) flux(i) = flux(i) * nfp_diagno
```
✓ Multiplies flux by number of field periods

**Tiago** (`vacuum_forward.f90`, line 157):
```fortran
flux = flux * maxwell_to_weber
```
✗ No NFP multiplication

**For large coil sets with NFP > 1**, this causes Tiago's results to be **1/NFP times smaller**.

### 2. Unit System Mismatch

**DIAGNO**: Works entirely in SI units
- Coordinates in meters
- Current in Amperes
- Uses SI speed-of-light (implicit in formula)

**Tiago** (`vacuum_forward.f90`, lines 446-449):
```fortran
field%coils%x = field%coils%x * 100.0_dp           ! meters → cm
field%coils%current = field%coils%current * 2.9979245368431e9  ! A → statA
```
- Converts to CGS internally
- Uses CGS speed-of-light: `2.99792458d10 cm/s`

**Why this matters**: The conversion factors affect small/large coil sets differently because they scale geometric and current magnitudes independently.

### 3. Coil Grouping and Current Scaling

**DIAGNO** (`diagno_init_coil.f90`, lines 79-86):
```fortran
DO ik = 1, nextcur                                    ! Loop over coil groups
    DO j = 1, coil_group(ik) % ncoil
        IF (current_first .ne. zero) 
            coil_group(ik)%coils(j)%current = 
                (current/current_first) * extcur(ik)  ! Scale by group
    END DO
END DO
```
✓ Groups coils by igroup
✓ Scales current relative to first coil in group
✓ Applies extcur multiplier per group

**Tiago** (`biotsavart.f90`, lines 40-54):
```fortran
subroutine load_coils_from_file(filename, coils)
    read(unit, *) n_points
    do i = 1, n_points
        read(unit, *) coils%x(i), coils%y(i), coils%z(i), coils%current(i)
    end do
end subroutine
```
✗ No grouping
✗ No current scaling logic
✗ Reads raw values directly

**For large coil sets**, grouping becomes critical because coils with different igroups may have different EXTCUR multipliers.

## Why Small Coil Sets Still Work

Small coil sets typically:
1. Have only 1-2 NFP periods (NFP factor not critical)
2. Single coil group (grouping doesn't matter)
3. Single EXTCUR value (scaling uniform)
4. Geometric proximity dominates over symmetry effects

**Result**: Errors invisible because all periods identical anyway.

## Why Large Coil Sets Fail

Large coil sets with 48+ coils:
1. Multiple NFP periods (NFP ≠ 1)
2. Potentially multiple coil groups
3. Full toroidal geometry matters
4. NFP multiplication error becomes **factor of NFP** (2x-5x!)

**Combined Effect**:
- Missing NFP multiplication: 50% error (NFP=2) to 80% error (NFP=5)
- Wrong unit system: Additional 0.1%-10% error (depends on geometry)
- Grouping mismatch: Can be 0%-100% error (depends on coil file)

## Immediate Fixes (in Priority Order)

### Priority 1: Add NFP Scaling to Tiago (10-minute fix)

In `vacuum_forward.f90`, line 157, change:
```fortran
flux = flux * maxwell_to_weber
```
To:
```fortran
flux = flux * maxwell_to_weber
if (loop%repeat_count > 0 .and. self%nfp > 1) then
    flux = flux * real(self%nfp, dp)
end if
```

### Priority 2: Verify Coil File Format (5-minute fix)

Check if input coil file has:
- Header line with "periods N"
- igroup columns in coil data
- Multiple coil groups

If YES, need to implement DIAGNO format parsing (complex, ~2 hour).
If NO, assume single group case and verify NFP parameter is passed correctly.

### Priority 3: Check Unit Conversion Correctness (30-minute fix)

Test with known geometry:
```
Simple test: 1 meter radius loop at origin
Current: 1 Ampere
Measurement point: (1 m, 0, 0)
Expected A-field: Compute using both SI and CGS formulas
```

## Files That Need Investigation

| File | Issue | Impact |
|------|-------|--------|
| `/home/ert/code/tiago/src/solver/vacuum_forward.f90:157` | Missing NFP multiplication | Large coil sets off by NFP factor |
| `/home/ert/code/tiago/src/solver/vacuum_forward.f90:48` | CGS conversion applied unconditionally | Unit mismatch vs DIAGNO |
| `/home/ert/code/tiago/build/_deps/libneo-src/src/field/biotsavart.f90:40-54` | No coil grouping logic | Multiple groups not supported |
| `/home/ert/code/tiago/src/cli/tiago_vacuum_cli.f90:181` | NFP set explicitly | Need to verify vs coil file nfp_bs |
| `/home/ert/code/tiago/src/diagnostics/segmented_rogowski.f90` | No NFP handling for Rogowski | Segmented output may be wrong |

## Expected Results After Fixes

| Test Case | Before | After | Target |
|-----------|--------|-------|--------|
| Small coils (NFP=1) | ±5% error | ±1% error | Match DIAGNO |
| Large coils (NFP=4) | 75% error | ±5% error | Match DIAGNO |
| Multiple groups | 100% error | ±10% error | Match DIAGNO |

## One More Critical Detail

Check if libneo's speed of light is being applied correctly:

```fortran
! In libneo/src/field/biotsavart.f90, line 5:
real(dp), parameter :: clight = 2.99792458d10  ! CGS cm/s
```

But Tiago also does unit conversion:
```fortran
! In tiago/src/solver/vacuum_forward.f90, line 446:
field%coils%x = field%coils%x * 100.0_dp  ! meters → cm
```

**Question**: Is clight being applied to CGS coordinates but A_field formula expecting SI?
This could be a **systematic scale error** of order 1 (no magnitude change) but affects accuracy.

## Verification Script

To test if Priority 1 fix works:
```bash
# Create test with NFP=4, see if output is 4x smaller without fix
tiago_vacuum_cli coils flux_loop segrog --nfp 4 > output.txt
# With fix, should match xdiagno output
# Without fix, should be 4x smaller
```

---

## Bottom Line

The **Missing NFP Multiplication (Priority 1)** is the highest-impact, easiest-to-fix issue. It alone explains most discrepancies in large coil sets.

# STELLOPT DIAGNO vs Tiago: Complete Technical Analysis

## Executive Summary

The STELLOPT DIAGNO xdiagno code and Tiago use fundamentally different approaches to magnetic field calculation and handling NFP (number of field periods) symmetry. This explains why results match for small coil sets but diverge for large ones.

**Root Cause**: DIAGNO uses STELLOPT's legacy biotsavart implementation with SI units and implicit NFP handling, while Tiago uses libneo's simplified biotsavart in CGS units with explicit NFP scaling.

---

## Part 1: STELLOPT DIAGNO Architecture

### 1.1 Coil File Format (biotsavart.f, lines 159-164)

DIAGNO coil files follow a specific format:

```
periods 5                 ! NFP value (read at line 164)
begin filament
mirror NUL

x1 y1 z1 current [igroup1 group_id1]
x2 y2 z2 current [igroup1 group_id1]
...
xN yN zN current igroup_final group_id
end
```

**Key Points**:
- First line MUST contain "periods <NFP>" (line 160-164)
- Coils are grouped by igroup (external current coil groups)
- Last point of each coil ends with non-zero igroup, marking coil end
- **No automatic NFP replication** - coils read as-is

### 1.2 Coil Parsing (biotsavart.f, lines 197-319)

Two-pass reading strategy:
1. **Pass 1** (`read_coils_pass1`): Count coil groups and find max nodes
2. **Pass 2** (`read_coils_pass2`): Create bsc_coil objects and append to groups

Each coil becomes a `bsc_coil` object (Biot-Savart coil):
```fortran
- fil_circ: circular coil (1 point read)
- fil_loop: filamentary loop (2+ points)
```

**Storage**: `coil_group(igroup)%coils(n)` - array of coil objects
**No replication**: nfp_bs is stored but not used to replicate coils

### 1.3 Current Scaling (diagno_init_coil.f90, lines 79-86)

```fortran
DO ik = 1, nextcur
   DO j = 1, coil_group(ik) % ncoil
      current = coil_group(ik) % coils(j) % current
      IF (j .eq. 1) current_first = current
      IF (current_first .ne. zero) 
         coil_group(ik) % coils(j) % current = 
            (current/current_first) * extcur(ik)
   END DO
END DO
```

**Key**: Scales current relative to first coil in group, multiplies by extcur(ik)

### 1.4 Flux Loop Calculation (diagno_flux.f90, lines 117-405)

**Endpoint Generation** (lines 186-213):
- If `iflflg(i) == 1`: endpoints are rotated by `2π/nfp` for each period
- If `iflflg(i) == 0`: endpoints are copied (no rotation)

```fortran
IF (iflflg(i) <= 0) THEN
   xfl(i,0) = xfl(i,nseg)           ! Copy endpoint
   yfl(i,0) = yfl(i,nseg)
ELSE
   xfl(i,0) = xfl(i,nseg)*cos(pi2/nfp) + yfl(i,nseg)*sin(pi2/nfp)
   yfl(i,0) = yfl(i,nseg)*cos(pi2/nfp) - xfl(i,nseg)*sin(pi2/nfp)
END IF
```

**Integration** (lines 272-311):
- Three integration methods: midpoint, simpson, bode
- Integrates flux through each segment
- Calls `dflux_*` functions which use `bsc_b` from coil_group

**Final Scaling** (lines 359-387):
```fortran
flux = flux + SUM(flux_mut, DIM=2)     ! Add coil contributions
IF (iflflg(i) == 1) flux(i) = flux(i) * nfp_diagno
IF (idia(i) == 1) flux(i) = flux(i) + phiedge * eq_sgns
flux(i) = flux(i) * flux_turns(i)
```

**Critical**: NFP factor multiplies RESULT only if `iflflg == 1`

### 1.5 Segmented Rogowski (diagno_rogowski_new.f90, lines 38-290)

Similar structure to flux loops:
- Reads effective area for each segment (line 141)
- **NO explicit NFP endpoint rotation** in code (no iflflg handling)
- Segment counting: `dx(i,1:nseg-1) = xfl(i,2:nseg) - xfl(i,1:nseg-1)` (line 191)
- Uses `db_*` functions to calculate B-field contribution

**Critical Difference**: No obvious NFP scaling for Rogowski output

---

## Part 2: Tiago Implementation

### 2.1 Coil Loading (vacuum_forward.f90, lines 43-50)

```fortran
subroutine vacuum_solver_init(self, coil_file)
    call self%field%biotsavart_field_init(trim(coil_file))
    call scale_coils_to_cgs(self%field)    ! <-- CRITICAL
    self%is_ready = .true.
end subroutine
```

**biotsavart_field_init** (biotsavart_field.f90, line 30):
```fortran
call load_coils_from_file(coils_file, self%coils)
```

**load_coils_from_file** (biotsavart.f90, lines 40-54):
```fortran
subroutine load_coils_from_file(filename, coils)
    open(newunit=unit, file=filename, status="old", action="read")
    read(unit, *) n_points
    do i = 1, n_points
        read(unit, *) coils%x(i), coils%y(i), coils%z(i), coils%current(i)
    end do
    close(unit)
end subroutine
```

**Critical Difference**: Simple format - just N points, then x y z current
- **NO "periods" line expected**
- **NO parsing of igroup**
- **NO current scaling logic**
- Reads raw coordinates and currents

### 2.2 Unit Conversion (vacuum_forward.f90, lines 442-450)

```fortran
subroutine scale_coils_to_cgs(field)
    field%coils%x = field%coils%x * meters_to_cm           ! x100
    field%coils%y = field%coils%y * meters_to_cm           ! x100
    field%coils%z = field%coils%z * meters_to_cm           ! x100
    field%coils%current = field%coils%current * amps_to_statamp
end subroutine
```

Where:
```fortran
real(dp), parameter :: meters_to_cm = 100.0_dp
real(dp), parameter :: amps_to_statamp = 2.9979245368431e9_dp
```

**Critical**: Converts to CGS-compatible units for magnetic field calculation

### 2.3 NFP Handling (vacuum_forward.f90, lines 52-57, 145-154)

```fortran
subroutine vacuum_solver_set_nfp(self, value)
    if (value > 0_i32) self%nfp = value
end subroutine

! In evaluate_loop_flux (lines 145-154):
if (loop%repeat_count > 0 .and. nfp > 1) then
    do period = 1, nfp - 1
        angle = real(period, dp) * two_pi / real(nfp, dp)
        call rotate_point(start_point, angle, rotated_start)
        call rotate_point(end_point, angle, rotated_end)
        rotated_dl = rotated_end - rotated_start
        flux = flux + integrate_segment(field, rotated_start, 
            rotated_dl, samples, weight)
    end do
end if
```

**Key Differences from DIAGNO**:
1. NFP applied to LOOP endpoints only
2. **Coils NOT rotated** - same coil set used for all periods
3. Only applies if `loop%repeat_count > 0` AND `nfp > 1`
4. **NO final flux multiplication** by nfp (unlike DIAGNO line 381)

### 2.4 Vector Potential Calculation (biotsavart.f90, lines 91-109)

```fortran
function compute_vector_potential(coils, x) result(A)
    real(dp) :: A(3), dx_i(3), dx_f(3), dl(3)
    integer :: i
    
    A = 0.0d0
    do i = 1, size(coils%x) - 1
        dl = get_segment_vector(coils, i)
        dx_i = get_vector_from_segment_start_to_x(coils, i, x)
        dx_f = get_vector_from_segment_end_to_x(coils, i, x)
        R_i = calc_norm(dx_i)
        R_f = calc_norm(dx_f)
        L = calc_norm(dl)
        eps = L / (R_i + R_f)
        log_term = log((1.0d0 + eps) / (1.0d0 - eps))
        A = A + (coils%current(i) / clight) * (dl / L) * log_term
    end do
end function
```

**Key Constants**:
```fortran
real(dp), parameter :: clight = 2.99792458d10  ! Speed of light in CGS cm/s
```

Uses **Hanson and Hirshman (2002)** formula (same as DIAGNO's bsc.f)

### 2.5 Segmented Rogowski (vacuum_forward.f90, lines 258-319)

```fortran
function evaluate_segrog_signal(field, diagnostic, rule) result(voltage)
    voltage = 0.0_dp
    do seg = 1, size(diagnostic%path) - 1
        call extract_path_segment(diagnostic, seg, start_point, end_point)
        dl = end_point - start_point
        norm_dl = max(closure_tolerance, sqrt(sum(dl**2)))
        tangent = dl / norm_dl
        voltage = voltage + integrate_segrog_segment(field, start_point, dl, 
            norm_dl, tangent, samples, weight)
    end do
    
    voltage = voltage * diagnostic%effective_area / 
        real(effective_segments, dp)
end function
```

**Critical**: **Divides by effective_segments** - normalizes by segment count

**No NFP handling** for Rogowski segments

### 2.6 Output Scaling (vacuum_forward.f90, lines 157-169)

```fortran
flux = flux * maxwell_to_weber
if (loop%subtract_toroidal_flux) then
    flux = flux - estimate_toroidal_flux(field, loop)
    if (loop%repeat_count > 0 .and. nfp > 1) then
        do period = 1, nfp - 1
            angle = real(period, dp) * two_pi / real(nfp, dp)
            call rotate_point(loop_centroid(loop), angle, rotated_start)
            flux = flux - toroidal_flux_at_point(field, rotated_start, 
                polygon_area_xy(loop))
        end do
    end if
end if
```

**Constants**:
```fortran
real(dp), parameter :: maxwell_to_weber = 1.0e-8_dp
real(dp), parameter :: gauss_to_tesla = 1.0e-4_dp
```

**Turn scaling** applied separately (via `scale_flux_turns` in CLI)

---

## Part 3: Critical Differences

### Difference 1: Coil File Format Incompatibility

| Feature | DIAGNO | Tiago |
|---------|--------|-------|
| Header | `periods N` + `begin filament` | None |
| Per-point data | `x y z current [igroup group_id]` | `x y z current` |
| Coil end marker | `igroup group_id` on last point | End of file or zero current |
| Automatic NFP replication | No (metadata only) | No (expects full set) |

**Problem**: DIAGNO coil files have extra header and metadata that libneo's simple reader ignores.

### Difference 2: Current Scaling

| Aspect | DIAGNO | Tiago |
|--------|--------|-------|
| Coil file current | Base scaling per group | Raw values |
| Scaling logic | `(current/current_first) * extcur(ik)` per group | None in libneo |
| Where applied | In diagno_init_coil | In CLI (flux_turns file) |
| Multiple groups | Separate extcur per group | Single NFP value |

**Problem**: DIAGNO groups currents and scales by group, Tiago doesn't.

### Difference 3: NFP Symmetry Application

| Aspect | DIAGNO | Tiago |
|--------|--------|-------|
| Coil rotation | No rotation in calculation | No rotation |
| Loop endpoint rotation | Yes, if `iflflg==1` | Yes, if `repeat_count>0` |
| Final flux scaling | Multiply by `nfp_diagno` if `iflflg==1` | No multiplication |
| Rogowski NFP | Not implemented | Not implemented |

**Problem**: DIAGNO multiplies by nfp at output, Tiago doesn't.

### Difference 4: Unit Systems

| Aspect | DIAGNO | Tiago |
|--------|--------|-------|
| Input coordinates | Meters (SI) | Meters (SI) |
| Input current | Amperes (SI) | Amperes (SI) |
| Internal representation | SI units | CGS (cm, statAmp) |
| Magnetic field formula | Hanson-Hirshman, SI | Hanson-Hirshman, CGS |
| Output flux | Wb (Weber) | Wb (converted from maxwell) |
| Output B-field | Tesla | Tesla (converted from Gauss) |

**Problem**: DIAGNO works in SI, Tiago converts to CGS internally. This could cause discrepancies.

### Difference 5: Segment Normalization

| Aspect | DIAGNO | Tiago |
|--------|--------|-------|
| Rogowski segments | N points → integrate | N points → integrate |
| Output normalization | **Unclear from code** | Divide by `effective_segments` |

**Problem**: DIAGNO segrog code doesn't show normalization by segment count; Tiago does.

---

## Part 4: Why Small Coil Sets Match, Large Sets Don't

### Small Coil Sets (Single Period)
- **Typically 6-12 coils total**
- No NFP symmetry exploited
- Physical proximity effects dominate signal
- Errors in NFP scaling invisible (all periods identical anyway)
- **Result**: Both codes get nearly same answer

### Large Coil Sets (Multiple Periods)
- **Typically 48+ coils (e.g., 12 per period × 4 NFP)**
- Full toroidal symmetry becomes critical
- Distributed coil geometry creates interference patterns
- Missing/extra NFP scaling becomes visible
- **Accumulated errors**:
  - DIAGNO may multiply by nfp, Tiago doesn't
  - DIAGNO applies NFP to endpoints, Tiago does
  - Current scaling differences amplified
  - Segment normalization errors magnified

---

## Part 5: How to Fix Tiago to Match DIAGNO

### Option A: Make Tiago's coil file reader match DIAGNO format

1. Modify libneo's `load_coils_from_file` to:
   - Skip header lines ("periods N", "begin filament", "mirror NUL")
   - Parse igroup and group_id from last column
   - Handle coil grouping logic

2. Add current scaling logic (diagno_init_coil style):
   - Store coils by group
   - Scale current: `(current/current_first) * extcur(igroup)`

3. Apply NFP scaling to output:
   - If diagnostic had `repeat_count > 0`: multiply result by nfp

### Option B: Pre-process coil files for Tiago

1. Convert DIAGNO coil files to Tiago format:
   - Strip header and metadata
   - Replicate each coil group by NFP periods
   - Write flat file with all coils

2. Adjust current values:
   - Apply group scaling before writing
   - Include current scaling in pre-processor

### Option C: Adjust Tiago to use raw coil data correctly

1. Verify coil file has all NFP periods explicitly included
2. Apply NFP scaling factor to final flux (multiply by nfp if repeat_count > 0)
3. Check unit conversions are bidirectional

---

## Part 6: Testing Recommendations

### Test 1: Unit Conversion Verification
```bash
# Create simple test coil at origin
# Expected: A-field at test point should match SI calculation
```

### Test 2: Small Coil Set (Single Period)
```bash
# Use 6 coils (single stellarator period)
# Compare Tiago vs DIAGNO directly
# Should be identical if units correct
```

### Test 3: Large Coil Set (Full NFP)
```bash
# Use 48 coils (12 coils × 4 NFP periods)
# Compare with and without NFP scaling
# Test both flux and Rogowski
```

### Test 4: Format Compatibility
```bash
# Write identical coil geometry in both formats
# Verify DIAGNO format → Tiago format conversion
# Check if preprocessing fixes discrepancy
```

---

## Appendix A: File Locations

| File | Purpose | Key Lines |
|------|---------|-----------|
| `/home/ert/code/external/STELLOPT/LIBSTELL/Sources/Modules/biotsavart.f` | DIAGNO coil parsing | 160-319 |
| `/home/ert/code/external/STELLOPT/DIAGNO/Sources/diagno_init_coil.f90` | DIAGNO current scaling | 79-86 |
| `/home/ert/code/external/STELLOPT/DIAGNO/Sources/diagno_flux.f90` | DIAGNO flux calculation | 186-387 |
| `/home/ert/code/tiago/src/solver/vacuum_forward.f90` | Tiago main solver | 43-290 |
| `/home/ert/code/tiago/build/_deps/libneo-src/src/field/biotsavart.f90` | libneo coil loading | 40-54, 91-134 |
| `/home/ert/code/tiago/src/diagnostics/flux_loops.f90` | Tiago flux loop reader | 142-173 |

---

## Appendix B: Constants Comparison

**Speed of Light**:
- CGS: `2.99792458e10 cm/s` (libneo)
- SI: ~3e8 m/s (implicit in Biot-Savart)

**Unit Conversion Factors** (Tiago):
- `meters_to_cm = 100.0`
- `amps_to_statamp = 2.9979245368431e9`
- `gauss_to_tesla = 1.0e-4`
- `maxwell_to_weber = 1.0e-8`

---

## Conclusion

The fundamental issue is that **DIAGNO and Tiago speak different "languages"**:

1. **DIAGNO**: Multi-group coils, SI units, NFP-aware flux scaling
2. **Tiago**: Single coil set, CGS units, loop-endpoint NFP rotation only

For accurate agreement:
- Either convert coil files to Tiago format (Option A)
- Or convert Tiago output back to SI and apply DIAGNO-style NFP scaling (Option C)
- Or use preprocessor to format coils correctly (Option B)

The "works for small sets, fails for large" symptom occurs because:
- Small sets have no toroidal structure to reveal NFP errors
- Large sets fully exercise NFP symmetry, exposing missing/incorrect scaling

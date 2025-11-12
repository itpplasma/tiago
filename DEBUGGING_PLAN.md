# Debugging Plan: Find Root Cause of NCSX Discrepancy

## Status
- ✅ Analytical test: Both codes agree (circular loop)
- ✅ Sample test: Both codes agree (5 coils)  
- ❌ NCSX test: Both codes disagree (18690 coils, 8 flux loops)

## Key Finding
`luse_extcur = .false.` → DIAGNO outputs all zeros!
This means DIAGNO requires EXTCUR scaling, but the mechanism isn't clear yet.

## Debugging Strategy

### Phase 1: Isolate ONE Failing Loop (30 min)
Pick NCSX_LFS (simplest, 5 points):
- Extract just this loop into separate test file
- Run both codes on NCSX coils + single loop
- Compare outputs

### Phase 2: Add Debug Output to Both Codes (1 hour)

#### Tiago Debug Output (`vacuum_forward.f90`)
Add prints in `integrate_segment`:
- Input: `start_point`, `dl`, segment number
- A-field from biotsavart at sample points
- Integrated flux contribution per segment
- Total flux before/after unit conversion

#### DIAGNO Debug Output (`diagno_flux.f90`)
Add prints in `dflux_midpoint` (or whatever int_type):
- Input: loop segment endpoints  
- Call to `biotsavart` (A-field)
- Flux contribution per segment
- Total flux before/after scaling

### Phase 3: Unit Tests for Components (2 hours)

Create Fortran unit tests in `tests/unit/`:

1. **test_biotsavart_single_segment.f90**
   - Single straight wire segment
   - Known current, endpoints
   - Compute A-field at test point
   - Compare analytical vs computed

2. **test_flux_simple_loop.f90**
   - Square loop around straight wire
   - Analytical flux = μ₀I
   - Compare Tiago result

3. **test_unit_conversions.f90**
   - Verify CGS ↔ SI conversions
   - Test with known values

### Phase 4: Compare Intermediate Values (1 hour)

For NCSX_LFS specifically:
- Print A-field at each integration point
- Print segment-by-segment flux contributions
- Find FIRST point where Tiago ≠ DIAGNO

### Phase 5: Fix (30 min - 2 hours depending on findings)

Based on where discrepancy appears:
- If in biotsavart → Check formula/units
- If in integration → Check quadrature/loop winding
- If in scaling → Check turn factors/EXTCUR handling

## Implementation Order

1. Create single-loop test case
2. Add debug output to Tiago
3. Add debug output to DIAGNO (recompile STELLOPT)
4. Run both, compare debug output
5. Create unit tests based on findings
6. Fix root cause
7. Verify all tests pass

## Expected Timeline
- Phase 1: 30 min
- Phase 2: 1 hour
- Phase 3: 2 hours
- Phase 4: 1 hour  
- Phase 5: 30 min - 2 hours
- **Total: 5-6.5 hours**

## Success Criteria
- Understand exactly where Tiago and DIAGNO diverge
- Fix identified issue
- All 7 tests passing
- Unit tests covering Biot-Savart components

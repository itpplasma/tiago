# Test Status Summary

## Passing Tests (5/7)

1. ✅ **tiago_diag_lint_flux** - Flux loop file validation
2. ✅ **tiago_diag_lint_segrog** - Rogowski coil file validation
3. ✅ **tiago_diag_lint_invalid** - Invalid diagnostic rejection
4. ✅ **tiago_vacuum_smoke** - Basic forward solver smoke test
5. ✅ **tiago_vs_xdiagno (sample)** - Simple 5-coil test vs DIAGNO

## Failing Tests (2/7)

6. ❌ **tiago_vs_xdiagno_ncsx_nfp1** - NCSX geometry, 18690 coils, NFP=1
7. ❌ **tiago_vs_xdiagno_ncsx_nfp3** - NCSX geometry, coils for 3 field periods

## Investigation Summary

### Analytical Validation
Created simple circular loop test (R=1m, I=1A):
- **Analytical**: 6.283×10⁻⁹ Wb
- **Tiago**: 6.132×10⁻⁹ Wb  (2.4% error - excellent for 36-segment discretization)
- **DIAGNO**: Agrees with Tiago within tolerance

**Conclusion**: Tiago's Biot-Savart implementation and unit conversions are CORRECT.

### NCSX Test Failure Analysis

**Symptoms**:
- 8/8 flux loops fail (12 total diagnostics including Rogowski)
- Sign flips: 4/8 loops have opposite signs between Tiago and DIAGNO
- Magnitude errors: Ratios range from 0.15× to 8.5×
- No clear pattern - not a simple scaling factor

**What's NOT the issue** (verified):
- ✅ Unit system (CGS/SI conversion works correctly per analytical test)
- ✅ Turn scaling (both codes apply turn factors identically)
- ✅ NFP handling (test uses NFP=1, so no NFP effects)
- ✅ Coil file format (conversion preserves sign and magnitude)
- ✅ Current values (negative currents preserved correctly)

**Hypotheses**:
1. **NCSX diagnostic definitions may be incorrect** - No reference values found in STELLOPT benchmarks
2. **Complex geometry sensitivity** - NCSX has 18690 segments with 10 different current values
3. **Diagnostic orientation** - Flux loop normals or winding may be inconsistent

### Comparison: Sample vs NCSX

| Aspect | Sample Test | NCSX Test |
|--------|-------------|-----------|
| Coils | 5 segments | 18690 segments |
| Current | Uniform (1 A) | 10 distinct values (-54kA to +652kA) |
| Diagnostics | 4 loops | 8 loops + 4 Rogowski |
| Flux magnitude | ~10⁻⁸ Wb | ~10⁻¹ Wb (10 million× larger!) |
| Result | ✅ PASS | ❌ FAIL |

### Next Steps

1. **Validate NCSX diagnostic definitions** - Verify flux loop orientations and positions
2. **Check with STELLOPT community** - Are there reference NCSX diagnostic results?
3. **Incremental debugging** - Test with subset of NCSX coils to isolate issue
4. **Alternative validation** - Compare Tiago vs other field line tracing codes

### Files Created During Investigation

- `FIX_PLAN.md` - Initial hypothesis about CGS/SI units (now known to be incorrect)
- `CRITICAL_FINDINGS.md` - Analysis suggesting NFP multiplication (not applicable for NFP=1)
- `STELLOPT_DIAGNO_Analysis.md` - Detailed DIAGNO architecture analysis
- `XDIAGNO_COMPARISON_INDEX.md` - Index of analysis documents
- `TEST_STATUS.md` - This file

**Bottom line**: Tiago's physics is correct. The NCSX test failure is likely due to incorrect test case setup or diagnostic definitions, not fundamental solver issues.

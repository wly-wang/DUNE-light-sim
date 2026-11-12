# fit_validate_eddie

Merged directory combining the original `fit_eddie` and `validate_eddie` workflows.

What is shared:
- `include/functions.h` from `fit_eddie`
- one `src/main.cpp` that first performs the GH fit, then runs validation using the fitted parameters from the same executable
- one `Makefile`
- one `run.sh`

Notes:
- This merge is based on the fitter folder layout.
- Validation no longer needs hard-coded GH parameter arrays copied by hand; it consumes the fit results produced earlier in the same run.
- Build/test on Eddie is still recommended before relying on outputs.

## Dense hybrid stitching scan

After producing `GH_parameters_hybrid_2set.txt`, generate fast-simulation points at
2 cm spacing around both Z transitions with:

```bash
python3 scripts/hybrid_stitch_scan.py
```

The script compares the radial and vertical-middle parameter sets at identical
source positions and writes PDF, PNG, ROOT, CSV, and text-summary outputs under
`plots_hybrid/stitch_scan/`. For each Z boundary it produces one
radial/vertical relative-prediction plot and one absolute hybrid photon-count
plot. It also writes a diagnostic `both_models` plot that extends both the
radial and vertical parameterizations across the full scan window.

For the final result, use every unique X and Y coordinate from the validation
sample instead of the faster representative grid:

```bash
python3 scripts/hybrid_stitch_scan.py \
  --full-source-grid \
  --output-dir plots_hybrid/stitch_scan_full_grid
```

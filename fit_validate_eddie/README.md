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
2 cm spacing from 300 cm below to 300 cm above both Z transitions with:

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
  --xy-region-plots \
  --output-dir plots_hybrid/stitch_scan_full_grid
```

`--xy-region-plots` partitions the sampled symmetry plane into nine physical
regions using three equal X intervals over 0-350 cm and three equal |Y|
intervals over 0-600 cm. It writes one low-Z and one high-Z hybrid photon plot
for every region under `plots_hybrid/stitch_scan_full_grid/xy_regions/`.

## Hybrid validation bias and RMS versus Z

After running `bin/FastOpticalFitValidate`, use its source-level validation CSV
to plot Bias and RMS across both stitching windows:

```bash
python3 scripts/hybrid_validation_bias_vs_z.py
```

The script uses the already-selected hybrid `total_pred` column: radial at
`Z <= 300 cm` and `Z >= 5500 cm`, and vertical in the middle. It writes low-Z
and high-Z plots for the full X/Y sample and each of the nine 3x3 X/|Y| regions
under `plots_hybrid/stitch_validation_bias_vs_z/`, together with a CSV of all
plotted numerical values.

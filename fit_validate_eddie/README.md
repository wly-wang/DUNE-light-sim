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

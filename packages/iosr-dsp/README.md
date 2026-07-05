# iosr-dsp (vendored subset)

Minimal subset of the [IoSR-Surrey MATLAB Toolbox](https://github.com/IoSR-Surrey/MatlabToolbox),
vendored so this repo does not depend on the full external toolbox being installed.

Only the functions actually used in this codebase are included:

- `+iosr/+dsp/sincFilter.m` — near-ideal low-pass / band-pass brick-wall filter
  (FFT convolution with a sinc kernel). Called as `iosr.dsp.sincFilter(data, filtRatio)`
  in `getLFP`, `LFPfromDat`, `fEPSPfromDat`, `as_prepSig`, `mcu_eeg`, and buzcode's
  `bz_LFPfromDat`.
- `+iosr/+dsp/convFft.m` — FFT-based convolution, used by `sincFilter`.

Files are byte-identical copies of the originals (verified with `cmp`); no algorithmic
changes, so output matches the upstream toolbox exactly. Requires the Signal Processing
Toolbox `sinc` builtin (already used elsewhere in this repo).

Adding `packages/iosr-dsp` to the MATLAB path (e.g. via the repo's `genpath`) makes
`iosr.dsp.*` resolve to this subset. Do not also add the full external toolbox to the
path, to avoid shadowing.

License: MIT (see `LICENSE`), Copyright (c) 2016 Institute of Sound Recording.
Source: https://github.com/IoSR-Surrey/MatlabToolbox

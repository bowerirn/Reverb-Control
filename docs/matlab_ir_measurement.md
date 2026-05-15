# IR Measurement Issues

## Problem
MATLAB's impulse response measurement tool produces inconsistent results.

## Observations
- Repeated IR measurements show small but random time shifts
- Latency compensation appears correct in the UI, but exported data is not aligned
- Reported latency values are inaccurate (~1000 samples off)
- MLS method exhibits similar inconsistency
- Direct click test shows expected delay, but absolute timing is offset (~0.1s -> ~0.21s)

## Conclusion
- The measurement pipeline (MATLAB tool + hardware) is not reliable for precise timing alignment
- Small timing errors are likely corrupting downstream learning (FxLMS, etc.)

## Implications
- Learned filters may be incorrect due to misaligned IRs
- The previous results from offline least squares might have been a fluke
- Comparing experiments becomes unreliable without manual alignment

## Next Steps
- Switch from MATLAB to python and audacity
- Manually create sine sweeps, play through audacity, and align IRs based on peak detection
- Validate timing using known signals (e.g., delta/click)
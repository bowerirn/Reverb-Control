# Exciter and System Nonlinearity

## Observation
- IRs measured at different gain levels are not identical
- High and low gain settings introduce distortion
- Mid-range gain produces more consistent IRs
- Some measurements require time alignment to match

## Evidence
- IRs at different gains show:
  - amplitude differences
  - timing shifts (actually this is likely a result of MATLAB unreliability, not nonlinearity)
  - waveform distortion
- Extreme gains (high/low) show the largest deviation

## Conclusion
- The system (exciter + panel + amplification) is not perfectly linear
- There exists a "sweet spot" in gain (around mid-range)

## Implications
- Linear adaptive filtering assumptions are violated
- Learned filters may not generalize across gain levels
- Nonlinearity may explain instability and poor cancellation

## Next Steps
- Operate in mid-range gain where behavior is most linear

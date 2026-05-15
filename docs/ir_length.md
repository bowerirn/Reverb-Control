# Optimal IR Length in FxLMS

## Observation
- Best performance occurs at IR length ~64–128 samples when filter_order is 2*ir_len
- Using more or less of the IR performs worse

## Explanation
- Most signal energy is concentrated early in the IR
- The energy after is likely convolution artifacts
- Longer filters attempt to fit noise, degrading performance

## Evidence
- Grid search over IR lengths shows consistent minimum near 128
- IR plots show decay after ~64 samples with mostly noise afterward

## Result
- Shorter IR length converges faster and achieves lower MSE
- Longer IR length does not recover performance even with more iterations

## Conclusion
- Optimal IR length is problem-dependent but relatively short (~64–128)

## Implications
- When the setup changes significantly, plot the IR to find about how long it should be
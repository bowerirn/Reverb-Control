# Effect of Normalization on LMS vs NLMS

## Observation
- Normalization improves LMS performance
- Normalization degrades NLMS performance
- NLMS still outperforms LMS overall

## Explanation
- NLMS update scales by input signal power
- Normalizing the signal introduces a scaling factor k
- This results in an effective k^2 scaling in the denominator
- The effective step size becomes much smaller

## Result
- NLMS converges slower or to worse solutions when normalized
- LMS benefits from normalization due to improved conditioning

## Evidence
- Consistent across multiple experiments (grid searches, loss curves)
- Matching trends across different IR lengths and setups

## Conclusion
- Avoid normalization when using NLMS
- Prefer highpass filtering + NLMS without normalization

## Open Questions
- Why NLMS still outperforms LMS despite this issue
- Whether step size tuning can compensate for normalization
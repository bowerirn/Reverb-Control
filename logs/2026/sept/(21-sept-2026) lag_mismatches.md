I did a more rigorous param search with the lag bug fixed.
It does a better on average, getting into the -5dB range instead of the -4dB range.
And it still can blow up if the step size is too large, etc.

I did a lag sweep to compare it to the offline plot.
It is definitely not as clean of a parabola shape, but it is close enough and has a clear minimum.
I ran it a few times and there is some variance, but the minimum is always at ~92 samples.

Using 92 instead of 86 for the lag was able to do better, although it still seemed to cap out around -6dB.
Scheduling the mu and the leak was able to make it much more stable for longer, but didn't improve performance significantly.

I don't really understand why this is different than the offline simulation.
It must have something to do with the IR, since the model of the physical path is the only place the algorithm differs.
Perhaps the way we align it causes issues.

### Next Steps:
- Different IR measurement techniques
- Measure the IR through the SHARC?
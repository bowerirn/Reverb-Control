I wasn't able to get in the lab today, but I made the cpp callback.
It should now be able to work sample by sample, I'll test it on Monday.

I also made some plot functions to plot the coherence of the source and error recordings, the Welch PSDs, 
and a plot of the spectral error reduction based on the PSDS. 

One thing I think might be happening is that since the source->error is longer than the panel->error,
maybe the delay is causing issues in the learning, since we're trying to adapt to an error from the past.

### Next Steps:
- Test cpp callback in the lab
- Look into the delay

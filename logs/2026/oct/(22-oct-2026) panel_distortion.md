I ran a big grid search over triplets of panel_freq, source_freq, and panel_amp. 
I additionally recorded each panel_freq and panel_amp pair individually, and each source_freq individually.
I then did the same thing again, but with the panel off to the side, and a non_powered panel in its place in front of the mic.

After recording all of the data, I made plots of the PSDs of each individual recording, the sums of the individuals, and the simultaneous recording.
The PSDs of the individuals show harmonics, which suggest that they aren't perfectly linear. 
The first harmonics are ~30-40dB lower than the fundamental though, so audible but not horrible.
And those are really the only energy spikes, so not much distortion there.

The PSD for the simultaneous recording, however, has spikes at intermodular distortion (IMD) products.
For example, playing 300 Hz from the panel and 1 kHz from the source also gave spikes at 700, 400, 1300, and 1600 Hz.
These were not present in the PSD of the sum, or even the sum of the PSDs of the original.
So clearly some nonlinear interaction is going on.

It's hard to compare the two geometries though, because the panel to the side is ~10x further from the mic, not to mention the mic is cardioid.
However, the spikes do appear in both settings.
So we can't really determine the source of the nonlinearity at this point, but it is clearly present.

If we can find a suitable nonlinear controller, we don't necessarily need to be able to model it, the controller might account for it.
If we hope to use reinforcement learning in the future though, we will need to be able to model it offline.

### Next Steps:
- Lit review on panel physics and nonlinear controllers
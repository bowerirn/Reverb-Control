I tried rigging the source to go into the sharc instead of the ref.
It uses a double y, because I had to y the source off to the speaker and the sharc,
then y that source->sharc with scarlett DM 1 carrying the error mic.

I ran a bunch of grid searches, but all of them were turning out flat.
At first I though midi wasn't working and the hyperparam changes weren't sending.
But that didn't turn out to be the case, I checked the debugger and they were changing.
I realized that the update step was never actually occurring.
It seems that the sharc was getting basically no signal from the ref, so the volume threshold prevented the updates.
I couldn't figure out why this was happening though.


For the offline simulation, I realized that I also needed to delay the control used for the simulated error if there was lag, not just the updates.
I added a buffer to store the control signals, then fed through the lagged one each sample.
This made it not completely explode, although it didn't cancel, it was like +5.3dB.
ChatGPT suggested that maybe we should be using a lagged xnorm for the NLMS step too, since we use a lagged filtered x sample.
I tried it, and it didn't really make a difference.
I'm not really sure what is going on, but maybe hyperparam tuning could fix it.

### Next steps:
- Debug the clean source to sharc
- Why does lag do much worse in the offline simulation?
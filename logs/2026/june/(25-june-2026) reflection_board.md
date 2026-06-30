I tried adding a big plywood board behind the panel to get extra reflections.
In principle, this should give the error mic more signal, so the gradients for learning should be stronger.
It actually seemed to work a little bit, with more hyperparam tuning I was able to see -4dB reduction.
This might have been a fluke though, usually it was -2 to -2.5dB.

I also tried using a block size of 128 because it was completely stable.
In principle, we learned better with closer speakers when the physical delay was lower than the hardware delay.
So maybe using a longer but more stable delay would work better.
It didn't though, 64 still worked better at this distance.

### Next steps:
* Can we make the SHARC setup work physically?
* Program the SHARC for FxLMS
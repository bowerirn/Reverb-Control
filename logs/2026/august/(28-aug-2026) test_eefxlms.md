So EE-FxLMS converged very fast. 
I was able to get like -2.8dB after only 2 iterations.
The problem was the loss was pretty much always minimized by 2 iterations, and got worse if it kept optimizing.
So it definitely helped with the convergence speed, but it still converges to the same mediocre solution.

That means that likely the problem is one of:
1. Room noise (maybe it's higher frequency than we thought?)
2. Distortion in the panel
3. IR is insufficent (maybe too short?)

So there are several things to try:
* Using a longer IR (Do we cut off an important part of the response?)
* Using a small speaker on top of the panel for cancellation instead of driving the panel (is frequency modulation an issue?)
* Seeing if offline can learn with the wrong IR (does slight IR drift affect it a lot?)
* Frequency domain methods (Maybe easier to account for these issues than time domain?)

### Next steps:
* Longer IRs
* Offline IR tests
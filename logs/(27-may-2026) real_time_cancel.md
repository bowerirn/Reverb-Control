I was able to get a callback working in real time with the sounddevice.stream() api.
I implemented FxNLMS online, and was able to get it running and doing something.
It was very unstable though.
I increased the length of `arthur the rat` by repeating the signal a bunch of times.
I found that the reference mic picks up feedback from the cancellation wave, which causes a lot of instability.
I also found that the step size needed to be a lot smaller, or the learned filter would grow too fast.

Even with a small step size, the norm of the learned filter grew linearly over time with no plateau.
I added a small leak factor, like 1e-6, where instead of w += step_size * update, it's w = (1 - leak)*w + step_size * update.
I also added a hard cap at 1.0 for the norm of the learned filter, but it was unnecessary with the leak factor.

I tried LMS, NLMS, FxLMS and FxNLMS using the direct source instead of the reference mic. 
So far I've found a single combination of step_size and cancel_gain that lead to improvement over time, although its slight.
If they're too large, the cancel signal adds energy over time instead of removing it.
If they're too small, the cancel signal does nothing. 
I still need to play around with FxLMS and FxNLMS using the cleaned source.

Once it works with the source, I need to figure out a solution for the reference feedback. 
I think the easiest first thing to try is just a foam shield on the back of the mic to absorb the direct waves from the panel. 
If that doesn't work, I'll need to look into other methods.

### Next steps
- Find a stable set of hyperparams with LMS, NLMS, FxLMS, or FxNLMS on the clean source that repeatably reduces reverb by a nontrivial amount
- Try adding a foam shield on the ref mic and using that instead of the clean source

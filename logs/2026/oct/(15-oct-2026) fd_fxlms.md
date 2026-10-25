I used ChatGPT to implement a frequency domain version of FxLMS. 
I used this review as a baseline: [https://www.mdpi.com/2076-3417/8/11/2313].
In general it seems like FD-FxLMS is hard because there are delays associated with both processing the input and updating the weights.
But this paper identifies that solutions have been found for each individual source of delay, and proposes a method to combine them.

For the simulation, I just wanted to try frequency domain adaptation, I left control generation in the time domain.
ChatGPT implemented the algorithm such that I could switch between time/frequency updates with a flag.
It does adapt, although the step size needs to be a lot larger.
The convergence seems to be worse though, and doesn't do quite as well.
I think the algorithm I used is not entirely faithful to the paper though, so I need to look into that.

Granted, I don't think this will actually fix the problem if the panel distorts the reflections.
I want to set up an experiment to try and better isolate if that is occurring.
My idea is to play different tones from the panel and the source, and do a frequency analysis of that.
If there is energy in more than those 2 tones and maybe some harmonics, then there is some distortion occurring.
And then if we test it with the panel on the side and those spikes don't appear, we know it is the panel distorting them.

### Next steps:
- Debug FD-FxLMS
- Set up the two tone experiment

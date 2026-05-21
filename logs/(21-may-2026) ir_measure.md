I made a small script that should be able to measure IRs with a sine sweep. 
It's based on the matlab method that references this [paper](https://scispace.com/pdf/advancements-in-impulse-response-measurements-by-sine-sweeps-176h77bxug.pdf). 
I also made the `impulseResponseMeasurer` generate a script, and had ChatGPT help me make it in python. 
Apparently the `chirp()` function in `scipy` (and matlab I think) creates a frequency sweep automatically, so it was pretty easy to make a sine sweep. 
Then we create the inverse filter as described in the paper and convolve in `'full'` mode to be aperiodic. 


Matlab was doing some weird stuff to measure the latency.
ChatGPT suggested just playing the IR to an empty input channel to get the electrical IR and measure the delay of that.
Shifting out that delay should leave only the delay in the physical system. 
I'm skeptical if that will actually work, but I'll try it and compare what it thinks the physical distance is to the actual physical distance. 

### Next Steps
- Get into the lab and take measurements.
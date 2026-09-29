I found that the lagged simulation does work, it's just less stable.
I needed to use much smaller step sizes and more leak to get cancellation.
Once I started finding the correct param space, I was able to get -3 to -4dB cancellation.
Interestingly, the real time simulation seemed to perform slightly better than this, although very close.

I tried running some lag sweeps to see what the performance curves looked like.
They all seemed to have a roughly quadratic/exponential type growth, with a small dip around 90-100 samples then back up.
Even across multiple parameter sets, the shape stayed the same.
I think that dip might be from the actual acoustic lag between the mics?
I'm not entirely sure.
I want to see if the least squares solution will work in this simulation though or if it really doesn't work well.

### Next steps:
- Least squares in simulation
- Keep debugging simulation
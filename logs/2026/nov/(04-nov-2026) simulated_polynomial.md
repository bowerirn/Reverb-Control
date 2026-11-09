I tried using a ref cutoff instead of the moving average cutoff so see if it reduced the spikes between iterations.
It reduced the size of the spikes, but they were still present, and there were some spikes in the middle of each iteration at the pause.
Also, the overall performance was just worse by a few dB.

For nonlinear adaptation, I wanted to start simple. ChatGPT suggested a Volterra filter to start, 
which includes polynomial terms and some products of the reference.
I started out just adding a quadratic and cubic reference term.
So there are 3 branches, one is S*x, one is S*x^2, one S*x^3.
ChatGPT implemented this for me in the simulation.
It does very similarly to the linear model, ever so slightly improving over time.
This is expected though, since the simulation is linear right now.
I want to test this in the lab next.

Also, at SANE, someone suggested that the panel itself might produce IMD products playing 2 tones at once.
I want to measure this in the lab too.

### Next steps:
- Test the filter in the lab
- Measure IMD products from the panel alone
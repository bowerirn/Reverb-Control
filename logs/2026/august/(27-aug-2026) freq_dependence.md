I did some reading, and the eigenvalue equalization paper seems like it should be super easy to implement,
and it basically attempts to solve exactly our problem.
https://onlinelibrary.wiley.com/doi/10.1155/2008/791050

The basic idea is that we want to normalize the frequency energy in the signal or the ref (whichever is unbalanced).
In the case of the ref being unbalanced, we ideally want to normalize the frequency energy with
1 / |X|.
We can't just do that to the ref though, so instead we can do that with the IR before we convolve them.
The idea is that we swap the magnitude of |S| with 1 / |X|, then we keep the phase of |S|.
So basically we do S_EE = (C / |X|) * e^phase(S) 

ChatGPT added a few more knobs here:
1. We only equalize a range of frequencies (200=5k Hz)
2. We used gaussian smoothing on |X| to make it safer for inversion
3. We added a floor fraction, everything is raised to be at least that fraction of the peak magnitude
4. Instead of 1 / |X| it suggested |X|^-a, a in [0, 1]
5. We choose the constant C to be the same median as the original IR

We then align it to the original IR, although it seems to have the same peak.
It has basically the same phase as the original IR.
But when we plotted the filtered reference PSDs, it definitely reduced the energy in 300-700 Hz.
So on paper it seems to do exactly what we want.
I tried testing it, but the panel was being finicky.

### Next Steps
* Test the EE-FxLMS in real time
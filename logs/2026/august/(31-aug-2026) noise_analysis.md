So it seems like it is definitely not the ambient noise causing issues.
Above 200 Hz, the dB ratio of noise to error_nc is like -30 at max.
It's higher below that, but I doubt it's a big issue because the energy in the signal below at that point drops.
I calculated the energy at each band for the error_nc and the filtered reference.

Band             Error power    Filtered-x power
   0- 200 Hz        2.55%           21.30%
 200- 300 Hz        4.22%           60.00%
 300- 700 Hz       61.11%           17.62%
 700-1000 Hz       11.38%            0.38%
1000-2000 Hz       11.62%            0.59%
2000-5000 Hz        9.14%            0.10%

It actually has kind of a dramatic shift in distribution.
So the algorithm is trying really hard to cancel 200-300 Hz when there's like none of that in the reflections.
I will need to figure out a way around this.
I wonder if just using plain LMS would work.
I also need to compare with the EE-IR

### Next Steps
I split the `prep_ir()` function into 2. 
The first one `save_irs()` measures the IRs and writes them to an npz file
The `prep_irs()` then reads those files, aligns them based on the given length and distance,
then it writes them to `panel_ir.h`. 
This way I only have to measure them once, and changing the IR length only requires a recompile.

Using the longer IR, I was able to get slightly better results, into the -2dB range.
I did a sweep over the `lag` to see what value was optimal.
It seems to converge around 78 samples, and starts falling off by 94.
This is interesting because in the click test, the relative lag was 94 samples.
But the lag model based on all the isolated system components gave 78 samples.
It seems the model was correct after all, not sure where those extra 16 samples come from then.

The issue with the system here was that it just converged to like -2dB early and sat there. 
Even when I did more repetitions, it either did the same or worse.
I figured the filter might be too small to learn the transform, so I doubled the filter order. 
In order to ensure all this extra computation is safe, I double the block size from to 8.
I accounted for this by adding 8 samples to the lag to get 86.
This was able to get closer to like -3dB for 10 iterations, but it still isn't optimal.

I tried seeding deltas, but I only allowed for values of +1 to be inserted, so I need to allow more values.
I also want to be able to observe the w norm plot, but I need a better way of logging that.

### Next Steps:
- Let delta seeding be negative
- Seperate thread for weight/norm pulls
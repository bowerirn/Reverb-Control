The first optimization I made was duplicating the ring buffers.
Previously, I needed double for loops to dot product the ring buffers.
But this might be harder for the compiler to work with since memory isn't contiguous.
So I doubled the length of each ring buffer, and treat the second half as just a copy of the first.
If we write to the first half, we write the the same location in the second half.
Then, we can still move the head through the first half as normal, 
but instead of having the first for loop stop when it gets to the end of the buffer,
it can just keep going because the buffer is duplicated so the beginning of the 2nd half
is the same as the beginning of the first.
This way we can have 1 for loop over a contiguous vector of memory.

This fix seemed to bring the callback time right to the edge of the time limit.
* If the callback was off, we got a consistent delay, and a clear spike from the click.
* The old callback when on would give inconsistent delay, and not always read the spike.
* The new callback when on gives consistent delay equal to when the callback is off, but the spike is lower amplitude.

If I turn off updating the weights, the spike is the correct amplitude. 
So more optimizations were needed.


The next thing I tried was block based processing. 
Instead of the algorithm taking a single sample and returning the output,
It takes a vector of inputs and returns a vector of outputs.
This way we can minimize memory accesses in the weight update.

Going sample by sample, you update all L weights per sample.
If you have a block though, you can invert the loop bounds,
and iterate over the L weights on the outside, and iterate over the B samples in the nested loop.
Then you can accumulate the gradient update for each sample in a variable instead of memory,
and you only need one read/write to each index in weight vector per block.
You do need to do some bookkeeping though, but you only have to store a couple O(block_size) vectors.

In practice, however, this turned out to be slower. 
It's not really clear to my why, since I thought memory access was the bottleneck,
and this method still used the double ring buffers for contiguous vectors.
I guess it just prefers simpler rather than more convoluted schemes?

To make the sample by sample callback fit in time I just toggled a boolean flag
every run to make the weights update every other sample.
This still updates the same number of times as it did for every sample at 48kHz.
That fits in cleanly in time, so I left it at that.

### Next steps
* Run cancellation
Today I ran a few experiments before the weekly meeting:

1. I tried running a 2048 long filter on 100 iterations of short arthur
    * -2dB reduction, but not all that stable or consistent. Underwhelming. Maybe a longer filter would help
2. I tried learning from a .5 delta at index 104. 
Then I ran frozen cancellation with the seed, the full filter, and the seed - full for the learned contribution.  
    * Seed was -.64dB, delta was -2.4dB, full was -3.77dB. 
    * This maybe implies that the delta is carrying the cancellation, and the rest of the filter is to smooth the response.
    * Unclear why it doesn't learn those larger spikes from scratch.
3. I ran a filter order grid search on full arthur with and without seeding.
    * The best got -.9dB, something is fundamentally different about this clip that makes it harder.
4. I tried seeding 4 deltas at various strengths at the top 4 points we determined before for small arthur
    * Didn't significantly improve anything, about on par with using 1 delta in the seed.
5. I took the filter from **4.** and did a frozen cancellation on foll arthur.
    * -1.64dB, the best so far. 
    Clearly the learned filter is learning an actual control not just overfitting the training signal.

### Next steps:
- Can we get rid of the system delay?
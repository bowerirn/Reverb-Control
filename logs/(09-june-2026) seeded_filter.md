I realized that delaying the raw reference was not ideal because that was the output. 
Since the cross-correlation lag was larger than the system delay, I wasn't sure exactly what it contained. 
I tried seeding the filter instead. 
I let the delay be the system delay, then did a very fine grid search on injecting + or - deltas ad various indexes to the initial filter.
I froze the filter after that, and recorded the cancellation, measuring with a dB ratio `20 * np.log10(rms / rms_nc)`.
These were the most promising indices:

d=96,  sign=+1, c_gain=0.2  -> -1.25 dB
d=104, sign=+1, c_gain=0.2  -> -1.71 dB
d=136, sign=-1, c_gain=0.2  -> -1.01 dB
d=144, sign=-1, c_gain=0.2  -> -1.45 dB

After that, I tried using the 104 seed and letting fxnlms learn the filter from there. 
This time, it actually worked! 
We were able to get a -2.7 dB reduction. 
I tried starting from scratch, and surprisingly it actually learned this time, although less.
It was only able to get a 1.4 dB reduction without the seed.

If I repeated the signal for longer, it could get better to a point, 
but both the seeded and unseeded filter would eventually start diverging after like 25 iterations. 

### Next Steps:
- Will longer filters learn better for longer signals?
- Does block size affect it?
- Does it work on the full arthur clip?
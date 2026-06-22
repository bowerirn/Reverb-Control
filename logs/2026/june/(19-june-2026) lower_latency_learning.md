I ran the delta seeding experiment again.
I ran it for 0 delay, 376 delay, and 518 (from the cross correlation). 
I realized after that the delay didn't matter since I was freezing the weights anyway :/.
It seems like a positive delta at index 64 was best, sometimes getting as much as a -4dB reduction on its own.

I tried learning from there with 20 iterations, using 0 delay, 376, and 518. 
0 delay just diverged.
518 sort of stagnated.
376 learned more and got to about -6.75dB reduction, which is the most yet!

I tried learning for 50 repetitions, but couldn't do better.
Whatever I tried, it seemed to start overfitting at some point and did worse.

### Next steps
- Can we remove all latency?
- Directly output audio to the speaker, not through the Scarlett?
- Directly receive mic input? Analogue to digital conversion though?
- Different DSP board like SHARC?

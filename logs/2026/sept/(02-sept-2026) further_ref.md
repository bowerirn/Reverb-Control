I figured out that it wasn't the IR that caused the massive shift in distribution of energy across frequencies,
it was the ref itself that was shifted.
I realized that it might be because the ref mic was like right inside the speaker cone which might cause distortion.
So I measured the ref again but I moved it up centered with the speaker, and then again but moved further back.
It seems that moving it back definitely made a big difference.

Band                 old error centered     error     further error      old ref    centered ref    further ref
   0- 200 Hz             1.86%              1.84%          1.80%         22.24%         21.11%         13.52%
 200- 300 Hz            11.88%             11.85%         12.17%         48.03%         46.05%         33.52%
 300- 700 Hz            62.39%             62.10%         62.89%         26.98%         28.68%         41.89%
 700-1000 Hz             8.77%              8.89%          7.74%          0.91%          1.16%          2.72%
1000-2000 Hz             8.49%              8.56%          8.64%          1.36%          2.11%          5.56%
2000-5000 Hz             6.62%              6.76%          6.76%          0.48%          0.89%          2.79%


I tried cancellation, and I was able to get slightly better results, but nothing significant.
I thought maybe the error mic was suffering from the same issue as the ref mic being really close to the panel.
I recorded the IR at 4.5cm, 11cm, and 27cm. 
The 27cm definitely had higher energy percentages at the higher frequencies, but it didn't seem to be significant.
I want to try normal LMS to confirm, but I don't think the extra delay is worth the small improvement.
It would be an additional 63 samples, about .65 ms.

So the current issue is still that we converge after like 2 iterations then it does worse.
I want to try a couple of things:

1. Normal LMS. Maybe the secondary path filter is too skewing to the energy distribution
2. EE-FxLMS again. Maybe with a better ref it will be better.
3. Freeze the weights after 2 iterations and see how they work for a longer cancellation
4. See the affect of different leak values on longer runs
5. Schedules for mu and leak

I already implemented schedules for mu and leak. 
I need to test them still though.

### Next steps
- Run those 5 experiments, they're all pretty cheap
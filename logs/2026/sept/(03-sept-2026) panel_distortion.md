Before I proceeded, I wanted to check the power distribution for a recording from the panel at different distances.
I played the arthur clip at 4.5cm, 11cm, and 29cm.
Moving further away definitely increased the relative energy in the higher frequencies, 
although it wasn't as drastic as with the ref mic.
11cm seemed like a good balance, it added about 18 samples of delay while helping a lot in the 700-1k Hz range.

Band               error 4.5cm     error 11cm    error 29 cm
   0- 200 Hz            10.64%          5.24%          2.77%
 200- 300 Hz            37.12%         30.68%         43.35%
 300- 700 Hz            46.77%         51.91%         38.05%
 700-1000 Hz             1.39%          3.93%         10.90%
1000-2000 Hz             3.66%          7.39%          4.44%
2000-5000 Hz             0.42%          0.85%          0.48%


I measured the IRs at 11cm and tried cancellation.
Worth noting that I left the source/ref at their further distance, although I didn't measure it.
It was definitely more stable for longer runs, but nevertheless got lower dB.
Best I could get was like -2dB even over like 10-20 iterations.
That makes it worse than 4.5cm regardless.
I did use 104 samples for lag and also swept the lag to find optimal. 
Sometimes it was like 90, sometimes 104.

I tried with IR len 256 and 512, both normal and EE IRs.
I can pretty confidently say that 256 is better than 512 now for IR length.
I did several grid searches in every setting and 256 always outperformed 512,
even accounting for extra time to learn given the extra information.
Also in this case I the EE-IR struggled a lot, I'm not sure exactly the reason for it.
But we know now that the error mic should be closer to the panel.
I still want to run those straightforward experiments now that the setup is (maybe) optimized.

### Next steps:
1. Normal LMS. Maybe the secondary path filter is too skewing to the energy distribution
2. EE-FxLMS again. Maybe with a better ref it will be better.
3. Freeze the weights after 2 iterations and see how they work for a longer cancellation
4. See the affect of different leak values on longer runs
5. Schedules for mu and leak
I ran the experiments from before:

- I tried using normal LMS instead of FxLMS, and it didn't really work. I was able to get a small degree of cancellation, but nothing like with Fx.

- I tried freezing the weights after 2 iterations. It still maintained a high level of cancellation.

- I did a sweep of leak values with 10 iterations. 1e-6 was still the best.

- Schedules for mu and leak definitely helped with stability, and slightly with performance too.

- I tried IRs of 196, 256, 384, and 512, 256 was the best.


I was able to get a 10 iteration run with like -4.3dB which is promising! 
but it still isn't what we want.
Dr. Duan recommended the following tests:

1. Try using a y cable to feed the clean source directly into the sharc instead of using the ref mic recording

2. Try moving the ref mic to the side of the speaker to see if the energy distribution skew is caused by reflections.
    - Possibly try smaller mics to reduce that?

3. Try other techniques to measure the IR

4. Make an offline simulation that matches the online algorithm exactly, and simulate that to see if it's the algorithm or the recordings themselves.

### Next Steps:
- Get an accurate offline simulation running
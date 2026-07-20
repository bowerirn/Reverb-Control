To test the delay from the SHARC I made the callback run fxlms but then just pass through the ref to the panel anyway.
Then I moved the speakers apart and turned them away from each other to test just the system delay, no acoustic delay.
I played a click from the source in python and recorded it from the error mic, 
then compared the ref and error mic peak indices.

Out of the gate, there was a lot of distortion. 
When I commented out the fxlms processing, it went away. 
Apparently there's a fixed time to process each buffer, 
and when you don't finish then there's distortion. 
I was able to turn on compiler optimizations, and that solved the problem

Then, the delay was like 5.3ms or 258 samples. 
Still better than before, but way slower than we expected. 
I found that I was able to turn off core 2 and make the buffer size smaller. 
I tried all combinations

| Block size | 32 | 16 | 8 |
|------------|---:|---:|--:|
| **2 cores** | 258 / 5.375 ms | 179 / 3.729 ms | 138 / 2.875 ms |
| **1 core**  | 193 / 4.021 ms | 146 / 3.042 ms | 122 / 2.542 ms |


So we can cut the default delay in half, but we were still hoping for less. 
Either way, this allows the speakers to be 1m apart safely, and the filter should be able to absorb the delay.

### Next steps:
- Try cancellation
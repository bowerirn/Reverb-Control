I looked into the buffer size instability for the no cancel run. 
It seems that anything less than a buffer size of 64 gets distorted, although 48 takes a bit longer than 16.
I listened to the learning run at a buffer size of 16, 
and realized that the distortion was from the source, not the cancellation wave being too loud. 

I thought maybe my callback was too slow, so I optimized it.
Most of the optimization here was preallocating blocks to copy the `indata` and `outdata` into each block.
This didn't fix it. 
I tried playing audio from spotify through the Scarlett instead of from Python. 
This also distorted with buffers smaller than 64, so the issue is the Scarlett I guess, or maybe my laptop is just slow. 
This means that the minimum buffer size I can use is 64, which is about 376 samples of latency (~8ms). 
Not ideal, but certainly better than 4312 samples.

### Next steps
- Try the delta seeding grid search
- Learn filters, see if less latency helps
I tried rerunning the sweep of block sizes, this time only with the ones allowed by the focusrite.
All of this was with 0.0 latency from `sounddevice`.
For whatever reason, my Focusrite drivers still bricked again and again after this sweep. 
My guess is that rapidly switching the buffer size then playing audio at that buffer size makes it upset. 
I added a pause of 10 seconds between each buffer size test, and also used the garbage collector between. 
That seemed to fix the issue, and I was able to make it run. 

As expected, latency went down with block size. 16 is the smallest available and had a latency of like 235 samples.
I tried learning some filters, and the audio was a little distorted when I ran the no cancel.  
I did the no cancel run at a block size of 64, and switched to 16 for learning and that seemed to be fine. 
I was able to get some results, but it seems like the cancellation wave gets too loud, so maybe I need to turn down the cancel gain.

### Next steps
- Why is the no cancel run distorting?
- Maybe run the delta seeding sweep again?

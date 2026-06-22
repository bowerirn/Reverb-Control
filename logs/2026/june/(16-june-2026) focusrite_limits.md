I found out that `sounddevice` has a latency setting that can be adjusted. 
The default setting has to be changed for playback, or it can be passed into streaming. 
It can be a float of the latency in seconds, or 'high' or 'low'. 
It defaults to 'high', which is about 90ms, exactly what we saw. 
I found that 'low' doesn't change anything from 'high', so I set it to 0.0.

I did a sweep over block sizes from 2048 down to 1, and also 0 because `sounddevice` mentioned that as an option.
The latency changed with block size this time, with lower block size having lower latency. 
After this sweep though, playback through the speakers became distorted. 
I did some isolation testing, and found that the issue was specifically when I played audio through the Scarlett. 
The amp and speakers were fine taking direct audio. 
The Scarlett audio was even messed up through the headphone out port, so the issue was software most likely. 
I reinstalled my Focusrite drivers, and that seemed to make the issue go away.
It seems like there are a set number of buffer sizes the Scarlett can use, and I probably bricked my driver trying to use a block size of 1.

### Next Steps
- Build in safeguards for different buffer sizes
- Explore how low I can make the latency
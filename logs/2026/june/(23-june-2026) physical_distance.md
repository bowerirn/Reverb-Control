I tried putting the panels as far apart as possible.
The idea is that in order for it to be causal, the source needs to take longer than the hardware latency to reach the error mic.
I'm not sure if I was able to get it quite far enough away, but it should at least be closer. 

I also found that a block size of 64 is not entirely stable.
Turning the safe mode off is definitely unstable for block size 64.
Block sizes 64-96 have clicking noises, although the audio doesn't devolve to total distortion.
128 seems to be the first completely stable block size.

I ran into some issues with panel distortion when increasing the gain on the focusrite.
For whatever reason, when I played a pure sine, it sounded like there were multiple tones
After some debugging, I just replaced the driver and that seemed to fix it.
The big panel was distorting too though, it might need a new driver as well.

Next steps:
* Delta seeding grid search
* Is cancellation better?
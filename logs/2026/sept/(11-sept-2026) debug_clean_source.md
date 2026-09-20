I tried playing around with the setup to give the clean source to the sharc.
Switching the direction of the y cable didn't do anything, 
and I also tried using a 1/4 to 1/8 adapter to use the line in port on sharc that I was using before instead of audio in.
None of that worked though, it still wasn't getting to the update step.

Next I tried putting the single end of the y into input 2 on the scarlett and playing stuff.
I tried it both ways, and it seemed to correctly pick up one of the inputs each time.
I'm not sure which is L or R, but correct signal is at least getting through 1 channel of the y.

Using audio in, I tried using the default out = in passthrough code on the sharc.
Playing some music through the source speaker, I was able to feel vibrations on the panel, so that seems to work.
But there still isn't enough signal strength for the update step.

Switching back to FxLMS, I made it keep track of the largest ref value it saw and inspected after.
Both ways for the y, the max value was like 1e-3. 
The max amplitude of the source digitally is like 0.8.
I tried line in with the adapter, and there was a single run where the max value was like 0.1.
But still no update happened, and I couldn't replicate it with either y direction.

Dr. Duan suggested y-ing in the source after the amp instead of after the laptop.
I still need to try this.
Dr. Heilemann suggested swapping the two mics to check if the energy distribution is still skewed, 
and he also gave me a well calibrated mic to test as well.

### Next steps
- Y the signal after the amp
- Debug the offline simulation
- Swap mics and test the good mics
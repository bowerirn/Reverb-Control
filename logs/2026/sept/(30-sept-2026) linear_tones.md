I implemented machinery on the SHARC to start and stop tones with given frequency and amplitude.
I didn't have an extra channel to send a new type of message, so I shared the sine sweep channel.
The idea is that they both need and on or off switch, 1 bit.
The sine sweep doesn't need anything else, so if everything else is 0, it's a sine sweep.
It doesn't make sense to play a tone with 0 frequency or amplitude.
Then for the tone, I used the 6 bits from the controller with the on switch for the amplitude, and the 7 bit value for the frequency.
The amplitude was just that value divided by 100, so 0-.63. The frequency was multiplied by 100, so 0-12.7kHz.

Using that I did a grid search of frequencies and amplitudes. 
The response amplitude was constant across frequency, and scaled with the input amplitude.
This strongly suggests that the response is linear, and an FIR is sufficient to represent the IR.

I also looked into the sine sweep we generated on the SHARC.
I had it doing a *= itself each sample to build the exponential, instead of explicitly computing the closed form exponential each sample.
They should be equivalent, but with the first way, float errors would accumulate and explode over time.
So I changed it to used the closed form calculation and keep track of its time index instead, which fixed the problem.

However, the SHARC IR didn't really seem to change, and was still negated and had the larger spike.
It also sounds different when I play it, so I need to figure out what's going on there.
I also think I'm not normalizing it the way I do the SCARLETT IR, so I need to look closer into my method here.

### Next steps:
- Make sure the SHARC IR is exactly equal to the SCARLETT IR
- Different IR measurement methods?
The only thing that seems inconsistent between the offline and online simulation is the actual control path.
So perhaps we aren't properly modeling the panel IR.

The first thing I tried was in simulation, shifting the update IR by a range of samples to see how it performed when the physical model IR was fixed.
The plot shows that the performance is very sensitive to these shifts.
So perhaps the way we align the IR by mic distance is wrong.
I'm not sure how else to do it though since the deconvolution puts the peak near the middle of the signal.
Perhaps another IR measurement method could work better for preserving alignment.

I also wanted to see if the Scarlett output and the SHARC output differed appreciably.
ChatGPT had my play a few pure tones at different frequencies through the SHARC one at a time and calculate the amplitudes from the error mic recordings.
It then converted those to dB scale, and compared it with the dB magnitude of the DTFT of the panel IR.
They differed by quite a bit, and particularly at 1kHz there was a lot more energy from the SHARC tone than the IR predicted.
So perhaps the Scarlett IR is not a good match for the SHARC output.

ChatGPT helped me make a streaming version of the sine sweep I use through the Scarlett so I could run it on the SHARC.
I don't think it's exactly the same as the scipy chirp function, but it produced an IR that looked quite similar.
There was larger magnitude overall, especially at the peak, and there were phase differences as well.
The polarity was reversed though.

I ran a 2x2 comparison in simulation of the update IR and the physical model IR using the 4 combinations of the SHARC and Scarlett IRs.
I had to negate the SHARC IR to be the same polarity as the Scarlett IR to make the mismatched cases not explode.
But once I did that, all 4 of them worked. 
They were all between -9.5dB and -11dB, so it seems even this mismatch shouldn't stop the real time algorithm.

I tried using the SHARC IR for real time cancellation, but it just didn't work.
I think I can put in some more work to verify that it's measuring correctly, but I don't think this alone will solve the issue

### Next steps:
* Poster for GIDS-AI presentation
* Different IR measurements


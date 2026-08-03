I was struggling a lot to get any results with physical cancellation, everything was inconsistent.
I realized that over time, the cancellation gets worse.
I did an experiment where I played the Arthur clip right after the program loaded, and it was like +2dB. 
Then I let it sit for a minute, and when I played it again it was +8dB. 
I let it sit one more time, and then it was +13dB.

To fix this, I used a threshold.
If I just did an if statement on the reference value, I would lose the pauses between words in speech.
I tried making a plot of the xnorm for NLMS, but it had the same issue.
512 is nothing on the scale of a 96kHz recording.

Instead I used an exponential moving average. 
`y[n] = a * x[n] + (1 - a) * y[x - 1]`
I set the characteristic time to be 100ms, which was `a = 0.999896`.
This was able to capture pretty much all the speech at a threshold of `3e-4`.

Using this, I was able to get cancelation to -1.82dB reduction, so it is finally starting to work!
Tuning it is painful though, I have to manually change everything between runs, and I can't inspect the `w`.

### Next Steps:
- Figure out how to send and receive data from Python to the SHARC
- Write a protocol for updating the params from Python
- Figure out how to log the weights
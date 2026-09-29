I resolved the lagged L2 regularized least squares with the current recordings and got the filter coefficients.
Freezing them in cancellation, it got exactly the same dB reduction as the LS predicted.
So at least the control generation path seems to be correct.
I tried learning from the LS filter though, and it started doing worse.
This suggested that the issue had something to do with the update.

I tried plotting the xf buffer after it had gotten worse, and plotted the expected xf buffer by filtering the ref sefgment with the IR.
This revealed that the last chunk of the xf buffer was not what it should have been.
Apparenltly I made the xf buffer only the filter length, so the lag portion looped around to the back, which caused the mismatch.

I fixed this lag bug by making the xf buffer be the filter length + lag to maintain the full history we need.
This immediately fixed the simulation, and it was able to get like -9.4dB with the same params.
I went and fixed the bug on the SHARC code too, although I have to test it still in real time.
I think the bug snuck in when I put the code on the SHARC because I added lag after getting the core working.

### Next steps:
- Test bugfix in lab
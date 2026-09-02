Ok so apparently it's possible to fry the ICE JTAG boards with the SHARC.
I'm not sure exactly how I did it, it most likely either overheated or something shorted.
From now on I'm gonna have it sitting on the bubble wrap bag it came in, 
and I'm not gonna leave it running straight for more than like an hour.

Once I got this working, I tested the parallel weight transfer, and it makes a beeping sound.
I realized 2 things:
1. Doing something like `while (!send_weight(w[i]));` 1024 times while also processing audio is probably throttling the core
2. When transferring the weights, I turn anc off which is definitely not what I want while cancelling.


I made a version that doesn't turn off anc, but I didn't address the nonblocking yet.
Instead, I made a request for only the wnorm, and I made that both blocking and nonblocking.

When testing the nonblocking versions, both timed out.
This was because of a fake spinlock I made to avoid pulling weights during an update.
I got rid of that and it stopped timing out.
While this means the weights/norm aren't always coherent, I hope it's good enough.
I might need to make it so that it copies the weights, then sends them off that copy.
For now I just used the wnorm request though.

The delta seeding experiment yeilded lackluster results. 
Perhaps this is because at double the sampling rate, each individual sample has a smaller effect.

I tried a 20 iteration long cancel run, and it got -3.35dB, with a converged filter norm of ~2.5.
The learned filter looked pretty reasonable, there was no clipping down to 1.
I froze it and ran nonadaptive cancelling for 5 iterations, and it got -3dB, so still pretty good.
It looks like it cancels the negative amplitudes really well but not the positives so much.
I don't really understand why it's skewed, or why we can only get in the -3s for dB reduction.

The next thing I tried was removing the panel and running cancellation to get a baseline for the best possible result in the room.
With the stand still there, it was 9.78dB. With the stand removed it was -9.92dB.
So clearly the panel induces a lot of reflections, and there should still be a lot of room to cancel.

I saved those, along with recordings from the error and ref mic during no_cancel so I can try to solve for the optimal filter.
Then I can try implementing something to send the weights from python to the sharc, freeze them, and use them.


### Next steps:
- Analyze the data offline
- Make a way to seed a filter on the sharc
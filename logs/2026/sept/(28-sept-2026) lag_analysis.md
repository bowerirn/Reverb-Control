Long story short, the lag relationship we're using right now is correct.
It doesn't make sense to delay the error and the ref more for the current control signal.
We use the current error mic to update the weights, and the ref is `system_delay` samples ahead of that.
So the way we were doing it is already correct, it's not the lag.

I did realize that I wasn't adding the ref feedback to the ref_nc signal in the offline simulation though.
I added that back in, then compared to without adding it, both with cleaning the ref and not.
They all did like within .1-.2dB of each other.
That makes sense because the ref IR is basically 0.

So now we're back to the model of the secondary path being the only culprit.

### Next steps:
* Compare pure tones with different amplitudes
* Feed the pure tones directly into the Scarlett
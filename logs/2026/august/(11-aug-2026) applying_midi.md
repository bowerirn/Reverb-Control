Ok a lot to fix with the midi setup.
First off, when the scarlett is plugged in first like after measuring the IRs,
the SHARC midi numbers are 1 and 2 instead of 0 and 1.
I should automate detecting this, but I just allowed passing in the specific devices.

Sometimes, the midi just doesn't work. 
It either doesn't send messages, or sometimes Windows can't recognize it.
In either case, I have to restart my laptop to get it working.
I used LED12 on the board to detect when it was receiving MIDI messages so I don't have to guess.

I was also having issues with the weight transfers. 
I could get the expected length, but not the weights.
I have limited lab time this week, so I decided to shelve that for now and skip weight logging.
In reality, I'll need a better solution for the weights.
If I pause each iteration to get the weights, then it doesn't sync well with the no_cancel baseline.
This is a future problem though.

Once I got all that working, the grid search seemed to work fine.
As the `cancel_gain` and `mu` increased though. 
If I did the same hyperparams for 5 runs, I might get 1-4 that fail catastrophically while the rest cancel.
I think the issue is that the ref mic gets feedback from the panel which causes it to go out of control.
My idea is to use a foam shield or something to block those waves.

### Next steps
- Make a foam shield for the mic
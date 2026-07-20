With help from ChatGPT I implemented FxLMS in c++ again for the SHARC. 
The code is mostly similar, but it takes just 1 sample, and the reference filtering is done in c++ now too. 
The reference filtering is an FIR which is basically just a convolution. 
We maintain another reversed ring buffer of the reference the size of the IR instead of the filter. 
Then it's just a dot product with the IR instead of with the filter. 
Kind of interesting that it's literally the same thing, that a convolution is just a dot product with the reverse of one signal.

I wanted to be able to at least turn cancellation on/off from python.
Me and ChatGPT went down a whole rabbit hole of trying to use a uart connection to sent bytes.
We were using the midi callback to read the uart port instead of midi. 
Apparently our SHARC board doesn't support uart connection though, so we had to scrap the idea.
We made all the knobs for fxlms be public so we can change them in the debugger.
It isn't ideal but it will suffice for now.

### Next steps:
- Test SHARC delay
- run fxlms
Instead of using the foam shield, 
I realized I could just use the IR from the panel to the ref mic that I measure for free.
The idea is basically I keep a buffer of IR_LEN with all the past control signals.
Right at the process() loop, we convolve the ref IR with the control buffer to estimate expected feedback for this sample.
Then we subtract that from the ref at this sample to get the cleaned reference.
We use that cleaned ref to get the new control, and put that in the buffer.
So it's basically the same process as with the other FIRs, 
except that we convolve before we update the buffer here, 
because we need the result of that FIR to clean the ref and generate the new control to put in the buffer.

This helped, although there would still occasionally be catastrophic failuers.
I did some logging, and realized that it was a result of the `xnorm` going negative on occasion.
The way I calculate the `xnorm` is `xnorm += new * new - old * old`.
This saves recomputing the whole dot product every time, 
but the subtraction allows float roundoff to go negative sometimes.
I just added an if statement to clip it to 0 if it's negative, and that fixed the issue.

This made everything stable, although performance did not really improve.
I realized that since I doubled the sample rate, I should also double the size of the IR and the filter.
Changing the IR length requires me resetting up for IR measurement though, I need to make a faster way to change it.


### Next Steps
- Make it easier to change IR length
- Double IR length and filter order
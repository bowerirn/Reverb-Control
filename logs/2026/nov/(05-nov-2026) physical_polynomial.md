I tested the IMD products with the panel itself, and it definitely produces some.
This means that an IR is not sufficient to model the secondary path, it is inherently nonlinear.

I used ChatGPT to implement the polynomial filter on the SHARC.
I tested it in the physical space, although results weren't great.
The original linear filter seemed to not be doing very well either.
It might have been a driver issue, although the IR looked fine, so I'm not sure.
I'll try with a fresh driver next time.

Polynomial did about the same as linear with the scuffed setup, maybe slightly worse.
This is to be expected, since it's the secondary path that is nonlinear, not the reference.
So somehow we're gonna need to replace the whole filter with some sort of nonlinear function of the ref and the current control.

### Next steps:
- Replace the driver and rerun tests
- Find other nonlinear secondary path ANC paradigms
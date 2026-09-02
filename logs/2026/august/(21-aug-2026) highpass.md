Dr. Heilemann suggested that perhaps the ambient room noise was causing issues.
He said that often times measurements from in the whisper room only look right after highpassing the signal.
The cutoff should be roughly 200 Hz.
He suggested something with a sharp cutoff like a 6th order Butterworth filter.

I implemented it in both python and c++ for the SHARC.
I stole the coefficients for the scipy biquad to use on the SHARC to make sure they're the same.
I started off by filtering the source and the IRs in python, the ref and error mic signals on the sharc,
and the ref and error logs back in python.

It didn't work. At best it did nothing, at worst it created feedback.
I realized that highpassing the IR and the ref meant doubling the filter, which was maybe problematic.
Since highpassing the IR changed the shape of the IR, I just used the plain IR.
This had less feedback, but still didn't do anything.
I need to debug and see if this is an implementation error or why the highpass causes problems.

### Next steps
- Debug the highpass filter
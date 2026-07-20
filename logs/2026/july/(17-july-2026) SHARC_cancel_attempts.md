Cancellation was very unstable at first. 
We got a lot of feedback and the speakers went crazy.
I turned down the step size, cancel gain, and increased the leak. 
I also swapped the L/R I had for the Scarlett direct monitor to ensure 
the panel played what came from the ref mic and not the error mic.

This fixed the problem, so I started trying to cancel audio. 
I used 122 samples of delay based on the table from yesterday.
I wasn't really able to get any results though.
Later I realized I had 32 buffer size dual core on, so the delay was wrong.
Nevertheless, we think there's room for improvement on the delay, so that's what we're gonna target next.
At the very least we can maybe do a buffer size of 4. 
But I'm gonna try and remove as much bloatware as I can.

### Next steps
- Measure the delay purely from the SHARC board
- Remove as much bloatware as possible from the SHARC
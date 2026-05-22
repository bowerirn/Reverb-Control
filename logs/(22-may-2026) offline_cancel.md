I gave up on the loopback delay estimation of the IR, and just measured the distance in cm to align it. 
I took 8 repeated measurements of the IR with the same settings and plotted them. 
They aren't perfect, but they seem to be pretty close, and the random Matlab time shifts are gone. 

I tried running through the whole cancellation routine. 
When I played the cancellation wave, it didn't sound as delayed as it did in Matlab. 
There was also not nearly as much energy added to the signal as with Matlab. 
That said, we did worse with the cancellation wave than without (i.e. we aren't cancelling it).

It might be because the adaptive filter is better at cancelling the end of the signal. 
The next step is probably to try testing it in an online setting, although I don't know how to do that in Python yet. 
I also might try the toeplitz least squares method we did before just to see if it was a fluke or not.

### Next steps
- Figure out how to modify the signal online with Python
- Let the adaptive filter learn in real time
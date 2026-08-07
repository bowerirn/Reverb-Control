## The IN and OUT labels on the midi cable refers to the device the USB is plugged into, not the midi plugs.
### So the IN plug goes into MIDI OUT on the SHARC, and OUT goes to MIDI IN on the SHARC.

That cost me more time than I care to admit lol.
Anyway, once I sorted that out, I was able to send messages.
ChatGPT helped me write a suite of python functions to implement the protocol using the `mido` package.
I turned off the LED11 heartbeat so I could use both 11 and 12 for debugging.
All of the receiving values worked properly.
I just used a conditional in each receive case to ensure it got the right value (within 1e-7 if it was a float).



### Next steps:
- Test receiving the weights
- Build a python class to manage the midi protocol

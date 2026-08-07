For a while I was unable to receive values from the board.
Instead of transferring weights, I was just sending a single message back when it got the request.
But it wasn't working.
I tried the MIDI THRU on the SHARC to see if it was maybe the output cable not working.
MIDI THRU should just bypass the code and send the input back.
That didn't work.
ChatGPT suggested that maybe the firmware doesn't support MIDI_THRU.

I tried restarting my laptop and replacing the midi callback with the original passthrough one.
This worked.
I went back to my custom code, and it started working.
I'm really not sure what was going wrong.

Transferring all the weights wasn't working though. 
It seems that you can't just do that much midi sending in the callback.
It queues it up in the buffer, and I guess that buffer is cleared outside the callback.
So it just waits for more space and is caught in deadlock.

Instead I set a flag when we request a weights reset or a weights transfer.
I wrote a backgound loop that just checks if the flags are set and handles them.
I call that loop from main() alongside the process_audio_background_loop().
This worked, and I was able to transfer the weights!

___

Next, I wanted to optimize the float representation.
Currently I use 7 bits to encode the positive integer, and 7 bits for a negative exponent. 
I don't actually need integer values that high though, scientific notation only needs 1-9.
In practice, I only need exponents down to about 10^-12, so 4 exponent bits are sufficient.

I used the other 10 bits to quantize the significand. 
The significand is represented as an integer in the range [100, 999],
corresponding to decimal values in [1.00, 9.99]. 
During decoding it is simply divided by 100.
This provides two decimal places of precision.
So I can use values like 1.25e-6  ->  significand = 125, exponent = 6
The all-zero payload is reserved to represent exactly 0.

controller: [negative exponent: 4][significand bits 9-7: 3]
value:      [significand bits 6-0: 7]

___

Lastly, I put all the python machinery into its own class that the SHARC class can use, 
so I should be ready to grid search in the lab.


### Next steps
* Grid search in the lab
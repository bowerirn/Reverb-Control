Today I designed the midi protocol and implemented it in c++.
There are a few different kinds of values we need to send:
1. Small floats in [0, 10), such as the step size or cancel gain
2. boolean flags (turning on/off updating, etc)
3. positive integers (lag, mavg characteristic time)
4. requests (resetting the weights, pulling weights from the board)
5. filter weights

I used midi control change messages. These are 3 bytes as follows:

        1101           |         xxxx
control change opcode          channel

0xxxxxxx: controller

0xxxxxxx: value


For sending updates to the SHARC, I use the 4 channel bits to determine what we want to update. 
That leaves 14 bits for the actual value.

For pulling data from the SHARC, I only used 4 messages.
1. Initiate weight transfer (with the number of weights to expect)
2. Weight value
3. End weight transfer
4. Status


### Encoding requests
This is the easiest, we don't need the controller or value

### Encoding booleans
This is easy, the value is just 0 or 1
* I treat the update sign as a boolean, and 0 just maps to -1 while 1 maps to 1

### Encoding positive integers
Generally none of the integers I use are big enough to use 16 bits, we just chech that it's less than 16383.
Then we split it into the high 7 bits and the low 7 bits.
The controller carries the high bits, the value carries the low bits.
In decoding, we just concatenate them and cast to an unsigned 16 bit int.


### Encoding small positive floats
I use the fact that most of these floats are values like 3e-6, or 1e-4, etc.
We have a small positive integer, and a negative exponent.
We don't need anything for the sign since we know them.
To encode, I just manually split the exponent and mantissa.
To decode, it's just mantissa * 10^-exponent  


* For the moving average weight, I send the characteristic time in ms as an int, and calculate the weight on the SHARC, since it has a lot of decimal places.



### Encoding filter weights
In general, none of the weights should be > 1. 
The norm of w should be <= 1, and in practice even when we seeded deltas they would drop lower than 1.
Knowing this, I just used a fixed point quantization.
Since the message only needs 2 channel bits, I used the additional 2 to get 16 total bits.
1 bit for the sign, 15 bits to quantize.
Basically we just do round(w * 32767) to encode, and divide by 32767 to decode.
This is technically not pure q15 because signd ints go from 32767 to -32768.
If we used 32768 though then +1 wouldn't be a true value.

channel:    [message type: 2][sign: 1][Q15 bit 14: 1]
controller: [Q15 bits 13-7]
value:      [Q15 bits 6-0]


## Next Steps
- Test with the LEDs
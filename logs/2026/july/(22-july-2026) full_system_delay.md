The SHARC audio processing works as follows:

// ADC reads to ping in buffer         |      CPU processes pong buffers
// DAC writes the ping out buffer      |

Once the ADC fills the ping buffer, 
context flips and the ADC/DAC work on pong, and CPU processes ping.

So in principle the processing should take:
* `block_size` samples while ADC reads in
* \<context flip\>
* `block_size` samples while CPU processes and ADC reads in other buffer
* \<context flip\>
* DAC starts writing it out

The total should be `2 * block_size`.
But the CPU stage looks like this:
* Write previous output buffer
* Read current input buffer
* Raise interrupt to signal callback

But by working with the previous output, we delay it a full `block_size` samples.
This is why we see `3 * bock_size`

I rewrote it to do:
* Read current input buffer
* Directly call callback
* Write current output buffer

As long as that all can finish in `block_size` samples worth of time, it saves a cycle.
Now it takes `2 * block_size` samples.


The next thing I wanted to measure was the direct monitor (DM) path for the Scarlett.
Supposedly it should just be a wire, 
but there are a lot of samples still unaccounted for after measuring the SHARC delay.

I had out1 wired to in1, the Scarlett set up to pass in1 through DM1, 
and I wired DM1 to in2.
At 48kHz the click was delayed by 29 samples through the DM.
I also tried it at 96kHz, and it was 19 samples.
I'm not sure if there are other ways to lower this.

I hooked up the whole system to measure the delay through the speakers again.
At 96kHz, the delay was about 94 samples.

* The SHARC with `block_size = 4` and 96kHz took 46 samples
* The Scarlett DM at 96kHz took 19 samples
* The error mic was 4.5cm away which at 96kHz is like 12-13 samples.

In total that's like 77-78 samples, so we still have 16-17 unaccounted for.
I didn't measure if the amp has any delay.
ChatGPT thinks perhaps the IRs are delaying the peak arrival by 10-15 samples,
and maybe that's the discrepancy we see. 
This is less than 1ms though which means 34.3cm is enough distance for the speakers.
Usually they're a bit further, so the delay is dominated by that distance now.

The FxLMS callback does not finish in time though, so I need to optimize that

### Next Steps
* Optimize the callback to fit in time
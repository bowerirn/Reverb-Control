To test the delay of the SHARC itself, 
I ran output 1 from the Scarlett through the SHARC, back into input 1 on the Scarlett,
And output 2 on the Scarlett to input 2 on the Scarlett, to compare the relative delays.

At 48kHz, the fitted model seemed to be:
`Delay samples = 50 + (3 * block_size)`

I did some digging into why.

The SHARC audio path looks like this:

Mic --> ADC            DAC --> Speaker
      |                    ^
      |                    |
      v                    |
   SPORT RX            SPORT TX
        |                 |
        |                 |
     DMA RX            DMA TX
        |                 ^
        v                 |
       RAM <----CPU----> RAM

According to the [datasheet](https://www.analog.com/media/en/technical-documentation/data-sheets/adau1761.pdf)
The ADC always takes ~23 samples
The DAC takes 25 samples @ 48kHz, or 11 samples @ 96kHz.


So the 50 samples was ADC + DAC + 2 unexplained samples.

Bumping it up to 96kHz reduced the number of samples, 
and made the samples each faster, making it a win-win overall.
The new model became:
`Delay samples = 23 + 11 + (3 * block_size) + 4`

So it looks like those unexplained samples are time based, ~42us.


Based on the fact that the Baremetal framework uses ping pong buffers, 
the `3 * block_size` term should in theory be able to be reduced to `2 * block_size`?
Also, people on the forums suggest that you can get a sample by sample system going 
if you don't use Baremetal and write your own DMA.
I'll look into the simple optimization first, at a block size of 4 we're dominated by the ADC and DAC.
We also still have a lot of unexplained samples of delay left, 
for block size 8 and 48kHz SHARC was only 74 out of like 122 samples.


### Next Steps
* Optimize the DMA handling on the SHARC
* Isolate the other delay components in the system

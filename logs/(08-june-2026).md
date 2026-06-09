The cpp extension works. 

To look at the delay I tried several things. 
First, I plotted the raw IR from the playback mode vs the raw IR from the streaming mode. 
They appeared to be the same, and the aligned versions were equivalent, so no delay there. 

ChatGPT suggested that the buffers for the input and output in the callback for streaming are not time aligned. 
I printed out the difference, and it is consistently off by about 4312 samples. 

I also plotted the coherence of the error recording with no cancellation to the source. 
Out of the box, it's basically 0, which explains why the filter didn't learn anything.
There wasn't really anything to learn. 
Plotting the coherence when delayed by 4312 samples yielded better results, but not perfect. 

I tried doing a cross-correlation between the error recording and the source to find the true delay, and it was about 4470 samples.
I think this extra delay is from the physical path, since ~160 samples aligns with our earlier delay experiment on the IRs for the reverb IR.
I can just calculate this lag when I run the no-cancel job at the start for reference, and use that lag with cancel.

I made the cpp extension support a delay on the filtered reference and the raw reference. 
The weights are learned off the filtered reference, but the output is generated off the raw source,
so there might be a different delay needed for each of them since the filtered reference is delayed from the IR. 
For now I'm leaving them both as the value from the cross correlation, just to see if anything happens, since it might be close enough.
but in the future the raw reference might need to be adjusted by the delay in the panel_ir.

I also made the cpp extension support a sign for the weight update to control whether it's a +=/-=, so I can test to see which is correct.

So far I haven't gotten any results yet, most of the nothing has been motivating the updates I made.
I just need to do more experiments now that I have more hyperparams to play with.


### Next Steps:
- Grid search the cancel_gain sign and the weight update sign to find the correct combination
- Search over different amounts of delay to see if the filter starts learning
- Try separate values for the raw source delay and the filtered x delay
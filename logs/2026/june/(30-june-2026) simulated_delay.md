I added simulated delay to the offline FxLMS program. 
I added 0 padding to the start of the filtered reference. 
In principle, the error mic sees 2 paths:
1. [source -> reflect off panel -> error mic]
2. [source -> ref mic -> hardware latency / processing -> panel IR -> error mic]

The first is the target. 
The second can be accomplished by adding delay to the panel IR.
I was worried about modifying the IR itself though, so I added the delay to the filtered reference.
This would be the same result as filtering with a delayed IR.
I also accounted for not delaying the IR when getting the prediction signal we use for MSE after every iteration. 

The results show that the delayed performance is an order of magnitude worse than not delayed.
I used a 1024 length filter for a 128 length IR and 376 samples of delay, so in principle the delay could have been learned.
The scales were off since I used old data because I'm out of town, but I think it's clear we need to confront this delay issue.

### Next steps:
- Can we get the SHARC set up physically?
- Program the SHARC to do FxLMS
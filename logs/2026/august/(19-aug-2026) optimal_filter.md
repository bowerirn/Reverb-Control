I solved for the optimal filter with least squares. 
It was able to get -18dB, but the norm was over 15,000
So I added on L2 regularization on the weights and swept the lambda coefficient on the regularizer.
At lambda .01 it was still just short of -18dB, but the norm was only 2.3.
I chose this filter since it was most similar to what the model learned.

I made it I can load weights into a seed file, and then load them into the sharc filter with a midi command.
I tested out this optimal filter, but got no special results.
I thought maybe I had missed a negative in the cancel gain so I tried both - and + gains, and swept from 0 to 1.
Negating it did worse, and the best it did normally was like -0.6dB which is nothing.
I don't really understand what is going on.
I did notice that the peak of the optimal filter is later in time than the peak of the learned filter.
I need to investigate this further.

### Next Steps
- Investigate lag and filter shape
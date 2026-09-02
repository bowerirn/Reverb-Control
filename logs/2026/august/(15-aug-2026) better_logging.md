I made a couple of changes based on the results from 2 days ago.
First, I allowed delta seeding to have more than just +1 value.
Since I don't expect to use a filter more than length 2048, I used only 11 bits for index.
That left me with 3, so I used 1 for sign and 2 for value.
The 4 values I chose are (.25, .5, .75, 1.0), and +/-.
I don't need to get super fine-grained for it, so that should be enough to play around with.

The next thing I did was spin up weight requests in a separate thread.
The idea is that this way I can just play the source to cancel, and pull weights on an interval.
The weights won't be exactly correct or coherent, but I hope they are close enough to get an idea.
I may end up needing to copy the vector before trying to send it to maximize coherency, but I'll test it first.

### Next steps:
- Test background weight logging
- Run a better delta sweep
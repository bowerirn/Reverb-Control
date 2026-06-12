I tried using longer filters on 50 iterations of short arthur. 
With a filter length of 2048, I was able to get -3.7dB reduction. 
It doesn't really look like the tail of the filter does anything though. 
Without seeding a delta at index 104, it only got -2.5dB reduction.
Still respectable, but not as much. 
I think the spike is doing most of the work, and the extra filter is to clean up the response, but I need to test for this.

I tried searching over a range of filter lengths and block sizes. 
A filter of 4096 does worse than 2048, so clearly longer filter only helps to a point before it is harder to learn.
The best block sizes seemed to vary with filter length, so I assume it's probably random, and filter length is the main factor in improvement.
Basically the lesson is around 1024-2048 filter length with many iterations.
I still need to push the limit and try like 100 iterations. 

I also tried on the full arthur clip.
It was able to learn, but not nearly as much.
I didn't try with the seeded filter though. 
It does still seem to learn a similar filter shape, which suggests that this is actually a useful control filter.
I want to try learning on short arthur, then using that filter and cancelling full arthur. 
Worth noting that full arthur is also noticably quieter than short arthur, which might affect things. 

### Next steps:
- Seed a filter on full arthur
- 100+ iterations on short arthur
- freeze a filter from short arthur and cancel full arthur
- Possibly test 3 frozen cancellations after learning a filter to see if the seed is carrying:
    - full learned filter
    - seed filter
    - full - seed (i.e. the learned component)
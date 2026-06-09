I looked into changing the block size for cancellation.
I started with a big grid search over the other hyperparams in addition to block size. 
It became apparent that stable params do not transfer across block size.
Nevertheless, I wasn't able to find a configuration for any block size that did meaningfully better than 1024. 
I went all the way down to 32, so the only thing to try next in this line would be sample by sample. 
I probably need to use cpp to do sample by sample, but I don't think that will fix it. 
It's likely still worth implementing anyway though.

I also tried using the bigger panel, and playing with the distance between the speakers. 
It didn't really change anything. 
I confirmed that with both panels, there is sufficient energy from reflections when the panel is there as opposed to foam. 
I switched back to the smaller panel because the larger one sounds less stable.

### Next Steps
- Implement the callback in cpp and try sample by sample improvement
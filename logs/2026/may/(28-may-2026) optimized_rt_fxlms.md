Yesterday in the lab there were sometimes when the audio would be cutting weirdly, 
so I decided to optimize the adaptive filter callback to ensure it wasn't a processing issue.
I made the following changes:
- preallocated log arrays to avoid a bunch of list appends
- separate paths for lms and nlms updates to avoid checking in the kernel every sample
- running dot product for the nlms norm
- I updated w in place instead of reassigning it


### Next steps
- Grid search in the lab to find good hyperparams

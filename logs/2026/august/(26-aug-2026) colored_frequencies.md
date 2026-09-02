The delayed optimal filter didn't really work.
It got like -1.62 which is better than the non-delayed version, but still not great.

I started looking into the uneven frequencies, and it seems like this is a known problem.
This paper finds that convergence for FxLMS for different filter lengths is frequency based
https://www.sciencedirect.com/science/article/pii/S0022460X25000628


I found some papers that attempt to solve this problem

https://www.sciencedirect.com/science/article/pii/S0022460X03001500
* Prewhitens the system to reduce correlation in the reference signal


https://pubmed.ncbi.nlm.nih.gov/18537375/   
https://onlinelibrary.wiley.com/doi/10.1155/2008/791050
* Modify the secondary path coefficients to reduce eigenvalue varience


https://saemobilus.sae.org/articles/modified-fxlms-algorithm-equalized-convergence-speed-active-control-powertrain-noise-2015-01-2217
* They also do something with the secondary path to help with frequency dependent conbergence


https://www.sciencedirect.com/science/article/pii/S1051200415001499
* Supposedly a frequency domain block algorithm can more easily account for spectral energy imbalances 


https://www.sciencedirect.com/science/article/pii/S0888327015002770
* They use a wavelet based FxLMS algorithm that converges better


I'm gonna read through some tomorrow and maybe try implementing something.


### Next Steps
* Read some papers
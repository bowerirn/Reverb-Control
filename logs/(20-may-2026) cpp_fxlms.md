I started porting over the system to Python. 
Since I had experience making c++ Python bindings from my parallel computing class, 
I decided to implement fxlms in c++. 
Everything seems to be working properly in test_fxlms.ipynb, and it's fast. 
I tried both LMS, NLMS, and a scheduled step function. 
I'm not sure if the mic recordings and IR measurements are from the same time though, but the algorithm seems to work. 
It works on Windows for me too, although other users might need to install Visual Studio c++ tools. 

<br>

### Next Steps
- First, I need to make a way to generate sine sweeps and measure IRs, ideally with latency compensation.
- Then, I want to try setting up Python to directly interface with ASIO. 
I think I can use `sounddevice`, and this will save lots of time exporting to wavs and running in audacity. 
I might make a class to encapsulate all this. 
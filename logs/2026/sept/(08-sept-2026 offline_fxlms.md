I built the offline simulation to basically follow the sharc code.
To simulate the error mic, I fed the control signal through the panel IR, then added the error recording.
The reflections + the panel out should be the error mic.
Other than that basically everything is the same.
I did it sample by sample in python, even though it's slower, because I wanted assurance that every part is accurate to the sharc.
The only thing I sped up is using numpy dot products for the FIR helpers.

With 0 lag, it got -9.4dB which is like exactly perfect.
It probably would have gotten the max score if it didn't have to learn the filter from nothing at the start.
When I added lag though, everything exploded. 
It was getting to like 1e37, then reaching inf.
I still need to look into why this is happening.

### Next steps:
- debug the simulation
- record the mics at different speaker locations
- Try feeding the source into the sharc instead of the ref mic
I tested the energy frequency distribution with the mics swapped, and also with a very good, well calibrated mic.
This was the result:

| Band       | Source | Error pg58 | Error sm58 | Ref pg58 | Ref sm58 | Good mic |
| ---------- | -----: | ---------: | ---------: | -------: | -------: | -------: |
| 0–200 Hz   |   3.91 |       3.72 |       1.93 |    11.79 |    15.37 |    11.89 |
| 200–300 Hz |  16.22 |      15.31 |      13.05 |    30.42 |    33.78 |    26.31 |
| 300–700 Hz |  67.62 |      59.08 |      62.20 |    44.71 |    39.21 |    50.57 |
| 700–1k Hz  |   3.75 |       6.99 |       7.80 |     3.20 |     2.56 |     2.98 |
| 1k–2k Hz   |   6.80 |       8.49 |       8.10 |     6.25 |     5.43 |     5.60 |
| 2k–5k Hz   |   1.70 |       6.40 |       6.93 |     3.62 |     3.65 |     2.65 |


So the good mic definitely helps a little, but really not by much.
Swapping the mics also didn't do very much.
This suggests that it may be an issue with the speakers.
I'm not sure how much of an issue this is if least squares can find a solution though.

I looked into y-ing the source after the amp to use the clean source in the algorithm.
ChatGPT said I could definitely break something if I fed amp level signal into a line level input though, and I definitely shouldn't do it.
So I think I probably need a better 1/4 to 1/8 adapter if I'm going to try this going forward.

### Next steps:
- Why doesn't the offline simulation work?
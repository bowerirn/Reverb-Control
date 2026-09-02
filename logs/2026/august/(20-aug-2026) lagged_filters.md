I tried doing a sweep of lag but recording the filters after each time and plotting them on each other.
As lag increased, the filter was shifted earlier in time.
For whatever reason, 0 lag was unstable.
But it seems like it was trending toward the peak with lag 0 being in the same place as the optimal filter.
Granted, the lag only gets applied when updating the weights, but we did a lag sweep to find the optimal lag.
I guess the difference is the offline recordings only capture the acoustic delay and not system delay?
And then when we learn in real time, we have to account for system delay.

Anyway, I tried delaying the error mic offline signal by different amounts and then doing regularized least squares again.
In all cases, it was able to still get -18dB.
I unfortunately did not have time to test that filter today though.

### Next steps
- Test the lagged optimal filter
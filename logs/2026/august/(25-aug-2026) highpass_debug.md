I investigated the frequency spectrum of the arthur signal both highpassed and not highpassed.
It seemed to have much less energy <200 Hz with the highpass, as expected.
It appears that the 6th order butterworth filter causes significant group delay near the cutoff.
Like on the order of 500 samples, or 5ms. 
That's huge, and likely why the highpass didn't work.
An FIR filter would introduce linear phase distortion, but then the whole spectrum is delayed in phase which isn't ideal.
I could in principle use a lower order IIR and have less phase delay but a shallower cutoff.

Doing some more experimenting, here is what I (and ChatGPT) found

- We divided the signal into different frequency bands and integrated the spectral power
    - The <200 Hz frequency band only makes up about 2% of the total energy of the arthur clip
        - It's probably not worth highpassing for that anyway
    - The most energy is between 300-700 Hz, which makes sense for speech
    
- We took a dB ratio of the PSDs of cancellation and no_cancel.
    - The ANC really only cancels in that 300-700 Hz range
    - <200 Hz is garbage, and >800 Hz is basically no cancellation
    - The main question is why the controller is only working in the main band

- We looked at the PSD of the no_cancel signal
    - There is definitely energy in the >800 Hz range, although less than the 300-700 Hz

- We plotted the magnitude-squared coherence between the ref and error mics during no_cancel
    - The coherence was basically 1 above 250 Hz through 2 kHz, so we definitely can cancel the energy above 800 Hz

- We plotted the magnitude spectrum of the panel to error IR.
    - It wasn't perfectly flat, but it didn't have any sharp drops at 700 Hz or above
    - The panel is certainly capable of producing these frequencies predictably

- We plotted the IR measurement repeated multiple times
    - The IR is very stable between runs in both magnitude and phase, so IR instability between runs isn't the issue here

- We plotted the PSD of the ref filtered with the IR
    - Again, the energy is strongest between 300-600, the gradient seems to drop off by 700, which could explain why we can't learn those frequencies.
    - With NLMS, the xnorm is dominated by the frequency bands with the most energy. This likely drives the convergence in only the main region.

- We plotted the wnorm with and without leakage
    - Leakege definitely stopped the growth of the norm and forced faster convergence
    - With leakege, performance was better though, letting the norm grow freely created worse performance over time.

It seems like the issue is that it's only cancelling the frequency band with the most energy.
I don't know how to fix this though.
I realized that I never actually tried the delayed optimal filter in the real time setting,
so I'll try that then do some research.

### Next Steps
* Try the delayed optimal filter
* Research uneven frequency energy with ANC

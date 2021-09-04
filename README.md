# Harmonic-percussive-residual separation plug-in

This work is a study on the plausibility of a sines-transients-noise decomposition inspired algorithm as a real-time plug-in application. [iPlug2](https://iplug2.github.io/) framework is used to provide the main scaffolding. It allows building on multiple platforms and as relevant plug-in formats (AU, VST, VST3, AAX, ...).

![GUI screenshot](hpr-facelift.png)

The current state is a prototype and thus is missing some of the practical implementations usually found in plug-in code. The most important one is that the frame size is not handled and has to be handled externally. This might require having an appropriate sound card for the computer to adequately manage the buffer size as well as using DAWs that do not override those settings.

Thesis work based on the following papers:

- J. Driedger, M. Müller, and S. Disch, Extending harmonic-percussive separation of audio signals, Jan. 2014. [Online]. Available: https://www.researchgate.net/publication/303667409_Extending_Harmonic-Percussive_Separation_of_Audio_Signals

- L. Fierro and V. Välimäki, “SiTraNo: A MATLAB App for sines-transient-noise decomposition of audio signals,” Vienna, Austria, Sep. 2021, p. 9. [Online]. Available: https://www.researchgate.net/publication/354076466_SiTraNo_A_Matlab_App_For_Sines-Transients-Noise_Decompositon_of_Audio_Signals

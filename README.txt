This repository contains files relating to:
Slator PJ, Cairo CW, Burroughs NJ (2015)
Detection of Diffusion Heterogeneity in Single Particle
Tracking Trajectories Using a Hidden Markov Model
with Measurement Noise Propagation. PLoS ONE
10(10): e0140759. doi:10.1371/journal.pone.0140759.

Paddy J. Slator, Warwick Systems Biology Centre,
email: paddyslator@gmail.com.


The files in the main directory are:

SimulateTwoStateDiffusion.m
—————————————————————
Simulate the two-state diffusion model either with or without measurement noise.


TwoStateDiffusionExactMCMC.m
—————————————————————
Run the MCMC algorithm for the two-state diffusion model with fixed measurement noise (exact likelihood version, i.e. Equation 39 in main text). 


InitTwoStateDiffusion.m
—————————————————————
An initialisation script which, for a set of “reasonable” parameters, simulates the two-state model with measurement noise, and then runs the MCMC algorithm (1000 samples) on the simulated trajectory. (Using the simulated noise variance as the variance in the MCMC algorithm.) Also plots the relevant MCMC output. This file takes around 30 seconds to run on my four year old MacBook Pro (4GB RAM).

Also included is the ApproxModel folder, which contains files for the approximate likelihood model MCMC and marginal likelihood calculations.

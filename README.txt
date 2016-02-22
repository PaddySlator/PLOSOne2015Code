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

----------------------------------
Also included is the ApproxModel directory, which contains files for the approximate likelihood model MCMC and marginal likelihood calculations.

The files in the ApproxModel directory are:

InitOneStateTwoStateModelSelection.m
Initialisation script which first runs one-state and two-state simulations (both with and without measurement noise). It then runs the approximate model MCMC algorithm for one-state and two-state models, and finally calculates the marginal likelihood for the one-state and two-state models.


ChenOneStateNoiseMH.m
Implementation of Chen's method for approx one-state model.

Reference: Chen MH. Computing marginal likelihoods from a single MCMC output. Statistica Neerlandica. 2005 Jan;59(1):16–29. doi: 10.1111/j.1467-9574.2005.00276.x.


ChenTwoStateNoiseMH.m
Implementation of Chen's method for approx two-state model.


ChibOneStateNoiseMH.m
Implementation of Chib's method for approx one-state model. (This model doesn't require additional MCMC runs.)

Reference: Chib S, Jeliazkov I. Marginal likelihood from the Metropolis–Hastings output. Journal of the American Statistical Association. 2001;96(453):270–281. doi: 10.1198/016214501750332848.


ChibTwoStateNoiseMH.m
Implementation of Chib's method for approx two-state model. (Requires additional MCMC runs.)


LogLikelihoodOneStateMH.m
Calculate log likelihood for approx one-state model.


LogLikelihoodTwoStateMH.m
Calculate log likelihood for approx two-state model (forward algorithm). 


OneStateNoiseMH.m
MCMC algorithm (Metropolis-Hastings and Gibbs moves) for inference of approx one-state model.


TwoStateNoiseMH.m
MCMC algorithm for inference of approx two-state model.









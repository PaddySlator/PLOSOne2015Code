%priors
%BM
prior.OneState.D_max=10^6;
prior.OneState.D_min=0;
prior.OneState.noise_max=5000;
prior.OneState.a_D=0;
prior.OneState.b_D=0000;
prior.OneState.sigma_U=1000;
prior.OneState.sigma2=38.4761;
prior.OneState.S_sigma2=12.7066;
prior.OneState.FixedNoise=41.09;

%binding
prior.TwoState.a_1=1;
prior.TwoState.b_1=1;
prior.TwoState.a_2=1;
prior.TwoState.b_2=1;
prior.TwoState.D_max=10^6;
prior.TwoState.D_min=0;
prior.TwoState.a_D_1=0;
prior.TwoState.b_D_1=0000;
prior.TwoState.a_D_2=0;
prior.TwoState.b_D_2=0000;
prior.TwoState.sigma_U=1000;
prior.TwoState.mu_U=[0 0];
prior.TwoState.noise_max=2000;
prior.TwoState.sigma2=38.4761;
prior.TwoState.S_sigma2=12.7066;
prior.TwoState.FixedNoise=41.09;

%BM
alg_parameters.OneState.MCMC_steps=200;
alg_parameters.OneState.burn_in=100;
alg_parameters.OneState.bins=50;
alg_parameters.OneState.thin=1;
alg_parameters.OneState.SaveEachTraj=0;
alg_parameters.OneState.sigma2Estimate=0;

%binding
alg_parameters.TwoState.MCMC_steps=200;
alg_parameters.TwoState.burn_in=100;
alg_parameters.TwoState.thin=1;
alg_parameters.TwoState.bins=50;
alg_parameters.TwoState.MHVariance=1000;
alg_parameters.TwoState.SaveEachTraj=1;
alg_parameters.TwoState.sigma2Estimate=0;
alg_parameters.TwoState.LogLikelihood=0;
alg_parameters.TwoState.swap=1;
alg_parameters.TwoState.overdisp=1;
alg_parameters.TwoState.GelmanCutoff=1.2;

onchains.OneState=[1 0];
onchains.TwoState=[1 1 1 1 1 0];

%Simulate Trajectories
Treatments={'TwoDiffSim','OneDiffSim','TwoDiffNoiseSim','OneDiffNoiseSim'};
NTreatments=length(Treatments);
parametersTwoDiff=[10^5 10^4, 0.01 0.01 1 0.001];
parametersOneDiff=[10^5 10^5 0.01 0.01 1 0.001];

noise_parameters.variance=41.09;
noise_parameters.type='Gaussian';
N=4000;

NTraj=1;

for i=1:NTreatments
    switch Treatments{i}
        case 'TwoDiffSim'
            Traj.(Treatments{i})=...
                SimulateTwoStateDiffusion(parametersTwoDiff,[],N,prior.TwoState,[]);
        case 'OneDiffSim'
            Traj.(Treatments{i})=...
                SimulateTwoStateDiffusion(parametersOneDiff,[],N,prior.TwoState,[]);
        case 'TwoDiffNoiseSim'
            Traj.(Treatments{i})=...
                SimulateTwoStateDiffusion(parametersTwoDiff,noise_parameters,N,prior.TwoState,[]);
        case 'OneDiffNoiseSim'
            Traj.(Treatments{i})=...
                SimulateTwoStateDiffusion(parametersOneDiff,noise_parameters,N,prior.TwoState,[]);
    end
end


for i=1:NTreatments
    MCMCOutput.(Treatments{i}).OneState=OneStateNoiseMH(...
        Traj.(Treatments{i}),...
        alg_parameters.OneState,...
        prior.OneState,...
        onchains.OneState,0,0);
      
    MarginalLikelihood.(Treatments{i}).OneState=ChibOneStateNoiseMH(MCMCOutput.(Treatments{i}).OneState);
        
    
    MCMCOutput.(Treatments{i}).TwoState=TwoStateNoiseMH(...
        Traj.(Treatments{i}),...
        alg_parameters.TwoState,...
        prior.TwoState,...
        onchains.TwoState,0,0);
    
    MAPParameters=MCMCOutput.(Treatments{i}).TwoState.MAP;
    
    MAPLikelihood=LogLikelihoodTwoStateMH(Traj.(Treatments{i}),...
        MAPParameters);
    

    MarginalLikelihood.(Treatments{i}).TwoState=ChenTwoStateNoiseMH(...
        MCMCOutput.(Treatments{i}).TwoState,...
        MAPLikelihood); 
end


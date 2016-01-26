% Example use of SimulateTwoStateDiffusion (simulate two-state with measurement noise model) 
% and TwoStateDiffusionExactMCMC (infer two-state with measurement noise
% model).
% see Slator et al., PLOS ONE, 2015
% Paddy Slator, Warwick Systems Biology Centre

%vector of parameters
%[D_0
...D_1
...p_01
...p_10
...timestep option (1=equal timesteps, 2=gamma distributed)
...if option=1, constant timestep
...if option=2, gamma parameters]
%(e.g. parameters(5)=1, parameters(6)=0.001 is equal timesteps, length 0.001;
%parameters(5)=2, parameters(7)=1, parameters(8)=10 is Gamma(1,10) distributed timesteps)
parameters=[10^5 2*10^4 0.01 0.02 1 10^-3];

%set variance to 0 for no noise
noise_parameters.variance=20;

%number of timesteps 
N=2000;

%if using a non-zero prior on the first timepoint in the simulation
prior=[];
%option for simulation plots
ploton=[];

%simulate two-state diffusion with measurement noise trajectory
Traj=SimulateTwoStateDiffusion(parameters,noise_parameters,N,prior,ploton);

%prior parameters for MCMC inference
prior.a_D_0=0;
prior.b_D_0=0;
prior.a_D_1=0;
prior.b_D_1=0;
prior.a_0=1;
prior.b_0=1;
prior.a_1=1;
prior.b_1=1;
prior.D_max=10^6;
prior.D_min=0;
prior.mu_U=[0 0];
prior.sigma_U=1000^2;
%fixed Gaussian measurement noise for inference
%set to simulated value
prior.FixedNoise=noise_parameters.variance;

%parameters for MCMC inference
alg_parameters.MCMC_steps=1000;
alg_parameters.burn_in=500;
alg_parameters.thin=1;
%option for overdispersed starting points for parameters (e.g. Gelman
%statistic)
alg_parameters.overdisp=0;

%which MCMC chains to include, for debugging etc.
%[D_0 D_1 p_01 p_10 noise U z]
onchains=[1 1 1 1 1 1 1];

%Run MCMC inference
%MCMCOutput contains all hidden state (z, U) chains 
%(which can take up a lot of memory), 
%whereas MCMCOutputSummary doesn't
[MCMCOutput,MCMCOutputSummary]=TwoStateDiffusionExactMCMC(Traj,alg_parameters,prior,onchains,[],[]);

%%%plot MCMC output
%%parameter histograms and chains
ParameterLabels=MCMCOutput.ParameterLabels;
figure;
%plot diffusion coefficients 
for i=1:2    
    subplot(1,2,1);hold on;
    histogram(MCMCOutput.ParameterChains(alg_parameters.burn_in+1:end,i))
    plot(Traj.parameters(i),0,'o')
    ylabel('Diffusion coefficient')
    ylabel('Frequency')
    legend([strcat(ParameterLabels(i),' MCMC samples'),...
        strcat(ParameterLabels(i),' simulation value')])
    subplot(1,2,2);hold on;
    plot(MCMCOutput.ParameterChains(:,i))
    xlabel('MCMC step')
    ylabel('Diffusion coefficient')
end
legend(MCMCOutput.ParameterLabels(1:2))

subplot(1,2,1);hold on; 
%not nice, is there a way to append legend entries?
legend([strcat(ParameterLabels(1),' MCMC samples'),...
       strcat(ParameterLabels(1),' simulation value'),...
       strcat(ParameterLabels(2),' MCMC samples'),...
       strcat(ParameterLabels(2),' simulation value')])
   
   
figure;
%plot transition probabilities
for i=3:4    
    subplot(1,2,1);hold on;
    histogram(MCMCOutput.ParameterChains(alg_parameters.burn_in+1:end,i))
    plot(Traj.parameters(i),0,'o')
    xlabel('Transition probability')
    ylabel('Frequency')
    legend([strcat(ParameterLabels(i),' MCMC samples'),...
        strcat(ParameterLabels(i),' simulation value')])
    subplot(1,2,2);hold on;
    plot(MCMCOutput.ParameterChains(:,i))
    xlabel('MCMC step')
    ylabel('Transition probability')
end
legend(MCMCOutput.ParameterLabels(3:4))
 
subplot(1,2,1);hold on; 
%not nice, is there a way to append legend entries?
legend([strcat(ParameterLabels(3),' MCMC samples'),...
       strcat(ParameterLabels(3),' simulation value'),...
       strcat(ParameterLabels(4),' MCMC samples'),...
       strcat(ParameterLabels(4),' simulation value')])
    
%%hidden states
%z
figure; hold on;
plot(Traj.Y(1:end-1,3),Traj.z)
plot(Traj.Y(1:end-1,3),MCMCOutput.z_mean)
legend('simulated z','mean inferred z')
xlabel('Time (s)')
ylabel('Confinement probability')

%U
figure;
subplot(1,2,1);hold on;
plot(Traj.Y(:,3),Traj.U(:,1));
plot(Traj.Y(:,3),MCMCOutput.Ux_mean);
legend('simulated U_1','mean inferred U_1')
xlabel('Time (s)')
ylabel('U_1')

subplot(1,2,2);hold on;
plot(Traj.Y(:,3),Traj.U(:,2));
plot(Traj.Y(:,3),MCMCOutput.Uy_mean);
legend('simulated U_2','mean inferred U_2')
xlabel('Time (s)')
ylabel('U_2')

%%trajectory coloured by inferred z 
figure;hold on;axis off;
X=Traj.Y(1:end-1,1)';
Y=Traj.Y(1:end-1,2)';
Z=zeros(size(X));
col=MCMCOutput.z_mean;
caxis([0 1])
surface([X;X],[Y;Y],[Z;Z],[col;col],'facecol','no','edgecol','interp','linew',.5);
colorbar('Location','SouthOutside','Ticks',[0 1],'TickLabels',MCMCOutput.Parameters(1:2))











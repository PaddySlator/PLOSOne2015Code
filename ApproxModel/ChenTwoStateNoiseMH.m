function MarginalLikelihood=ChenTwoStateNoiseMH(MCMCOutput,LogLikelihood)
% Calculate Marginal likelkihood  for two-state diffusion with measurement
% noise (approximate likelihood) model using Chen's method
% 
% see Slator et al., PLOS ONE, 2015, and also
% Chen MH. Computing marginal likelihoods from a single MCMC output. 
% Statistica Neerlandica. 2005 Jan;59(1):16?29. doi: 10.1111/j.1467-9574.2005.00276.x.
%
% Paddy Slator, Warwick Systems Biology Centre

burn_in=MCMCOutput.alg_parameters.burn_in;
thin=MCMCOutput.alg_parameters.thin;
MCMC_steps=MCMCOutput.alg_parameters.MCMC_steps;
Traj=MCMCOutput.Traj;

N=length(Traj.Y)-1;
X=[Traj.Y(1:end-1,1), Traj.Y(1:end-1,2)];
dX=diff([Traj.Y(:,1),Traj.Y(:,2)]);
Dt=diff(Traj.Y(:,3));
R=diff(Traj.Y(:,1)).^2 + diff(Traj.Y(:,2)).^2;

D_1_MAP=MCMCOutput.MAP(1);
D_2_MAP=MCMCOutput.MAP(2);
p_12_MAP=MCMCOutput.MAP(3);
p_21_MAP=MCMCOutput.MAP(4);
noise_MAP=MCMCOutput.MAP(5);

a_1=MCMCOutput.prior.a_1;
b_1=MCMCOutput.prior.b_1;
a_2=MCMCOutput.prior.a_2;
b_2=MCMCOutput.prior.b_2;
D_max=MCMCOutput.prior.D_max;

NPostSamples=MCMC_steps/thin-burn_in/thin;

%samples from the posterior distriubtion
D_1Post=MCMCOutput.ParameterPosteriorSamples(:,1);
D_2Post=MCMCOutput.ParameterPosteriorSamples(:,2);
p_12Post=MCMCOutput.ParameterPosteriorSamples(:,3);
p_21Post=MCMCOutput.ParameterPosteriorSamples(:,4);
noisePost=MCMCOutput.ParameterPosteriorSamples(:,5);
zPost=MCMCOutput.zPosteriorSamples;

ChenTerm=zeros(NPostSamples,1);

for i=1:NPostSamples
    
    z=zPost(i,:);
    n_12 = 0;
    n_11 = 0;
    n_21 = 0;
    n_22 = 0;
    
    for j=1:N-1
        if z(j) == 1
            if z(j+1) == 2
                n_12 = n_12 + 1;
            else
                n_11 = n_11 + 1;
            end
        else
            if z(j+1) == 1
                n_21 = n_21 + 1;
            else
                n_22 = n_22 + 1;
            end
        end
    end
    
    ChenTerm(i)=...
        ...numerator
        ...\pi(X|theta^*,z_i)
        ...z==1 part
        -sum(log(2*pi*(2*D_1_MAP*Dt(zPost(i,:)==1) + 2*noise_MAP)))...
        -sum(sum(dX(zPost(i,:)==1,:).^2,2)./(4*noise_MAP+4*D_1_MAP*Dt(zPost(i,:)==1)))...
        ...z==2 part
        -sum(log(2*pi*(2*D_2_MAP*Dt(zPost(i,:)==2) + 2*noise_MAP)))...
        -sum(sum(dX(zPost(i,:)==2,:).^2,2)./(4*noise_MAP+4*D_2_MAP*Dt(zPost(i,:)==2)))...
        ...\pi(z_i|theta^*)
        +log(betapdf(p_12_MAP,n_12+1,n_11+1))...
        +log(betapdf(p_21_MAP,n_21+1,n_22+1))...
        ...denominator
            -(...
        ...\pi(X|theta_i,z_i)
        ...z==1 part
        -sum(log(2*pi*(2*D_1Post(i)*Dt(zPost(i,:)==1) + 2*noisePost(i))))...
        -sum(sum(dX(zPost(i,:)==1,:).^2,2)./(4*noisePost(i)+4*D_1Post(i)*Dt(zPost(i,:)==1)))...
        ...z==2 part
        -sum(log(2*pi*(2*D_2Post(i)*Dt(zPost(i,:)==2) + 2*noisePost(i))))...
        -sum(sum(dX(zPost(i,:)==2,:).^2,2)./(4*noisePost(i)+4*D_2Post(i)*Dt(zPost(i,:)==2)))...
        ...\pi(z_i|theta_i)
        +log(betapdf(p_12Post(i),n_12+1,n_11+1))...
        +log(betapdf(p_21Post(i),n_21+1,n_22+1))...
            );
               
end

A=max(ChenTerm);

MarginalCorrection=-log(1/NPostSamples) ...
    +A + log(sum(exp(ChenTerm-A)));


if isempty(LogLikelihood)
    LogLikelihood=LogLikelihoodBindingMH(Traj,MCMCOutput.MAP);
end

MarginalLikelihood=LogLikelihood - MarginalCorrection;


end
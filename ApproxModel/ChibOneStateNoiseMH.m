function MarginalLikelihood = ChibOneStateNoiseMH(MCMCOutput)

burn_in=MCMCOutput.alg_parameters.burn_in;
thin=MCMCOutput.alg_parameters.thin;
MCMC_steps=MCMCOutput.alg_parameters.MCMC_steps;
Traj=MCMCOutput.Traj;

N=length(Traj.Y)-1;
X=[Traj.Y(1:end-1,1), Traj.Y(1:end-1,2)];
dX=diff([Traj.Y(:,1),Traj.Y(:,2)]);
Dt=diff(Traj.Y(:,3));
R=diff(Traj.Y(:,1)).^2 + diff(Traj.Y(:,2)).^2;

D_MAP=MCMCOutput.MAP(1);
noise_MAP=MCMCOutput.MAP(2);

D_max=MCMCOutput.prior.D_max;

NPostSamples=(MCMC_steps/thin-burn_in/thin);

MH_D_SD=MCMCOutput.MH_D_SD;
MHNoise_SD=MCMCOutput.MHNoise_SD;

%samples from posterior
DPost=MCMCOutput.ParameterPosteriorSamples(:,1);
noisePost=MCMCOutput.ParameterPosteriorSamples(:,2);

%samples from jumping distribution q(D*,sigma2* -> D,sigma2)
DJump=zeros(NPostSamples,1);
noiseJump=zeros(NPostSamples,1);

LogAlphaPost=zeros(NPostSamples,1);

Logq=zeros(NPostSamples,1);

LogAlphaJump=zeros(NPostSamples,1);

LogLikelihoodMAP = -sum(log(2*pi*(2*D_MAP*Dt + 2*noise_MAP)))...
    -sum(sum(dX.^2,2)./(2*2*noise_MAP+4*D_MAP*Dt));
    
sigma2_Michalet=MCMCOutput.sigma2_Michalet;

for i=1:NPostSamples
    
    %log \alpha(DPost^(m),D*|X)
    LogLikelihoodPost = -sum(log(2*pi*(2*DPost(i)*Dt + 2*noisePost(i))))...
        -sum(sum(dX.^2,2)./(2*2*noisePost(i)+4*DPost(i)*Dt));
        
    LogAlphaPost(i)=min(log(1),...
        LogLikelihoodMAP-LogLikelihoodPost);
                
    %log q(DPost^(m),D*)
    Logq(i)=log(normpdf(DPost(i),D_MAP,MH_D_SD));
        
        
    %sample DJumps's from jump distribution q(D*,D) (Chib,2001,JASA) (choosing D*=D_MAP)
    DJump(i)= D_MAP + randn*sqrt(MH_D_SD);
    %sample noiseJump's from jump dist
    noiseJump(i)= noise_MAP + randn*sqrt(MHNoise_SD);
    
    %log \alpha(D*, DJump^(j)|X)
    LogLikelihoodJump = -sum(log(2*pi*(2*DJump(i)*Dt + 2*noiseJump(i))))...
        -sum(sum(dX.^2,2)./(2*2*noiseJump(i)+4*DJump(i)*Dt));
    
    LogAlphaJump(i)=min(log(1),...
        LogLikelihoodJump-LogLikelihoodMAP);
end


MarginalCorrection=-log(NPostSamples)...
    +log(sum(exp(LogAlphaPost+Logq)))...
    +log(NPostSamples)...
    -log(sum(exp(LogAlphaJump)));

LogLikelihood=LogLikelihoodOneStateMH(Traj,[D_MAP noise_MAP])
LogPrior=-log(D_max)-log(normpdf(noise_MAP,sigma2_Michalet,MH_D_SD))

MarginalLikelihood=LogLikelihood...
    +LogPrior...
    -MarginalCorrection;



MarginalCorrection



end
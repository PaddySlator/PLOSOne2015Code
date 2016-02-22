function MCMCOutput = OneStateNoiseMH(Traj,alg_parameters,prior,onchains,initial_values,figures)
% MCMC algorithm (Metropolis-Hastings and Gibbs moves) for inference of approx one-state model.
% see Slator et al., PLOS ONE, 2015
% Paddy Slator, Warwick Systems Biology Centre



MCMC_steps=alg_parameters.MCMC_steps;
burn_in=alg_parameters.burn_in;
bins=alg_parameters.bins;
thin=alg_parameters.thin;

D_max=prior.D_max;
D_min=prior.D_min;

N=length(Traj.Y)-1;
X=[Traj.Y(1:end-1,1), Traj.Y(1:end-1,2)];
dX=diff([Traj.Y(:,1),Traj.Y(:,2)]);
Dt=diff(Traj.Y(:,3));

DChain=zeros(MCMC_steps/thin,1);
noiseChain=zeros(MCMC_steps/thin,1);
LogLikelihood=zeros(MCMC_steps/thin,1);


MH_D_SD=30;
MHNoise_SD=10;

onD=onchains(1);
onnoise=onchains(2);

%MSDAnalysis=MichaletAnalysis(O,0);

if alg_parameters.sigma2Estimate
    prior.sigma2=MSDAnalysis.sigma2;
    prior.S_sigma2=MSDAnalysis.S_sigma2;
end 

%D_min=max(0,-MSDAnalysis.sigma2/mean(Dt))
prior.D_min=D_min;
sigma2_Michalet=prior.sigma2;
S_sigma2_Michalet=prior.S_sigma2;

if onD
    if isfield(initial_values,'D')
        D=initial_values.D;
    else
        D=D_min+rand*(D_max-D_min);
        %D=MSDAnalysis.D;
        while D<D_min || D>D_max
            D=D_min+rand*(D_max-D_min);
        end
    end
else
    D=Traj.parameters(1);
end

N_SD=2;
noise_min=sigma2_Michalet-N_SD*S_sigma2_Michalet
noise_max=sigma2_Michalet+N_SD*S_sigma2_Michalet

if onnoise
    if isfield(initial_values,'noise')
        noise=initial_values.noise;
    else        
        noise=sigma2_Michalet+randn*MHNoise_SD;
        while noise < noise_min || noise > noise_max
            noise=sigma2_Michalet+randn*MHNoise_SD;
        end
    end
else
    noise=prior.FixedNoise;
end
DChain(1)=D;
noiseChain(1)=noise;


TuneFreq=100;
TuneProportion=0.1;

Moves=0;
MovesD=0;
MovesNoise=0;
TotalMoves=0;


Li = -sum(log(2*pi*(2*D*Dt + 2*noise)))...
    -sum(sum(dX.^2,2)./(2*2*noise+4*D*Dt));

% Li = -sum(log(2*pi*(2*D*Dt + 2*noise - (2/3)*D*0.001 )))...
%      -sum(sum(dX.^2,2)./(2*2*noise - 2*(2/3)*D*0.001 +4*D*Dt));

LogLikelihood(1)=Li;

if onnoise
    Prior=log(normpdf(noise,sigma2_Michalet,MHNoise_SD));
else 
    Prior=0;
end

SingleMHUpdate=1;
JointMHUpdate=0;

if JointMHUpdate
   MH_D_SD=1;
   MHNoise_SD=1; 
end


for i=2:MCMC_steps
    
    if JointMHUpdate
        DProp=D+randn*MH_D_SD;
        
        noiseProp=noise+randn*MHNoise_SD;
        
        LiProp = -sum(log(2*pi*(2*DProp*Dt + 2*noiseProp)))...
            -sum(sum(dX.^2,2)./(4*noiseProp+4*DProp*Dt));
        
        PriorProp=log(normpdf(noiseProp,sigma2_Michalet,S_sigma2_Michalet));
        
        if  (DProp > D_min) && (DProp < D_max) && (noiseProp > noise_min) && (noiseProp < noise_max)
            if (LiProp + PriorProp) > (Li + Prior)
                D=DProp;
                noise=noiseProp;
                Li=LiProp;
                Prior=PriorProp;
                Moves=Moves+1;
                TotalMoves=TotalMoves+1;
            elseif log(rand) < (LiProp+PriorProp-Li-Prior)
                D=DProp;
                noise=noiseProp;
                Li=LiProp;
                Prior=PriorProp;
                Moves=Moves+1;
                TotalMoves=TotalMoves+1;
            end
        end
        
        DChain(i)=D;
        noiseChain(i)=noise;
        
        LogLikelihood(i)=Li;
                
        %     %Tune Metropolis-Hastings jumping variance
        if i < burn_in
            if rem(i,TuneFreq)==0
                AcceptanceRate=Moves/TuneFreq;
                if AcceptanceRate < 0.2
                    MH_D_SD=MH_D_SD*(1-TuneProportion);
                end
                if AcceptanceRate > 0.3
                    MH_D_SD=MH_D_SD*(1+TuneProportion);
                end
                Moves=0;
            end
        end
    end
 
    if SingleMHUpdate
        if onD
            DProp=D+randn*MH_D_SD;
%             while DProp < D_min || DProp > D_max
%                 DProp=D+randn*MH_D_SD
%             end


%                 LiProp = -sum(log(2*pi*(2*DProp*Dt + 2*noise - (2/3)*DProp*0.001 )))...
%                     -sum(sum(dX.^2,2)./(2*2*noise - 2*(2/3)*DProp*0.001 +4*DProp*Dt));

            LiProp = -sum(log(2*pi*(2*DProp*Dt + 2*noise)))...
                -sum(sum(dX.^2,2)./(4*noise+4*DProp*Dt));

            
            LiProp = LogLikelihoodOneStateMH(Traj,[DProp noise]);
            
            D_min=max(0,-noise/mean(Dt));
                        
            
            if onnoise
                PriorProp=log(normpdf(noise,sigma2_Michalet,S_sigma2_Michalet));
            else
                PriorProp=0;
            end
            
            if (DProp > D_min) && (DProp < D_max)
                if isreal(LiProp)
                    if (LiProp + PriorProp) > (Li + Prior)
                        D=DProp;
                        Li=LiProp;
                        Prior=PriorProp;
                        MovesD=MovesD+1;
                        TotalMoves=TotalMoves+1;
                    elseif log(rand) < (LiProp+PriorProp-Li-Prior)
                        D=DProp;
                        Li=LiProp;
                        Prior=PriorProp;
                        MovesD=MovesD+1;
                        TotalMoves=TotalMoves+1;
                    end
                end
            end
        end
    
        if onnoise
            noiseProp=noise+randn*MHNoise_SD;
                        
                        %LiProp = -sum(log(2*pi*(2*D*Dt + 2*noiseProp - (2/3)*D*0.001 )))...
            %         -sum(sum(dX.^2,2)./(4*noiseProp - (4/3)*D*0.001 +4*D*Dt));
            
            LiProp = -sum(log(2*pi*(2*D*Dt + 2*noiseProp)))...
                -sum(sum(dX.^2,2)./(4*noiseProp+4*D*Dt));
            
            LiProp = LogLikelihoodBrownianMH(Traj,[D noiseProp]);
            
            PriorProp=log(normpdf(noiseProp,sigma2_Michalet,S_sigma2_Michalet));
            
            if (noiseProp > noise_min) && (noiseProp < noise_max)
                if (LiProp + PriorProp) > (Li + Prior)
                    noise=noiseProp;
                    Li=LiProp;
                    Prior=PriorProp;
                    MovesNoise=MovesNoise+1;
                    TotalMoves=TotalMoves+1;
                elseif log(rand) < (LiProp+PriorProp-Li-Prior)
                    noise=noiseProp;
                    Li=LiProp;
                    Prior=PriorProp;
                    MovesNoise=MovesNoise+1;
                    TotalMoves=TotalMoves+1;
                end
            end
            
        end
        
        if rem(i,thin)==0
            DChain(i/thin)=D;
            noiseChain(i/thin)=noise;
            
            LogLikelihood(i/thin)=Li;
        end
        
        
        %Tune Metropolis-Hastings jumping variance
        if i < burn_in
            if rem(i,TuneFreq)==0
                AcceptanceRateD=MovesD/TuneFreq;
                AcceptanceRatenoise=MovesNoise/TuneFreq;
                
                if AcceptanceRateD < 0.2
                    MH_D_SD=MH_D_SD*(1-TuneProportion);
                end
                if AcceptanceRateD > 0.3
                    MH_D_SD=MH_D_SD*(1+TuneProportion);
                end
                if AcceptanceRatenoise < 0.2
                    MHNoise_SD=MHNoise_SD*(1-TuneProportion);
                end
                if AcceptanceRatenoise > 0.3
                    MHNoise_SD=MHNoise_SD*(1+TuneProportion);
                end
%                 AcceptanceRateD
%                 MH_D_SD
%                 AcceptanceRatenoise
%                 MHNoise_SD
                MovesD=0;
                MovesNoise=0;
            end
        end
    end
        

end


MCMCOutput.sigma2_Michalet=sigma2_Michalet;
MCMCOutput.MH_D_SD=MH_D_SD;
MCMCOutput.MHNoise_SD=MHNoise_SD;
MCMCOutput.TotalMoves=TotalMoves;

MCMCOutput.DChain=DChain;
MCMCOutput.noiseChain=noiseChain;
MCMCOutput.LogLikelihood=LogLikelihood;

MCMCOutput.MAP(1)=mean(DChain(burn_in/thin+1:end));
MCMCOutput.MAP(2)=mean(noiseChain(burn_in/thin+1:end));

MAPLikelihood=LogLikelihoodOneStateMH(Traj,[MCMCOutput.MAP noise]);
MCMCOutput.MAPLikelihood=MAPLikelihood;
%BIC
NParameters=1;
BIC=-2*MAPLikelihood+NParameters*log(N);
MCMCOutput.BIC=BIC;
%AIC
AIC=-2*MAPLikelihood-2*NParameters;
MCMCOutput.AIC=AIC;

MCMCOutput.Parameters={'D','sigma2'};

MCMCOutput.alg_parameters=alg_parameters;
MCMCOutput.prior=prior;

MCMCOutput.FixedNoise=noise;

MCMCOutput.Traj=Traj;

MCMCOutput.ParameterPosteriorSamples=[DChain(burn_in/thin+1:end) noiseChain(burn_in/thin+1:end)];


MarginalLikelihood=ChibOneStateNoiseMH(MCMCOutput);

MCMCOutput.MarginalLikelihood=MarginalLikelihood;



end
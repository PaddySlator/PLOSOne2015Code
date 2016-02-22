function [MCMCOutput,MCMCOutputSummary] = TwoStateNoiseMH(Traj,alg_parameters,prior,onchains,initial_values,figures)
% MCMC algorithm for inference of approx two-state model.
% see Slator et al., PLOS ONE, 2015
% Paddy Slator, Warwick Systems Biology Centre


MCMC_steps=alg_parameters.MCMC_steps;
burn_in=alg_parameters.burn_in;
bins=alg_parameters.bins;
thin=alg_parameters.thin;

D_max=prior.D_max;
D_min=prior.D_min;

a_1=prior.a_1;
b_1=prior.b_1;
a_2=prior.a_2;
b_2=prior.b_2;

N=length(Traj.Y)-1;
X=[Traj.Y(1:end-1,1), Traj.Y(1:end-1,2)];
dX=diff([Traj.Y(:,1),Traj.Y(:,2)]);
Dt=diff(Traj.Y(:,3));

D_1_chain=zeros(MCMC_steps/thin,1);
D_2_chain=zeros(MCMC_steps/thin,1);
noise_chain=zeros(MCMC_steps/thin,1);
p_12_chain=zeros(MCMC_steps/thin,1);
p_21_chain=zeros(MCMC_steps/thin,1);

z_chain=zeros(MCMC_steps/thin,N);

LogLikelihood=zeros(MCMC_steps/thin,1);

onD_1=onchains(1);
onD_2=onchains(2);
onp_12=onchains(3);
onp_21=onchains(4);
onz=onchains(5);
onnoise=onchains(6);



D_min=prior.D_min;

NoiseMean=prior.sigma2;
MHNoise_SD=prior.S_sigma2;

if isfield(initial_values,'TunedMHSD')
    MHD_1_SD=initial_values.TunedMHSD.D_1
    MHD_2_SD=initial_values.TunedMHSD.D_2
else
    MHD_1_SD=5000
    MHD_2_SD=2000
end

N_SD=2;
noise_min=NoiseMean-N_SD*MHNoise_SD;
noise_max=NoiseMean+N_SD*MHNoise_SD;

if onnoise
    if isfield(initial_values,'noise')
        noise=initial_values.noise;
    else
        noise=NoiseMean+randn*MHNoise_SD;
        while noise<noise_min || noise>noise_max
            noise=NoiseMean+randn*MHNoise_SD;
        end
    end
else
    noise=prior.FixedNoise;  
end



if onp_12
    if isfield(initial_values,'p_12')
        p_12=initial_values.p_12;
    else
        p_12=betarnd(a_1,b_1);
        if alg_parameters.overdisp
            p_12=rand;
        end
    end
else
    p_12=Traj.parameters(3);
end
if onp_21
    if isfield(initial_values,'p_21')
        p_21=initial_values.p_21;
    else
        p_21=betarnd(a_2,b_2);
        if alg_parameters.overdisp
            p_21=rand;
        end
    end
else
    p_21=Traj.parameters(4);
end

InitialLogLikelihood=sqrt(-1);
while ~isreal(InitialLogLikelihood)
    if onD_1
        if isfield(initial_values,'D_1')
            D_1=initial_values.D_1;
        else
            D_1=D_min+rand*(D_max-D_min);
            D_1=D_min+betarnd(0.1,0.1)*rand*D_max;
            
            while D_1<D_min || D_1>D_max
                D_1=D_min+rand*(D_max-D_min);
                if alg_parameters.overdisp
                   D_1=D_min+betarnd(0.1,0.1)*D_max;
                end
            end
        end
    else
        D_1=Traj.parameters(1);
    end
    if onD_2
        if isfield(initial_values,'D_2')
            D_2=initial_values.D_2;
        else
            D_2=D_min+rand*(D_max-D_min);
            D_2=D_min+betarnd(0.1,0.1)*rand*(D_max-D_min)
            while D_2<D_min || D_2>D_max
                D_2=D_min+rand*(D_max-D_min);
                if alg_parameters.overdisp
                    D_2=D_min+betarnd(0.1,0.1)*(D_max-D_min);
                end
            end
        end
    else
        D_2=Traj.parameters(2);
    end
    D_min
    D_max
    D_1 
    D_2 
    p_12 
    p_21 
    noise
    InitialLogLikelihood=LogLikelihoodTwoStateMH(Traj,[D_1 D_2 p_12 p_21 noise])
    if alg_parameters.LogLikelihood
        LogLikelihood(1)=InitialLogLikelihood;
    end
end


SingleMHUpdate=1;
JointMHUpdate=0;
onjoint_D_1_D_2=0;


if JointMHUpdate
    if onD_1 && onD_2
        onjoint_D_1_D_2=1;
        onD_1=0;
        onD_2=0;
        onnoise=0;
    else
        onjoint_D_1_D_2=0;
    end
end


% onjoint_D_1_D_2=0;
% onD_1=1;
% onD_2=1;

if onz
    if isfield(initial_values,'z')
        z=initial_values.z;
    else
        pi_1=p_21/(p_21 + p_12);
        pi_2= 1 - pi_1;
        
        z=ones(N,1);
        
        if rand < pi_1
            z(1) = 1;
        else
            z(1) = 2;
        end
        
        change12=[];
        change21=[];
        for i=2:N
            if z(i-1) == 1 && rand < p_12
                % binding event
                z(i) = 2;
                change12=[change12,i];
            elseif z(i-1) == 2 && rand < p_21
                % unbinding event
                z(i) = 1;
                change21=[change21,i];
            else
                % no event
                z(i)=z(i-1);
            end
        end
    end
else
    z=Traj.z;
end




%noise=O.noise_parameters.variance;
i=1;
D_1_chain(i)=D_1;
D_2_chain(i)=D_2;
noise_chain(i)=noise;
p_12_chain(i)=p_12;
p_21_chain(i)=p_21;
z_chain(i,:)=z;




MH_SD=1000;


%"likelihoods" for MH
LikelihoodD_1= -sum(log(2*pi*(2*D_1*Dt(z==1) + 2*noise)))...
    -sum(sum(dX(z==1,:).^2,2)./(4*noise+4*D_1*Dt(z==1)));

LikelihoodD_2= -sum(log(2*pi*(2*D_2*Dt(z==2) + 2*noise)))...
    -sum(sum(dX(z==2,:).^2,2)./(4*noise+4*D_2*Dt(z==2)));


%LogLikelihood(i)=LikelihoodD_2;

TuneFreq=50;
TuneProportion=0.1;

MovesD_1=0;
TotalMovesD_1=0;

MovesD_2=0;
TotalMovesD_2=0;

Movesnoise=0;
TotalMovesnoise=0;

Moves=0;
TotalMoves=0;



for i=2:MCMC_steps
    
    if onD_1
        D_1_Prop=D_1+randn*MHD_1_SD;
        %         while D_1_Prop < D_min || D_1_Prop > D_max
        %             D_1_Prop=D_1+randn*MHD_1_SD;
        %         end
        
        %D_min=max(0,-noise/mean(Dt));

        LikelihoodD_1= -sum(log(2*pi*(2*D_1*Dt(z==1) + 2*noise)))...
            -sum(sum(dX(z==1,:).^2,2)./(4*noise+4*D_1*Dt(z==1)));
        
        LikelihoodD_1_Prop= -sum(log(2*pi*(2*D_1_Prop*Dt(z==1) + 2*noise)))...
            -sum(sum(dX(z==1,:).^2,2)./(4*noise+4*D_1_Prop*Dt(z==1)));
                    
        if (D_1_Prop > D_min) && (D_1_Prop<D_max)
            if isreal(LikelihoodD_1_Prop)
                if LikelihoodD_1_Prop > LikelihoodD_1
                    D_1=D_1_Prop;
                    LikelihoodD_1=LikelihoodD_1_Prop;
                    MovesD_1=MovesD_1+1;
                    TotalMovesD_1=TotalMovesD_1+1;
                elseif log(rand) < (LikelihoodD_1_Prop-LikelihoodD_1)
                    D_1=D_1_Prop;
                    LikelihoodD_1=LikelihoodD_1_Prop;
                    MovesD_1=MovesD_1+1;
                    TotalMovesD_1=TotalMovesD_1+1;
                end
            end
        end
    end
    
    if onD_2
        D_2_Prop=D_2+randn*MHD_2_SD;
        %         while D_2_Prop < D_min || D_2_Prop > D_max
        %             D_2_Prop=D_2+randn*MHD_2_SD;
        %         end
        
%         D_min=max(0,noise/mean(Dt));
                        
        LikelihoodD_2= -sum(log(2*pi*(2*D_2*Dt(z==2) + 2*noise)))...
            -sum(sum(dX(z==2,:).^2,2)./(4*noise+4*D_2*Dt(z==2)));
        
        LikelihoodD_2_Prop= -sum(log(2*pi*(2*D_2_Prop*Dt(z==2) + 2*noise)))...
            -sum(sum(dX(z==2,:).^2,2)./(4*noise+4*D_2_Prop*Dt(z==2)));
                
        if (D_2_Prop > D_min) && (D_2_Prop < D_max)
            if isreal(LikelihoodD_2_Prop)
                if LikelihoodD_2_Prop > LikelihoodD_2
                    D_2=D_2_Prop;
                    LikelihoodD_2=LikelihoodD_2_Prop;
                    MovesD_2=MovesD_2+1;
                    TotalMovesD_2=TotalMovesD_2+1;
                elseif log(rand) < (LikelihoodD_2_Prop-LikelihoodD_2)
                    D_2=D_2_Prop;
                    LikelihoodD_2=LikelihoodD_2_Prop;
                    MovesD_2=MovesD_2+1;
                    TotalMovesD_2=TotalMovesD_2+1;
                end
            end
        end
    end
    
    if onnoise
       noiseProp=noise+randn*MHNoise_SD; 
       
       Likelihoodnoise=-sum(log(2*pi*(2*D_1*Dt(z==1) + 2*noise)))...
            -sum(sum(dX(z==1,:).^2,2)./(4*noise+4*D_1*Dt(z==1)))...
            -sum(log(2*pi*(2*D_2*Dt(z==2) + 2*noise)))...
            -sum(sum(dX(z==2,:).^2,2)./(4*noise+4*D_2*Dt(z==2)));
       
       Prior=log(normpdf(noise,NoiseMean,prior.S_sigma2));
       
       LikelihoodnoiseProp=-sum(log(2*pi*(2*D_1*Dt(z==1) + 2*noiseProp)))...
            -sum(sum(dX(z==1,:).^2,2)./(4*noiseProp+4*D_1*Dt(z==1)))...
            -sum(log(2*pi*(2*D_2*Dt(z==2) + 2*noiseProp)))...
            -sum(sum(dX(z==2,:).^2,2)./(4*noiseProp+4*D_2*Dt(z==2)));
       
       PriorProp=log(normpdf(noiseProp,NoiseMean,prior.S_sigma2));
       
       if (noiseProp > noise_min) && (noiseProp < noise_max)
           if isreal(LikelihoodnoiseProp)
               if (LikelihoodnoiseProp + PriorProp) > (Likelihoodnoise + Prior)
                   noise=noiseProp;
                   Movesnoise=Movesnoise+1;
                   TotalMovesnoise=TotalMovesnoise+1;
               elseif log(rand) < (LikelihoodnoiseProp + PriorProp - Likelihoodnoise - Prior)
                   noise=noiseProp;
                   Movesnoise=Movesnoise+1;
                   TotalMovesnoise=TotalMovesnoise+1;
               end
           end
       end
                   
    end
    

   
    if onp_12 || onp_21
        %Update p01 and p10
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
    end
    
    if onp_12
        p_12=betarnd(a_1+n_12,b_1+n_11);
    end
    
    if onp_21
        p_21=betarnd(a_2+n_21,b_2+n_22);
    end
    
    
    if onz
        P_t = [ 1 - p_12, p_12
            p_21, 1 - p_21];
        
        j=1;
        P_1=(p_21/(p_12+p_21))*...
            exp(-(sum(dX(j,:).^2))/(4*(D_1*Dt(j)+noise)))/((4*pi*(D_1*Dt(j)+noise)))...
            *P_t(1,z(j+1));
        
    
        P_2=(p_12/(p_12+p_21))*...
            exp(-(sum(dX(j,:).^2))/(4*(D_2*Dt(j)+noise)))/((4*pi*(D_2*Dt(j)+noise)))...
            *P_t(2,z(j+1));
        
        if rand < P_1/(P_1+P_2)
            z(j)=1;
        else
            z(j)=2;
        end
        
               
        for j=2:N-1            
            P_1=P_t(z(j-1),1)...
                *exp(-(sum(dX(j,:).^2))/(4*(D_1*Dt(j)+noise)))/((4*pi*(D_1*Dt(j)+noise)))...
                *P_t(1,z(j+1));            
            
            P_2=P_t(z(j-1),2)...
                *exp(-(sum(dX(j,:).^2))/(4*(D_2*Dt(j)+noise)))/((4*pi*(D_2*Dt(j)+noise)))...
                *P_t(2,z(j+1));
            
            if rand < P_1/(P_1+P_2)
                z(j)=1;
            else
                z(j)=2;
            end
        end
        
        j=N;
        P_1=P_t(z(j-1),1)...
            *exp(-(sum(dX(j,:).^2))/(4*(D_1*Dt(j)+noise)))/((4*pi*(D_1*Dt(j)+noise)));
        P_2=P_t(z(j-1),2)...
            *exp(-(sum(dX(j,:).^2))/(4*(D_2*Dt(j)+noise)))/((4*pi*(D_2*Dt(j)+noise)));

        
        if rand < P_1/(P_1+P_2)
            z(j)=1;
        else
            z(j)=2;
        end                
    end
     
    %Tune Metropolis-Hastings jumping variance
    if i < burn_in
        if rem(i,TuneFreq)==0
            
            if onjoint_D_1_D_2
                AcceptanceRate=Moves/TuneFreq
                if AcceptanceRate < 0.2
                    MH_SD=MH_SD*(1-TuneProportion)
                    if AcceptanceRateD_1 < 0.02
                        MH_SD=MH_SD*TuneProportion
                    end
                end
                if AcceptanceRate > 0.3
                    MH_SD=MH_SD*(1+TuneProportion)
                end                
                Moves=0;
            end
            
            
            if onD_1
                AcceptanceRateD_1=MovesD_1/TuneFreq;
                if AcceptanceRateD_1 < 0.2
                    MHD_1_SD=MHD_1_SD*(1-TuneProportion);                    
                end
                if AcceptanceRateD_1 > 0.3
                    MHD_1_SD=MHD_1_SD*(1+TuneProportion);
                end
                MovesD_1=0;
            end
            if onD_2
                AcceptanceRateD_2=MovesD_2/TuneFreq;
                if AcceptanceRateD_2 < 0.2
                    MHD_2_SD=MHD_2_SD*(1-TuneProportion);                    
                end
                if AcceptanceRateD_2 > 0.3
                    MHD_2_SD=MHD_2_SD*(1+TuneProportion);
                end
                MovesD_2=0;
            end            
            if onnoise
                AcceptanceRatenoise=Movesnoise/TuneFreq;
                if AcceptanceRatenoise < 0.2
                    MHNoise_SD=MHNoise_SD*(1-TuneProportion);                    
                end
                if AcceptanceRatenoise > 0.3
                    MHNoise_SD=MHNoise_SD*(1+TuneProportion);
                end
                Movesnoise=0;
            end
        end
    end
    
    if rem(i,thin)==0
        D_1_chain(i/thin)=D_1;
        D_2_chain(i/thin)=D_2;
        noise_chain(i/thin)=noise;
        p_12_chain(i/thin)=p_12;
        p_21_chain(i/thin)=p_21;
        
        z_chain(i/thin,:)=z;
        
        if alg_parameters.LogLikelihood
            LogLikelihood(i/thin)=LogLikelihoodTwoStateMH(Traj,[D_1 D_2 p_12 p_21 noise]);
        end
    end
       
end

AcceptanceRate=Moves/(MCMC_steps-burn_in);
MCMCOutput.AcceptanceRate=AcceptanceRate;

if alg_parameters.swap
    %swap chains to make D_2 the bound state
    if mean(D_1_chain(burn_in/thin+1:end)) < mean(D_2_chain(burn_in/thin+1:end))
        [D_2_chain,D_1_chain]=deal(D_1_chain,D_2_chain);
        [p_21_chain,p_12_chain]=deal(p_12_chain,p_21_chain);
        
        [MHD_2_SD,MHD_1_SD]=deal(MHD_1_SD,MHD_2_SD);
            
        for i=2:MCMC_steps/thin
            for j=1:N
                if z_chain(i,j)==1
                    z_chain(i,j)=2;
                else
                    z_chain(i,j)=1;
                end
            end
        end
    end
end

MCMCOutput.D_1_chain=D_1_chain;
MCMCOutput.D_2_chain=D_2_chain;
MCMCOutput.noise_chain=noise_chain;
MCMCOutput.p_12_chain=p_12_chain;
MCMCOutput.p_21_chain=p_21_chain;
MCMCOutput.z_chain=z_chain;

z_mean=mean(MCMCOutput.z_chain(burn_in/thin+1:end,:));


MCMCOutput.z_mean=z_mean;
MCMCOutput.z_post=z_mean-1;

MCMCOutput.LogLikelihood=LogLikelihood;


MCMCOutput.MAP=[mean(D_1_chain(burn_in/thin+1:end))...
    mean(D_2_chain(burn_in/thin+1:end))...
    mean(p_12_chain(burn_in/thin+1:end))...
    mean(p_21_chain(burn_in/thin+1:end))...
    mean(noise_chain(burn_in/thin+1:end))];

MCMCOutput.FinalState.D_1=D_1_chain(end);
MCMCOutput.FinalState.D_2=D_2_chain(end);
MCMCOutput.FinalState.p_12=p_12_chain(end);
MCMCOutput.FinalState.p_21=p_21_chain(end);
MCMCOutput.FinalState.noise=noise_chain(end);
MCMCOutput.FinalState.z=z_chain(end,:);


MCMCOutput.alg_parameters=alg_parameters;
MCMCOutput.prior=prior;

MCMCOutput.FixedNoise=noise;

MCMCOutput.Traj=Traj;

MCMCOutput.TotalMovesD_1=TotalMovesD_1;
MCMCOutput.TotalMovesD_2=TotalMovesD_2;

MCMCOutput.ParameterPosteriorSamples=[D_1_chain(burn_in/thin+1:end)...
    D_2_chain(burn_in/thin+1:end)...
    p_12_chain(burn_in/thin+1:end)...
    p_21_chain(burn_in/thin+1:end)...
    noise_chain(burn_in/thin+1:end)];

MCMCOutput.zPosteriorSamples=z_chain(burn_in/thin+1:end,:);

MCMCOutput.TunedMHSD.D_1=MHD_1_SD;
MCMCOutput.TunedMHSD.D_2=MHD_2_SD;
MCMCOutput.TunedMHSD.Noise=MHNoise_SD;

MCMCOutputSummary.TunedMHSD=MCMCOutput.TunedMHSD;

%Log likelihood for MAP parameters
MAPLikelihood=LogLikelihoodTwoStateMH(Traj,MCMCOutput.MAP);
MCMCOutput.MAPLikelihood=MAPLikelihood;
%BIC
NParameters=4;
BIC=-2*MAPLikelihood+NParameters*log(N);
MCMCOutput.BIC=BIC;
%AIC
AIC=-2*MAPLikelihood-2*NParameters;
MCMCOutput.AIC=AIC;
%Chib method for marginal likelihood
% if alg_parameters.MarginalLikelihood
%     MarginalLikelihood=ChibBindingNoiseMH(P,MAPLikelihood);
%     P.MarginalLikelihood=MarginalLikelihood;
% end

MarginalLikelihood=ChenTwoStateNoiseMH(MCMCOutput,MAPLikelihood);
MCMCOutput.MarginalLikelihood=MarginalLikelihood;
MCMCOutputSummary.MarginalLikelihood=MarginalLikelihood;

MCMCOutputSummary.TotalMovesD_1=TotalMovesD_1;
MCMCOutputSummary.TotalMovesD_2=TotalMovesD_2;


MCMCOutputSummary.D_1_chain=D_1_chain;
MCMCOutputSummary.D_2_chain=D_2_chain;
MCMCOutputSummary.noise_chain=noise_chain;
MCMCOutputSummary.p_12_chain=p_12_chain;
MCMCOutputSummary.p_21_chain=p_21_chain;

MCMCOutputSummary.Parameters={'D_0','D_1','p_01','p_10','sigma2'};

MCMCOutputSummary.z_mean=z_mean;
MCMCOutputSummary.z_post=z_mean-1;

MCMCOutputSummary.LogLikelihood=LogLikelihood;


MCMCOutputSummary.MAP=MCMCOutput.MAP;

MCMCOutputSummary.FinalState=MCMCOutput.FinalState;

MCMCOutputSummary.ParameterPosteriorSamples=MCMCOutput.ParameterPosteriorSamples;


MCMCOutputSummary.alg_parameters=alg_parameters;
MCMCOutputSummary.prior=prior;

MCMCOutputSummary.FixedNoise=noise;

MCMCOutputSummary.Traj=Traj;













end
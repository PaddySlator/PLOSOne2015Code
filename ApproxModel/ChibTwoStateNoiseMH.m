function [MarginalLikelihood,PReduced]=ChibTwoStateNoiseMH(P,LogLikelihood)
% Calculate Marginal likelkihood  for two-state diffusion with measurement
% noise (approximate likelihood) model using Chibs's method
% 
% (This is mentioned in Slator et al., PLOS ONE, 2015, but not used in the paper.)
%
% Reference: Chib S, Jeliazkov I. Marginal likelihood from the Metropolis?Hastings output.
% Journal of the American Statistical Association. 2001;96(453):270?281. 
% doi: 10.1198/016214501750332848.
%
% Paddy Slator, Warwick Systems Biology Centre


burn_in=P.alg_parameters.burn_in;
thin=P.alg_parameters.thin;
MCMC_steps=P.alg_parameters.MCMC_steps;
O=P.O;

N=length(O.Y)-1;
X=[O.Y(1:end-1,1), O.Y(1:end-1,2)];
dX=diff([O.Y(:,1),O.Y(:,2)]);
Dt=diff(O.Y(:,3));
R=diff(O.Y(:,1)).^2 + diff(O.Y(:,2)).^2;

D_1_MAP=P.MAP(1);
D_2_MAP=P.MAP(2);
p_12_MAP=P.MAP(3);
p_21_MAP=P.MAP(4);
noise_MAP=P.MAP(5);

a_1=P.prior.a_1;
b_1=P.prior.b_1;
a_2=P.prior.a_2;
b_2=P.prior.b_2;
D_max=P.prior.D_max;

NPostSamples=MCMC_steps/thin-burn_in/thin;

%samples from the posterior distriubtion
D_1Post=P.ParameterPosteriorSamples(:,1);
D_2Post=P.ParameterPosteriorSamples(:,2);
p_12Post=P.ParameterPosteriorSamples(:,3);
p_21Post=P.ParameterPosteriorSamples(:,4);
noisePost=P.ParameterPosteriorSamples(:,5);
zPost=P.zPosteriorSamples;

MHD_1_SD=P.TunedMHSD.D_1;
MHD_2_SD=P.TunedMHSD.D_2;
MHNoise_SD=P.TunedMHSD.Noise;

%First term
%\pi(D_1|X)

%Reduced sampler with D_1=D_1*
% %set initial values as MAP values
% initial_values.D_1=D_1_MAP;
% initial_values.D_2=D_2_MAP;
% initial_values.p_12=p_12_MAP;
% initial_values.p_21=p_21_MAP;
% initial_values.noise=noise_MAP;
% initial_values.z=round(zPost);
initial_values=[];

OReduced=O;
OReduced.parameters=D_1_MAP;
onchains=[0 1 1 1 1 1];

k=1;
alg_parameters.Reduced=P.alg_parameters;
alg_parameters.Reduced.MarginalLikelihood=0;
alg_parameters.Reduced.swap=0;
PReduced{k}=BindingNoiseMH(OReduced,alg_parameters.Reduced,P.prior,onchains,initial_values,[]);

D_2Reduced=PReduced{k}.ParameterPosteriorSamples(:,2);
p_12Reduced=PReduced{k}.ParameterPosteriorSamples(:,3);
p_21Reduced=PReduced{k}.ParameterPosteriorSamples(:,4);
noiseReduced=PReduced{k}.ParameterPosteriorSamples(:,5);
zReduced=PReduced{k}.zPosteriorSamples;
k=k+1;

%at each step of reduced MCMC run draw D_1 from q(D_1* -> D_1)
D_1Reduced=zeros(NPostSamples,1);
for i=1:NPostSamples
    D_1Reduced(i)=D_1_MAP+randn*MHD_1_SD;
end

Numerator=zeros(NPostSamples,1);
Denominator=zeros(NPostSamples,1);

for i=1:NPostSamples
       
    %numerator
    LogLikelihoodMAP=-sum(log(2*pi*(2*D_1_MAP*Dt(zPost(i,:)==1) + 2*noisePost(i))))...
        -sum(sum(dX(zPost(i,:)==1,:).^2,2)./(4*noisePost(i)+4*D_1_MAP*Dt(zPost(i,:)==1)));
    
    LogLikelihoodPosterior=-sum(log(2*pi*(2*D_1Post(i)*Dt(zPost(i,:)==1) + 2*noisePost(i))))...
        -sum(sum(dX(zPost(i,:)==1,:).^2,2)./(4*noisePost(i)+4*D_1Post(i)*Dt(zPost(i,:)==1)));
    
        
    q=normpdf(D_1_MAP,D_1Post(i),MHD_1_SD);
    
    Numerator(i)=...
    ...acceptance probability
        min(0,LogLikelihoodMAP-LogLikelihoodPosterior)...    
    ...jumping probability
        +log(q);
        
    
    
    %denominator
    LogLikelihoodReducedPosterior=-sum(log(2*pi*(2*D_1Reduced(i)*Dt(zReduced(i,:)==1) + 2*noiseReduced(i))))...
        -sum(sum(dX(zReduced(i,:)==1,:).^2,2)./(4*noiseReduced(i)+4*D_1Reduced(i)*Dt(zReduced(i,:)==1)));
    
    LogLikelihoodReducedMAP=-sum(log(2*pi*(2*D_1_MAP*Dt(zReduced(i,:)==1) + 2*noiseReduced(i))))...
        -sum(sum(dX(zReduced(i,:)==1,:).^2,2)./(4*noiseReduced(i)+4*D_1_MAP*Dt(zReduced(i,:)==1)));
            
    Denominator(i)=min(0,LogLikelihoodReducedPosterior-LogLikelihoodReducedMAP); 
      
end

PosteriorOrdinate(1)=log(sum(exp(Numerator)))-log(sum(exp(Denominator)));

PosteriorOrdinate(1)

ANum=max(Numerator);
ADen=max(Denominator);

PosteriorOrdinate(1)=ANum+log(sum(exp(Numerator-ANum)))...
    -(ADen+log(sum(exp(Denominator-ADen))));

PosteriorOrdinate(1)

PosteriorOrdinate(1)

%second term
%\pi(D_2*|D_1*,X)
%set previous reduced run (D_1=D_1*) as new posterior samples
D_2Post=D_2Reduced;
p_12Post=p_12Reduced;
p_21Post=p_21Reduced;
noisePost=noiseReduced;
zPost=zReduced;

%Reduced sampler with D_1=D_1*,D_2=D_2*
OReduced=O;
OReduced.parameters=[D_1_MAP D_2_MAP];
onchains=[0 0 1 1 1 1];

alg_parameters.Reduced=P.alg_parameters;
alg_parameters.Reduced.MarginalLikelihood=0;
alg_parameters.Reduced.swap=0;
PReduced{k}=BindingNoiseMH(OReduced,alg_parameters.Reduced,P.prior,onchains,initial_values,[]);

p_12Reduced=PReduced{k}.ParameterPosteriorSamples(:,3);
p_21Reduced=PReduced{k}.ParameterPosteriorSamples(:,4);
noiseReduced=PReduced{k}.ParameterPosteriorSamples(:,5);
zReduced=PReduced{k}.zPosteriorSamples;
k=k+1;

%at each step of reduced MCMC run draw D_2 from q(D_2* -> D_2)
D_2Reduced=zeros(NPostSamples,1);
for i=1:NPostSamples
    D_2Reduced(i)=D_2_MAP+randn*MHD_2_SD;
end

Numerator=zeros(NPostSamples,1);
Denominator=zeros(NPostSamples,1);

for i=1:NPostSamples
       
    %numerator
    LogLikelihoodMAP=-sum(log(2*pi*(2*D_2_MAP*Dt(zPost(i,:)==2) + 2*noisePost(i))))...
        -sum(sum(dX(zPost(i,:)==2,:).^2,2)./(4*noisePost(i)+4*D_2_MAP*Dt(zPost(i,:)==2)));
        
    LogLikelihoodPosterior=-sum(log(2*pi*(2*D_2Post(i)*Dt(zPost(i,:)==2) + 2*noisePost(i))))...
        -sum(sum(dX(zPost(i,:)==2,:).^2,2)./(4*noisePost(i)+4*D_2Post(i)*Dt(zPost(i,:)==2)));
        
    q=normpdf(D_2_MAP,D_2Post(i),MHD_2_SD);
    
    Numerator(i)=...
    ...acceptance probability
        min(0,LogLikelihoodMAP-LogLikelihoodPosterior)...    
    ...jumping probability
        +log(q);
        
    
    %denominator
    LogLikelihoodReducedPosterior=-sum(log(2*pi*(2*D_2Reduced(i)*Dt(zReduced(i,:)==2) + 2*noiseReduced(i))))...
        -sum(sum(dX(zReduced(i,:)==2,:).^2,2)./(4*noiseReduced(i)+4*D_2Reduced(i)*Dt(zReduced(i,:)==2)));
    
    LogLikelihoodReducedMAP=-sum(log(2*pi*(2*D_2_MAP*Dt(zReduced(i,:)==2) + 2*noiseReduced(i))))...
        -sum(sum(dX(zReduced(i,:)==2,:).^2,2)./(4*noiseReduced(i)+4*D_2_MAP*Dt(zReduced(i,:)==2)));
            
    Denominator(i)=min(0,LogLikelihoodReducedPosterior-LogLikelihoodReducedMAP);            
end

PosteriorOrdinate(2)=log(sum(exp(Numerator)))-log(sum(exp(Denominator)));

ANum=max(Numerator);
ADen=max(Denominator);

PosteriorOrdinate(2)=ANum+log(sum(exp(Numerator-ANum)))...
    -(ADen+log(sum(exp(Denominator-ADen))));

%third term
%\pi(sigma2*|D_1*,D_2*,X)
%set previous reduced run (D_1=D_1*,D_2=D_2*) as new posterior samples
p_12Post=p_12Reduced;
p_21Post=p_21Reduced;
noisePost=noiseReduced;
zPost=zReduced;

%Reduced sampler with D_1=D_1*,D_2=D_2*,sigma2=sigma2*
OReduced=O;
OReduced.parameters=[D_1_MAP D_2_MAP];
OReduced.noise_parameters.variance=noise_MAP;
onchains=[0 0 1 1 1 0];

alg_parameters.Reduced=P.alg_parameters;
alg_parameters.Reduced.MarginalLikelihood=0;
alg_parameters.Reduced.swap=0;
PReduced{k}=BindingNoiseMH(OReduced,alg_parameters.Reduced,P.prior,onchains,initial_values,[]);

p_12Reduced=PReduced{k}.ParameterPosteriorSamples(:,3);
p_21Reduced=PReduced{k}.ParameterPosteriorSamples(:,4);
zReduced=PReduced{k}.zPosteriorSamples;
k=k+1;

%at each step of reduced MCMC run draw sigma2 from q(sigma2* -> sigma2)
noiseReduced=zeros(NPostSamples,1);
for i=1:NPostSamples
    noiseReduced(i)=noise_MAP+randn*MHNoise_SD;
end

Numerator=zeros(NPostSamples,1);
Denominator=zeros(NPostSamples,1);

for i=1:NPostSamples
       
    %numerator
    LogLikelihoodMAP=...z==1
        -sum(log(2*pi*(2*D_1_MAP*Dt(zPost(i,:)==1) + 2*noise_MAP)))...
        -sum(sum(dX(zPost(i,:)==1,:).^2,2)./(4*noise_MAP+4*D_1_MAP*Dt(zPost(i,:)==1)))...
    ...z==2
        -sum(log(2*pi*(2*D_2_MAP*Dt(zPost(i,:)==2) + 2*noise_MAP)))...
        -sum(sum(dX(zPost(i,:)==2,:).^2,2)./(4*noise_MAP+4*D_2_MAP*Dt(zPost(i,:)==2)));
        
    LogLikelihoodPosterior=...z==1
        -sum(log(2*pi*(2*D_1_MAP*Dt(zPost(i,:)==1) + 2*noisePost(i))))...
        -sum(sum(dX(zPost(i,:)==1,:).^2,2)./(4*noisePost(i)+4*D_1_MAP*Dt(zPost(i,:)==1)))...
    ...z==2
        -sum(log(2*pi*(2*D_2_MAP*Dt(zPost(i,:)==2) + 2*noisePost(i))))...
        -sum(sum(dX(zPost(i,:)==2,:).^2,2)./(4*noisePost(i)+4*D_2_MAP*Dt(zPost(i,:)==2)));
        
    q=normpdf(noise_MAP,noisePost(i),MHNoise_SD);
    
    Numerator(i)=...
    ...acceptance probability
        min(0,LogLikelihoodMAP-LogLikelihoodPosterior)...    
    ...jumping probability
        +log(q);
        
    
    %denominator
        LogLikelihoodReducedPosterior=...z==1
        -sum(log(2*pi*(2*D_1_MAP*Dt(zReduced(i,:)==1) + 2*noiseReduced(i))))...
        -sum(sum(dX(zReduced(i,:)==1,:).^2,2)./(4*noiseReduced(i)+4*D_1_MAP*Dt(zReduced(i,:)==1)))...
    ...z==2
        -sum(log(2*pi*(2*D_2_MAP*Dt(zReduced(i,:)==2) + 2*noiseReduced(i))))...
        -sum(sum(dX(zReduced(i,:)==2,:).^2,2)./(4*noiseReduced(i)+4*D_2_MAP*Dt(zReduced(i,:)==2)));
    
    LogLikelihoodReducedMAP=...z==1
        -sum(log(2*pi*(2*D_1_MAP*Dt(zReduced(i,:)==1) + 2*noise_MAP)))...
        -sum(sum(dX(zReduced(i,:)==1,:).^2,2)./(4*noise_MAP+4*D_1_MAP*Dt(zReduced(i,:)==1)))...
    ...z==2
        -sum(log(2*pi*(2*D_2_MAP*Dt(zReduced(i,:)==2) + 2*noise_MAP)))...
        -sum(sum(dX(zReduced(i,:)==2,:).^2,2)./(4*noise_MAP+4*D_2_MAP*Dt(zReduced(i,:)==2)));
                   
    Denominator(i)=min(0,LogLikelihoodReducedPosterior-LogLikelihoodReducedMAP);            
end

PosteriorOrdinate(3)=log(sum(exp(Numerator)))-log(sum(exp(Denominator)));




%fourth term
%\pi(p_12|D_0*,D_1*,sigma2*,X)

%set previous reduced run (D_1=D_1*,D_2=D_2*,sigma2=sigma2*) as new posterior samples
p_12Post=p_12Reduced;
p_21Post=p_21Reduced;
zPost=zReduced;

p_12ConditionalDensity=zeros(NPostSamples,1);

%zReduced are samples from \pi(z|D_0*,D_1*,sigma2*,X)
for i=1:NPostSamples
    z=zPost;
    
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
    
    p_12ConditionalDensity(i)=betapdf(p_12_MAP,a_1+n_12,b_1+n_11);
end
A=max(p_12ConditionalDensity);
PosteriorOrdinate(4)=A+log(sum(exp(p_12ConditionalDensity-A)));

%fifth term
%\pi(p_12|D_0*,D_1*,sigma2*,X)


%Reduced sampler with D_1=D_1*,D_2=D_2*,sigma2=sigma2*,p_01=p_01^*
OReduced=O;
OReduced.parameters=[D_1_MAP D_2_MAP p_12_MAP];
OReduced.noise_parameters.variance=noise_MAP;
onchains=[0 0 0 1 1 0];

alg_parameters.Reduced=P.alg_parameters;
alg_parameters.Reduced.MarginalLikelihood=0;
alg_parameters.Reduced.swap=0;
PReduced{k}=BindingNoiseMH(OReduced,alg_parameters.Reduced,P.prior,onchains,initial_values,[]);

p_21Reduced=PReduced{k}.ParameterPosteriorSamples(:,4);
zReduced=PReduced{k}.zPosteriorSamples;
k=k+1;

p_21ConditionalDensity=zeros(NPostSamples,1);

%zReduced are samples from \pi(z|D_0*,D_1*,sigma2*,p_01^*,X)
for i=1:NPostSamples
    z=zReduced;
    
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
    
    p_21ConditionalDensity(i)=betapdf(p_21_MAP,a_2+n_21,b_2+n_22);
end

A=max(p_21ConditionalDensity);
PosteriorOrdinate(5)=A+log(sum(exp(p_21ConditionalDensity-A)));

PosteriorOrdinate

if isempty(LogLikelihood)
    LogLikelihood=LogLikelihoodBindingMH(O,P.MAP)
end

LogPrior=-2*log(D_max)...
    + log(betapdf(p_12_MAP,a_1,b_1))...
    + log(betapdf(p_21_MAP,a_2,b_2)) 

MarginalLikelihood=LogLikelihood+LogPrior-sum(PosteriorOrdinate);

end
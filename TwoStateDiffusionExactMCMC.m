function [MCMCOutput,MCMCOutputSummary] = TwoStateDiffusionExactMCMC(Traj,alg_parameters,prior,onchains,initial_values,figures)
% MCMC inference for two-state diffusion with measurement noise model
% (using the exact likelihood)
% see Slator et al., PLOS ONE, 2015
% Paddy Slator, Warwick Systems Biology Centre

%unpack trajectory data
N=length(Traj.Y)-1;
X=[Traj.Y(1:end,1), Traj.Y(1:end,2)];
dX=diff([Traj.Y(:,1),Traj.Y(:,2)]);
R=dX(:,1).^2+dX(:,2).^2;
Dt=diff(Traj.Y(:,3));

%which MCMC chains are on/off (useful for debugging etc.) 
onD_0=onchains(1);
onD_1=onchains(2);
onp_01=onchains(3);
onp_10=onchains(4);
onnoise=onchains(5);
onU=onchains(6);
onz=onchains(7);

MCMC_steps=alg_parameters.MCMC_steps;
burn_in=alg_parameters.burn_in;
thin=alg_parameters.thin;

if burn_in >= MCMC_steps
    error('burn in longer than total MCMC steps')
end

%prior parameters
D_max=prior.D_max;
D_min=prior.D_min;
a_D_0=prior.a_D_0;
b_D_0=prior.b_D_0;
a_D_1=prior.a_D_1;
b_D_1=prior.b_D_1;
a_0=prior.a_0;
b_0=prior.b_0;
a_1=prior.a_1;
b_1=prior.b_1;


if onnoise    
    noise=prior.FixedNoise;
end

%MCMC chains
D_0_chain=zeros(MCMC_steps/thin,1);
D_1_chain=zeros(MCMC_steps/thin,1);
p_01_chain=zeros(MCMC_steps/thin,1);
p_10_chain=zeros(MCMC_steps/thin,1);
noise_chain=zeros(MCMC_steps/thin,1);
z_chain=zeros(MCMC_steps/thin,N);
Ux_chain=zeros(MCMC_steps/thin,N+1);
Uy_chain=zeros(MCMC_steps/thin,N+1);

%initialise MCMC chains, either using given initial_values, 
%or from the prior, or overdispered (e.g. Gelman stat), or using simulated
%value if that MCMC chain is turned off.

if onD_0    
    if isfield(initial_values,'D_0') %given initial value
        D_0=initial_values.D_0;
    else
        if alg_parameters.overdisp %overdispersed for Gelman
            D_0 = D_min+betarnd(0.1,0.1)*D_max;
        else %from uniform prior
            D_0 = D_max*rand;
        end
    end
else %use simulated value (if parameter turned off)    
    D_0=Traj.parameters(1);
end

if onD_1
    if isfield(initial_values,'D_1')
        D_1=initial_values.D_1;
    else
        if alg_parameters.overdisp        
            D_1=D_min+betarnd(0.1,0.1)*(D_max-D_min);
        else
            D_1 = D_max*rand;
        end
    end
else
    D_1=Traj.parameters(2);
end

if onp_01
    if isfield(initial_values,'p_01')
        p_01=initial_values.p_01;
    else
        if alg_parameters.overdisp
            p_01=rand;
        else
            p_01 = betarnd(a_0, b_0);
        end
    end
else
    p_01 = Traj.parameters(3);
end

if onp_10
    if isfield(initial_values,'p_10')
        p_10=initial_values.p_10;
    else
        if alg_parameters.overdisp
            p_10 = rand;
        else
            p_10 = betarnd(a_1,b_1);
        end
    end
else
    p_10=Traj.parameters(4);
end

%initialise hidden state z
if onz %simulate a MC using p_01,p_10
    if isfield(initial_values,'z')
        z=initial_values.z;
    else
        pi_0=p_10/(p_10 + p_01);
        pi_1= 1 - pi_0;
                       
        z=zeros(N,1);
        
        if rand < pi_0
            z(1) = 0;
        else
            z(1) = 1;
        end
        
        for i=2:N
            if z(i-1) == 0 && rand < p_01
                %0->1 event
                z(i) = 1;
            elseif z(i-1) == 1 && rand < p_10
                %1->0 event
                z(i) = 0;                
            else
                %no event
                z(i)=z(i-1);
            end
        end 
    end
else %use simulated value (if z chain turned off)
    z=Traj.z;
end

%initialise hidden state U
if onU
    if isfield(initial_values,'U') %if given initial value
        U=initial_values.U;
    else %start at U=X
        U=X;
    end
else %if U chain off 
    if isfield(Traj,'U') %use simulated value
        U=Traj.U;
    else
        U=X;
    end
end


step=1;
D_0_chain(step)=D_0;
D_1_chain(step)=D_1;
p_01_chain(step)=p_01;
p_10_chain(step)=p_10;
z_chain(step,:)=z;
Ux_chain(step,:)=U(:,1)';
Uy_chain(step,:)=U(:,2)';


if onnoise
    noise_chain(step)=noise;    
end


for step=2:MCMC_steps
    
    %%parameter updates
    
    if onD_0
        dU=diff(U);        
        %calculate time-scaled square displacements 
        R=sum(dU.^2,2)./Dt;               
                
        if (sum(z==0)+a_D_0==0) || (b_D_0+0.25*sum(R(z==0))==0)
            D_0=rand*D_max;
        else
            D_0 = 1/(gamrnd(sum(z==0)+1+a_D_0, 1/(b_D_0+0.25*sum(R(z==0)))));
            while D_0>D_max
                D_0 = 1/(gamrnd(sum(z==0)+1+a_D_0, 1/(b_D_0+0.25*sum(R(z==0)))));                       
            end
        end
    end
    
    if onD_1
        dU=diff(U);
        
        %calculate time-scaled square displacements 
        R=sum(dU.^2,2)./Dt;
        
        if ( sum(z==1)+a_D_1==0) || (b_D_1+0.25*sum(R(z==1))==0)
            D_1=rand*D_max;
        else
            D_1 = 1/(gamrnd(sum(z==1)+1+a_D_1, 1/(b_D_1+0.25*sum(R(z==1)))));
            while D_1>D_max
                D_1 = 1/(gamrnd(sum(z==1)+1+a_D_1, 1/(b_D_1+0.25*sum(R(z==1)))));
            end
        end
        
    end
    
    if onp_01 || onp_10       
        %calculate number of each possible switching/non-switching event
        n_01=sum(z(1:end-1)==0 & z(2:end)==1);
        n_00=sum(z(1:end-1)==0 & z(2:end)==0);      
        n_10=sum(z(1:end-1)==1 & z(2:end)==0);       
        n_11=sum(z(1:end-1)==1 & z(2:end)==1);        
    end
    
    if onp_01
        p_01=betarnd(a_0+n_01,b_0+n_00);
    end
    
    if onp_10
        p_10=betarnd(a_1+n_10,b_1+n_11);
    end                   
    
    %%hidden state updates
    
    if onz
        P_t = [ 1 - p_01, p_01
            p_10, 1 - p_10];
        
        dU=diff(U);
        
        i=1;        
        P_0=exp(-(sum(dU(i,:).^2))/(4*D_0*Dt(i)))/((4*pi*D_0*Dt(i)))...
            *P_t(1,z(i+1)+1);
                    
        P_1=exp(-(sum(dU(i,:).^2))/(4*D_1*Dt(i)))/((4*pi*D_1*Dt(i)))...
            *P_t(2,z(i+1)+1);
        
        if rand < P_0/(P_0+P_1)
            z(i)=0;
        else
            z(i)=1;
        end
        
        for i=2:N-1                        
            P_0=P_t(z(i-1)+1,1)...
                *exp(-(sum(dU(i,:).^2))/(4*D_0*Dt(i)))/((4*pi*D_0*Dt(i)))...
                *P_t(1,z(i+1)+1);                                   
            
            P_1=P_t(z(i-1)+1,2)...
                *exp(-(sum(dU(i,:).^2))/(4*D_1*Dt(i)))/((4*pi*D_1*Dt(i)))...
                *P_t(2,z(i+1)+1);
            
            if rand < P_0/(P_0+P_1)
                z(i)=0;
            else
                z(i)=1;
            end           
        end
        
        i=N;       
        P_0=P_t(z(i-1)+1,1)...
            *exp(-(sum(dU(i,:).^2))/(4*D_0*Dt(i)))/((4*pi*D_0*Dt(i)));
                    
        P_1=P_t(z(i-1)+1,2)...
            *exp(-(sum(dU(i,:).^2))/(4*D_1*Dt(i)))/((4*pi*D_1*Dt(i)));
        
        if rand < P_0/(P_0+P_1)
            z(i)=0;
        else
            z(i)=1;
        end
        
    end
    
    if onU
        D(1)=D_0;
        D(2)=D_1;
        i=1;
        tau= 1/prior.sigma_U + 1/(2*D(z(i)+1)*Dt(i)) + 1/noise;
        mu=(1/tau)*(...
            prior.mu_U/prior.sigma_U...
            + U(i+1,:)./(2*D(z(i)+1)*Dt(i))...
            + X(i,:)/noise...
            );
        U(i,:)= mu + sqrt(1/tau)*randn(1,2);
        
        for i=2:N
            tau=1/(2*D(z(i-1)+1)*Dt(i-1)) + 1/(2*D(z(i)+1)*Dt(i)) + 1/noise;
            mu=(1/tau)*(...
                U(i-1,:)./(2*D(z(i-1)+1)*Dt(i-1))...
                + U(i+1,:)./(2*D(z(i)+1)*Dt(i))...
                + X(i,:)/noise...
                );
            U(i,:)= mu + sqrt(1/tau)*randn(1,2);
            
        end
        i=N+1;
        tau=1/(2*D(z(i-1)+1)*Dt(i-1)) + 1/noise;
        mu=(1/tau)*(...
            U(i-1,:)./(2*D(z(i-1)+1)*Dt(i-1))...
            + X(i,:)/noise...
            );
        U(i,:)= mu + sqrt(1/tau)*randn(1,2);
    end     
        
    %store variable values at the sampling rate
    if mod(step,thin)==0
        D_0_chain(step/thin)=D_0;
        D_1_chain(step/thin)=D_1;
        p_01_chain(step/thin)=p_01;
        p_10_chain(step/thin)=p_10;
        noise_chain(step/thin)=noise;
        if onnoise
            noise_chain(step/thin)=noise;
        end 
        z_chain(step/thin,:)=z;
        Ux_chain(step/thin,:)=U(:,1)';
        Uy_chain(step/thin,:)=U(:,2)';        
    end
 
end %end of MCMC loop

%force D_0<D_1 by swapping chains if necessary
%swap 0,1 in chains if mean(D_0) > mean(D_1)
if mean(D_0_chain(burn_in/thin:end)) < mean(D_1_chain(burn_in/thin:end))
    [D_1_chain,D_0_chain]=deal(D_0_chain,D_1_chain);
    [p_10_chain,p_01_chain]=deal(p_01_chain,p_10_chain);
    
    for i=2:MCMC_steps/thin
        for j=1:N
            if z_chain(i,j)==0
                z_chain(i,j)=1;
            else
                z_chain(i,j)=0;
            end
        end
    end
end

MCMCOutput.ParameterChains=[D_0_chain D_1_chain p_10_chain p_01_chain];

MCMCOutput.ParameterPosteriorSamples=[D_0_chain((burn_in/thin+1):end)...
        D_1_chain((burn_in/thin+1):end)...
        p_01_chain((burn_in/thin+1):end)...
        p_10_chain((burn_in/thin+1):end)];

%posterior probability of being in state z=1
z_mean= mean(z_chain((burn_in/thin):end,:));

%mean U posterior positions
Ux_mean=mean(Ux_chain((burn_in/thin):end,:));
Uy_mean=mean(Uy_chain((burn_in/thin):end,:));


MCMCOutput.Ux_chain=Ux_chain;
MCMCOutput.Uy_chain=Uy_chain;

MCMCOutput.z_chain=z_chain;

MCMCOutput.z_mean=z_mean;
MCMCOutput.Ux_mean=Ux_mean;
MCMCOutput.Uy_mean=Uy_mean;


MCMCOutput.Traj=Traj;
MCMCOutput.onchains=onchains;

MCMCOutput.alg_parameters=alg_parameters;
MCMCOutput.prior=prior;

MCMCOutput.Parameters={'D_0','D_1','p_01','p_10'};
%useful for plot axes
MCMCOutput.ParameterLabels={'D_0','D_1','p_{01}','p_{10}'};

if onnoise
    MCMCOutput.Model='MCMC run for two-state diffusion model with Gaussian noise';
else
    MCMCOutput.Model='MCMC run for two-state diffusion model: no noise';
end

%these are useful for continuing an MCMC run
MCMCOutput.FinalState.D_0=D_0_chain(end);
MCMCOutput.FinalState.D_1=D_0_chain(end);
MCMCOutput.FinalState.p_01=p_01_chain(end);
MCMCOutput.FinalState.p_10=p_10_chain(end);
MCMCOutput.FinalState.noise=noise_chain(end);
MCMCOutput.FinalState.z=z_chain(end,:);
MCMCOutput.FinalState.U=[Ux_chain(end,:)' Uy_chain(end,:)'];

%summary of output (less memory)
MCMCOutputSummary=MCMCOutput;
%these are expensive to save!!!
MCMCOutputSummary=rmfield(MCMCOutputSummary,'z_chain');
MCMCOutputSummary=rmfield(MCMCOutputSummary,'Uy_chain');
MCMCOutputSummary=rmfield(MCMCOutputSummary,'Ux_chain');


   





end
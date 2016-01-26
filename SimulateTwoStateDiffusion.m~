function Traj = SimulateTwoStateDiffusion(parameters,noise_parameters,N,prior,ploton)
% 2D Brownian motion with two diffusion coefficents (e.g. Das et al., Plos Comp Bio., 2009)
% but with the addition of Gaussian noise (e.g. Slator et al., PLOS ONE, 2015)
% Paddy Slator, Warwick Systems Biology Centre

if isempty(noise_parameters)
    onnoise=0;
    noise=0;
else
    onnoise=1;
    noise=noise_parameters.variance;
end

D_0 = parameters(1);
D_1 = parameters(2); 
p_01 = parameters(3); 
p_10 = parameters(4);

% Time steps
switch parameters(5)
    case 1
        Dt=parameters(6)*ones(N,1);
    case 2 % Gamma distribution.
        Dt=gamrnd(parameters(7)*ones(N,1),1/parameters(8));
        if ploton
            figure;hist(Dt,100);
            xlabel('Dt');
            title('Distribution of displacement times');
        end
end

% Compute stationary probabilities for the 2-state MC
pi_1=p_10/(p_10 + p_01);
pi_2= 1 - pi_1;

z=zeros(N,1);

%simulate first timepoint based on stationary probs
if rand < pi_1
    z(1) = 0;
else
    z(1) = 1;
end

%also note the switchpoints - useful for plotting 
switch12=[];
switch21=[];
for i=2:N
    if z(i-1) == 0 && rand < p_01   
        % 0->1 event
        z(i) = 1;
        switch12=[switch12,i];
    elseif z(i-1) == 1 && rand < p_10
        % 1->0 event
        z(i) = 0;
        switch12=[switch12,i];
    else
        % no event
        z(i)=z(i-1);
    end
end


U=zeros(N+1,2);

if isempty(prior)
    U(1,:)=[0 0];
else
    U(1,:)=prior.mu_U + prior.sigma_U*randn(1,2);
end

D=[D_0 D_1];

%simulate particle path using hidden state sequence z
for i=1:N    
    U(i+1,:)=U(i,:) + sqrt(2*D(z(i)+1)*Dt(i))*randn(1,2);
end

%add measurement noise (Gaussian) to particle path
if onnoise
    X=zeros(N+1,2);
    for i=1:N+1
        X(i,:) = U(i,:)  + sqrt(noise).*randn(1,2);
    end
    Traj.U=U;
else
    X=U;
end


Traj.noise_parameters=noise_parameters;

%trajectory path
%(x coordinate, y coordinate, time)
Traj.Y=[X [0; cumsum(Dt)]];

Traj.z=z;

Traj.parameters=parameters;

Traj.type='Simulation';
if onnoise
        Traj.Model='Two-state diffusion model with Gaussian noise';    
else
    Traj.Model='Two-state diffusion model without noise';
end


    


end






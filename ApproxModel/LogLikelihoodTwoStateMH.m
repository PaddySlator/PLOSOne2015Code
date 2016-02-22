function LogLikelihood = LogLikelihoodTwoStateMH(Traj,parameters)
% Calculate log likelihood for approx two-state model (forward algorithm).
% see Slator et al., PLOS ONE, 2015
% also Das, Cairo, Coombs, PLOS Computational Biology, 2009
% Paddy Slator, Warwick Systems Biology Centre

D_1 = parameters(1);
D_2 = parameters(2);
p_12 = parameters(3);
p_21 = parameters(4);
noise = parameters(5);

N=length(Traj.Y)-1;
X=[Traj.Y(1:end-1,1), Traj.Y(1:end-1,2)];
dX=diff([Traj.Y(:,1),Traj.Y(:,2)]);
Dt=diff(Traj.Y(:,3));

pi_1=p_21/(p_21 + p_12);
pi_2= 1 - pi_1;
p_11= 1 - p_12;
p_22= 1 - p_21;

LL=zeros(N,2);

LL(:,1)=-log(2*pi*(2*D_1*Dt + 2*noise))...
    -sum(dX.^2,2)./(4*noise+4*D_1*Dt);

LL(:,2)=-log(2*pi*(2*D_2*Dt + 2*noise))...
    -sum(dX.^2,2)./(4*noise+4*D_2*Dt);

LogAlpha=zeros(N,2);
LogAlpha(1,:)=[log(pi_1) + LL(1,1),log(pi_2) + LL(1,2)];

for j=2:N
    b_1=max(LogAlpha(j-1,1) + log(p_11) + LL(j,1)...
        ,LogAlpha(j-1,2) + log(p_21) + LL(j,1));
    
    LogAlpha(j,1)=b_1...
        +log(exp(LogAlpha(j-1,1) + log(p_11)+ LL(j,1)-b_1)...
        +exp(LogAlpha(j-1,2) + log(p_21) + LL(j,1)-b_1));
    
    b_2=max(LogAlpha(j-1,1) + log(p_12)+ LL(j,2)...
        ,LogAlpha(j-1,2) + log(p_22) + LL(j,2));
    
    LogAlpha(j,2)=b_2...
        +log(exp(LogAlpha(j-1,1) + log(p_12)+ LL(j,2)-b_2)...
        +exp(LogAlpha(j-1,2) + log(p_22) + LL(j,2)-b_2));
end

b=max(LogAlpha(N,1),LogAlpha(N,2));

LogLikelihood=b+log(exp(LogAlpha(N,1)-b)+exp(LogAlpha(N,2)-b));

end
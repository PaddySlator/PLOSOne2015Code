function MarginalLikelihood = ChenOneStateNoiseMH(P,LogLikelihood)

burn_in=P.alg_parameters.burn_in;
thin=P.alg_parameters.thin;
MCMC_steps=P.alg_parameters.MCMC_steps;
O=P.O;

N=length(O.Y)-1;
X=[O.Y(1:end-1,1), O.Y(1:end-1,2)];
dX=diff([O.Y(:,1),O.Y(:,2)]);
Dt=diff(O.Y(:,3));
R=diff(O.Y(:,1)).^2 + diff(O.Y(:,2)).^2;

D_MAP=P_BM.MAP(1);
noise_MAP=P_BM.MAP(2);





end
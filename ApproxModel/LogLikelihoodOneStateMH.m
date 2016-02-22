function LogLikelihood = LogLikelihoodOneStateMH(O,parameters)

D=parameters(1);
noise=parameters(2);

N=length(O.Y)-1;
X=[O.Y(1:end-1,1), O.Y(1:end-1,2)];
dX=diff([O.Y(:,1),O.Y(:,2)]);
Dt=diff(O.Y(:,3));


% for i=1:N
%     LogLikelihood(i)=log(prod(normpdf(dX(i,:),0,2*(D*Dt(i)+noise))));
% end
% 
% sum(LogLikelihood)

LogLikelihood = -sum(log(2*pi*(2*D*Dt + 2*noise)))...
    -sum(sum(dX.^2,2)./(2*2*noise+4*D*Dt));



end
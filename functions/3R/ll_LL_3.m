function [logLik]=ll_LL_3(teta,times,sigmas)

B=reshape(teta(1:9),3,3);

L2=diag(teta(10:12));

L3=diag(teta(13:15));

K1 = pinv(B);
K2 = pinv(B*L2*B');
K3 = pinv(B*L3*B');
  
T1=times(1);
T2=times(2);
T3=times(3);
T=T1+T2+T3;

M=size(B,1);

logLik=-(-0.5*T*M*(log(2*pi))...
    +0.5*T1*log((det(K1))^2)-0.5*T1*trace(K1'*K1*sigmas{1})...
    +0.5*T2*log((det(K2))^2)-0.5*T2*trace(K2'*K2*sigmas{2})...
    +0.5*T3*log((det(K3))^2)-0.5*T3*trace(K3'*K3*sigmas{3}));   

end
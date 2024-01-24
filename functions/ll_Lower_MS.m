function [logLik]=ll_Lower_MS(teta,times,sigmas)

B=[teta(1) 0       0;
   teta(2) teta(3) 0;
   0       0       teta(4)]; % new 0 in (1,2)

Q2=[teta(5) 0       teta(8);
    teta(6) teta(7) 0;
    0       0       teta(9)]; % new 0 in (3,1)

Q3=[0        0        teta(12);
    teta(10) teta(11) teta(13);
    0        0        teta(14)];

Q4=[teta(15),0,0;
    teta(16),teta(17),teta(18);
    0,0,0]; % new 0 in (3,1) and in (1,2)


K1 = pinv(B);
K2 = pinv(B+Q2);
K3 = pinv(B+Q2+Q3);
K4 = pinv(B+Q2+Q3+Q4);
  
T1=times(1);
T2=times(2);
T3=times(3);
T4=times(4);
T=T1+T2+T3+T4;

M=size(B,1);

logLik=-(-0.5*T*M*(log(2*pi))...
    +0.5*T1*log((det(K1))^2)-0.5*T1*trace(K1'*K1*sigmas{1})...
    +0.5*T2*log((det(K2))^2)-0.5*T2*trace(K2'*K2*sigmas{2})...
    +0.5*T3*log((det(K3))^2)-0.5*T3*trace(K3'*K3*sigmas{3})...
    +0.5*T4*log((det(K4))^2)-0.5*T4*trace(K4'*K4*sigmas{4}));    

end
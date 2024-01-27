clear
clc
cd '/Users/samueleborsini/Library/Mobile Documents/com~apple~CloudDocs/Università/Economics and econometrics/II anno/Structural Macroeconometrics/Paper/final'
addpath functions/4R % functions folder
data=readtable("Data_Set.csv"); % data

%% estimation
lags=4; % number of lags of the VAR
DataSet=[table2array(data(:,"UM1")),table2array(data(:,"IP")),table2array(data(:,"UF1"))]; % selected variables of interest
M=size(DataSet,2); % number of variables

VAR_PI=cell(1,4);
for i = 1:lags
    VAR_PI{i} = nan(3);
end
VAR_c=nan(3,1);
VAR = varm('Constant',VAR_c,'AR',VAR_PI); % VAR setup

HorizonIRF=60; % time horizon of the IRF
repetitions=999; % bootstrap reptitions

%% dates
date=table2array(data(:,"date")); % dates

%% Dataset for each regime
TBs=[285,570,716]; % break dates
DataSets={DataSet(1:TBs(1),:),DataSet(TBs(1)+1-lags:TBs(2),:),DataSet(TBs(2)+1-lags:TBs(3),:),DataSet(TBs(3)+1-lags:end,:)}; % cell with the 4 different datasets

%% duplication matrix
DuplicationMatrix = zeros(M^2,0.5*M*(M+1));
DuplicationMatrix(1,1)=1;
DuplicationMatrix(2,2)=1;
DuplicationMatrix(3,3)=1;
DuplicationMatrix(4,2)=1;
DuplicationMatrix(5,4)=1;
DuplicationMatrix(6,5)=1;
DuplicationMatrix(7,3)=1;
DuplicationMatrix(8,5)=1;
DuplicationMatrix(9,6)=1;
mDD=pinv(DuplicationMatrix'*DuplicationMatrix)*DuplicationMatrix';
mNN=DuplicationMatrix*mDD;

%% VAR overall period
[EstVAR,EstSE,All.Lk,All.res] = estimate(VAR,DataSet); % VAR estimation
All.T = size(DataSet,1)-lags; % T

All.Const=EstVAR.Constant; % constant
All.Pi=EstVAR.AR; % Standard errors of the reduced form parameters (autoregressive parameters) (big PI) are available from the built-in estimation
All.Pi_se=EstSE.AR; % Pi s.e.
All.Const_se=EstSE.Constant; % const s.e.
for i=1:lags
    All.Pi_t_stat{i}=All.Pi{i}./All.Pi_se{i}; %t-stats of the coefficients of the i-th lag
end
All.Const_t_stat=All.Const./All.Const_se; % t-stats of the constant term
All.Companion_Matrix=[cell2mat(EstVAR.AR); eye(M*(lags-1),M*(lags-1)) zeros(M*(lags-1),M)]; % companion matrix

All.Sigma=EstVAR.Covariance; %var-cov matrix
All.Sigma_se=reshape(DuplicationMatrix*sqrt(diag(2/All.T*((mDD*kron(EstVAR.Covariance,EstVAR.Covariance)*(mDD)')))),3,3); % Standard errors of the reduced form parameters (covariance matrix)
All.Sigma_VC=(2/All.T*((mDD*kron(EstVAR.Covariance,EstVAR.Covariance)*(mDD)'))); %this is the 6x6 matrix (taking also the covariances and not only the elements on the diagonal
All.Sigma_Corr_Matrix= pinv(sqrt(diag(diag(EstVAR.Covariance))))*EstVAR.Covariance*pinv(sqrt(diag(diag(EstVAR.Covariance)))); % sigma expressed as correlation matrix
All.Sigma_t_stat = All.Sigma./All.Sigma_se; %t-stats (divide estimated coeff by their s.e.)

%% different regimes estimation
regimes=cell(1,4);
for i=1:4
    [EstVAR,EstSE,regimes{i}.Lk,regimes{i}.res] = estimate(VAR,DataSets{i}); % VAR estimation
    regimes{i}.T = size(DataSets{i},1)-lags; % T
    
    regimes{i}.Const=EstVAR.Constant; % constant
    regimes{i}.Pi=EstVAR.AR; % Standard errors of the reduced form parameters (autoregressive parameters) (big PI) are available from the built-in estimation
    regimes{i}.Pi_se=EstSE.AR; % Pi s.e.
    regimes{i}.Const_se=EstSE.Constant; % const s.e.
    for j=1:lags
        regimes{i}.Pi_t_stat{j}=regimes{i}.Pi{j}./regimes{i}.Pi_se{j}; %t-stats of the coefficients of the i-th lag
    end
    regimes{i}.Const_t_stat=regimes{i}.Const./regimes{i}.Const_se; % t-stats of the constant term
    regimes{i}.Companion_Matrix=[cell2mat(EstVAR.AR); eye(M*(lags-1),M*(lags-1)) zeros(M*(lags-1),M)]; % companion matrix
    
    regimes{i}.Sigma=EstVAR.Covariance; %var-cov matrix
    regimes{i}.Sigma_se=reshape(DuplicationMatrix*sqrt(diag(2/regimes{i}.T*((mDD*kron(EstVAR.Covariance,EstVAR.Covariance)*(mDD)')))),3,3); % Standard errors of the reduced form parameters (covariance matrix)
    regimes{i}.Sigma_VC=(2/regimes{i}.T*((mDD*kron(EstVAR.Covariance,EstVAR.Covariance)*(mDD)'))); %this is the 6x6 matrix (taking also the covariances and not only the elements on the diagonal
    regimes{i}.Sigma_Corr_Matrix= pinv(sqrt(diag(diag(EstVAR.Covariance))))*EstVAR.Covariance*pinv(sqrt(diag(diag(EstVAR.Covariance)))); % sigma expressed as correlation matrix
    regimes{i}.Sigma_t_stat = regimes{i}.Sigma./regimes{i}.Sigma_se; %t-stats (divide estimated coeff by their s.e.)
end

%% sigmas and times
times=[regimes{1}.T,regimes{2}.T,regimes{3}.T,regimes{4}.T]; % n obs each regime
sigmas={regimes{1}.Sigma,regimes{2}.Sigma,regimes{3}.Sigma,regimes{4}.Sigma}; %sigma of each regime

%% test on sigma
%test H0 (no change in sigma and Pi)
df=(3*13*4+6*4)-(3*13+6); %under H1 we have estimated 4 different PIs (3x13 each) and 4 different sigmas (symmetric 3x3 each); under H0 we have estimated a single Pi (3x13) and a single sigma (3x3 symmetric) 
pvalH0=1-chi2cdf(2*(regimes{1}.Lk+regimes{2}.Lk+regimes{3}.Lk+regimes{4}.Lk-All.Lk),df); %under H0, the VAR estimated in All is the specification
disp(['Testing sigma and Pi constant across periods gives a pvalue of ' num2str(pvalH0)]) 

%test H0' (no change in sigma given no change in Pi)
%separating the residuals in each volatility bloc
residuals_regime{1}=All.res(1:TBs(1)-lags,:);
residuals_regime{2}=All.res(TBs(1)-lags+1:TBs(2)-lags,:);
residuals_regime{3}=All.res(TBs(2)-lags+1:TBs(3)-lags,:);
residuals_regime{4}=All.res(TBs(3)-lags+1:end,:);

%estimating the sigma in each block
sigmas_test=cell(1,4);
for i=1:4
    sigmas_test{i}=(1/size(residuals_regime{i},1))*residuals_regime{i}'*residuals_regime{i};
end

%computing the loglikelihood (the same formula of the maximized log-likelihood, yet with 4 sigmas)
lk_H0_1=-0.5*All.T*M*log(2*pi)-0.5*times(1)*log(det(sigmas_test{1}))-0.5*times(2)*log(det(sigmas_test{2}))-0.5*times(3)*log(det(sigmas_test{3}))-0.5*times(4)*log(det(sigmas_test{4}))-0.5*All.T*M;

df=(3*13+6*4)-(3*13+6); %under H1' we have estimated a single Pi (3x13) and 4 different sigmas (symmetric 3x3 each); under H0' we have estimated a single Pi (3x13) and a single sigma (3x3 symmetric) 
pvalH0_1=1-chi2cdf(2*(lk_H0_1-All.Lk),df); %under H0', the VAR estimated in All is the specification
disp(['Testing sigma constant across periods given constant Pi gives a pvalue of ' num2str(pvalH0_1)]) 

%% SVAR (upper)
disp(' ')
disp('Upper (macro uncertainty endogenous and bidirectionality between macro and financial uncertainty):')
initialvalues=[0.0112;-0.1377;0.0003;0.7540;-0.0259;-0.0037;0.0470;-0.0047;-0.2927;0.0040;0.0535;-0.0318;0.1355;0.0014;-0.0839;-0.0066;0;0;0;0;0;0]; %initial values (the first 16 are the estimates obtain by ABCF, we set the new parameters initial value equal to 0)

options = optimset('MaxFunEvals',200000,'TolFun',1e-500,'MaxIter',200000,'TolX',1e-100,'Display', 'off'); % optimization settings
fun = @(teta) ll_Upper_MS(teta,times,sigmas); % function to optimization
[SVAR.params,nLk,~,~,~,SVAR.Hessian] = fminunc(fun,initialvalues,options); % optimization

SVAR.Lk=-nLk; %invert the sign of the likelihood
SVAR.SE=diag(pinv(SVAR.Hessian)).^0.5; %standard errors from the Hessian
SVAR.t_stat=SVAR.params./SVAR.SE; %t-stat

PVarl_LRTest = 1-chi2cdf(2*(regimes{1}.Lk+regimes{2}.Lk+regimes{3}.Lk+regimes{4}.Lk-SVAR.Lk),24-size(SVAR.params,1)); % overidentification test
disp(['The p-value of the overidentification restrictions is ' num2str(PVarl_LRTest)]) 
%crucial result above: we can reject the model in which macro uncertainty is endogeneous and we have bidirectionality

SVAR.B=[SVAR.params(1),SVAR.params(3),0;
        SVAR.params(2),SVAR.params(4),0;
        0,0,SVAR.params(5)];

SVAR.Q2=[SVAR.params(6),0,SVAR.params(10);
         SVAR.params(7),SVAR.params(9),0;
         SVAR.params(8),0,SVAR.params(11)];

SVAR.Q3=[0,0,SVAR.params(14);
         SVAR.params(12),SVAR.params(13),SVAR.params(15);
         0,0,SVAR.params(16)]; 

SVAR.Q4=[SVAR.params(17),SVAR.params(20),0;
        SVAR.params(18),SVAR.params(21),SVAR.params(22);
        SVAR.params(19),0,0];  %new matrix (last volatility period)

%recall SVAR.B is for the first period (GI)
SVAR.B2=SVAR.B+SVAR.Q2; %second period (GM)
SVAR.B3=SVAR.B2+SVAR.Q3; %third period (GR)
SVAR.B4=SVAR.B3+SVAR.Q4; %fourth period (COVID)

SVAR.V=pinv(SVAR.Hessian); % var-cov matrix
SVAR.SE_sums=SE_sums(SVAR); % s.e. of B, B2, B3 and B4

%sign normalization:
for i = 1:M
    if SVAR.B(i,i)<0
        SVAR.B(:,i)=-SVAR.B(:,i);
    end
end

for i = 1:M
    if SVAR.B2(i,i)<0
        SVAR.B2(:,i)=-SVAR.B2(:,i);
    end
end

for i = 1:M
    if SVAR.B3(i,i)<0
        SVAR.B3(:,i)=-SVAR.B3(:,i);
    end
end

for i = 1:M
    if SVAR.B4(i,i)<0
        SVAR.B4(:,i)=-SVAR.B4(:,i);
    end
end

%obtain Q2, Q3 and Q4 again after the sign normalization
SVAR.Q2=SVAR.B2-SVAR.B; 
SVAR.Q3=SVAR.B3-SVAR.B2;
SVAR.Q4=SVAR.B4-SVAR.B3;

%Rank condition check
V1=kron(SVAR.B,eye(M));
V2=kron(SVAR.B+ SVAR.Q2, eye(M));
V3=kron(SVAR.B+ SVAR.Q2+ SVAR.Q3, eye(M));
V4=kron(SVAR.B+ SVAR.Q2+ SVAR.Q3+ SVAR.Q4, eye(M));

% Calculates the matrix for checking the rank condition (full column rank)
RankMatrix=kron(eye(4),mDD)*[V1,   zeros(M^2,M^2),  zeros(M^2,M^2),    zeros(M^2,M^2);
                             V2,    V2,           zeros(M^2,M^2),    zeros(M^2,M^2);
                             V3,    V3,            V3,               zeros(M^2,M^2);   
                             V4,    V4,            V4,                  V4]; 
 
% Selection matrix for extracting the structural parameters                         
HSelection=zeros(4*M^2,size(SVAR.params,1));

HSelection(1,1)=1;
HSelection(2,2)=1;
HSelection(4,3)=1;
HSelection(5,4)=1;
HSelection(9,5)=1;

HSelection(10,6)=1;
HSelection(11,7)=1;
HSelection(12,8)=1;
HSelection(14,9)=1;
HSelection(16,10)=1;
HSelection(18,11)=1;

HSelection(20,12)=1;
HSelection(23,13)=1;
HSelection(25,14)=1;
HSelection(26,15)=1;
HSelection(27,16)=1;

HSelection(28,17)=1;
HSelection(29,18)=1;
HSelection(30,19)=1;
HSelection(31,20)=1;
HSelection(32,21)=1;
HSelection(34,21)=1;
HSelection(35,22)=1;

Jacobian=RankMatrix*HSelection;
Jacobian=2*Jacobian; % multiplicative scalar 2, see Fanelli 2015 OBES

% Report the rank of the matrix for checking the identification
disp(['The rank of the Jacobian is ' num2str(rank(Jacobian)) ' and the number of parameters is ' num2str(size(SVAR.params,1))])  %full col rank (22)

SVAR_U=SVAR; % storing the results

%% SVAR (lower)
disp(' ')
disp('Lower (macro uncertainty exogenous and unidirectionality from financial to macro uncertainty):')
initialvalues=[0.0112;-0.1377;0.7540;-0.0259;-0.0037;0.0470;-0.2927;0.0040;0.0535;-0.0318;0.1355;0.0014;-0.0839;-0.0066;0;0;0;0]; %initial values (the first 16 are the estimates obtain by ABCF, we set the new parameters initial value equal to 0)

options = optimset('MaxFunEvals',200000,'TolFun',1e-500,'MaxIter',200000,'TolX',1e-100,'Display', 'off'); % optimization settings
fun = @(teta) ll_Lower_MS(teta,times,sigmas);
[SVAR.params,nLk,~,~,~,SVAR.Hessian] = fminunc(fun,initialvalues,options);

SVAR.Lk=-nLk; %invert the sign of the likelihood
SVAR.SE=diag(pinv(SVAR.Hessian)).^0.5; %standard errors from the Hessian
SVAR.t_stat=SVAR.params./SVAR.SE; %t-stat

PVarl_LRTest = 1-chi2cdf(2*(regimes{1}.Lk+regimes{2}.Lk+regimes{3}.Lk+regimes{4}.Lk-SVAR.Lk),24-size(SVAR.params,1)); % overidentification test
disp(['The p-value of the overidentification restrictions is ' num2str(PVarl_LRTest)]) 
%crucial result above: we can reject the model in which macro uncertainty is exogenous but we have unidirectionality

SVAR.B=[SVAR.params(1),0,0;
        SVAR.params(2),SVAR.params(3),0;
        0,0,SVAR.params(4)];

SVAR.Q2=[SVAR.params(5),0,SVAR.params(8);
         SVAR.params(6),SVAR.params(7),0;
         0,0,SVAR.params(9)];

SVAR.Q3=[0,0,SVAR.params(12);
         SVAR.params(10),SVAR.params(11),SVAR.params(13);
         0,0,SVAR.params(14)];

SVAR.Q4=[SVAR.params(15),0,0;
        SVAR.params(16),SVAR.params(17),SVAR.params(18);
        0,0,0]; % new 0's in (3,1) and (1,2)

%recall Svar B is for the first period (GI)
SVAR.B2=SVAR.B+SVAR.Q2; %second period (GM)
SVAR.B3=SVAR.B2+SVAR.Q3; %third period (GR)
SVAR.B4=SVAR.B3+SVAR.Q4; %fourth period (COVID)

SVAR.V=pinv(SVAR.Hessian); % var-cov matrix
SVAR.SE_sums=SE_sums(SVAR); % s.e. of B, B2, B3 and B4

%sign normalization:
for i = 1:M
    if SVAR.B(i,i)<0
        SVAR.B(:,i)=-SVAR.B(:,i);
    end
end

for i = 1:M
    if SVAR.B2(i,i)<0
        SVAR.B2(:,i)=-SVAR.B2(:,i);
    end
end

for i = 1:M
    if SVAR.B3(i,i)<0
        SVAR.B3(:,i)=-SVAR.B3(:,i);
    end
end

for i = 1:M
    if SVAR.B4(i,i)<0
        SVAR.B4(:,i)=-SVAR.B4(:,i);
    end
end

%obtain Q2, Q3 and Q4 again after the sign normalization
SVAR.Q2=SVAR.B2-SVAR.B; 
SVAR.Q3=SVAR.B3-SVAR.B2;
SVAR.Q4=SVAR.B4-SVAR.B3;

%Rank condition check
V1=kron(SVAR.B,eye(M));
V2=kron(SVAR.B+ SVAR.Q2, eye(M));
V3=kron(SVAR.B+ SVAR.Q2+ SVAR.Q3, eye(M));
V4=kron(SVAR.B+ SVAR.Q2+ SVAR.Q3+ SVAR.Q4, eye(M));

% Calculates the matrix for checking the rank condition (full column rank)
RankMatrix=kron(eye(4),mDD)*[V1,   zeros(M^2,M^2),  zeros(M^2,M^2),    zeros(M^2,M^2);
                             V2,    V2,           zeros(M^2,M^2),    zeros(M^2,M^2);
                             V3,    V3,            V3,               zeros(M^2,M^2);   
                             V4,    V4,            V4,                  V4]; 
 
% Selection matrix for extracting the structural parameters                         
HSelection=zeros(4*M^2,size(SVAR.params,1));

HSelection(1,1)=1;
HSelection(2,2)=1;
HSelection(5,3)=1;
HSelection(9,4)=1;

HSelection(10,5)=1;
HSelection(11,6)=1;
HSelection(14,7)=1;
HSelection(16,8)=1;
HSelection(18,9)=1;

HSelection(20,10)=1;
HSelection(23,11)=1;
HSelection(25,12)=1;
HSelection(26,13)=1;
HSelection(27,14)=1;

HSelection(28,15)=1;
HSelection(29,16)=1;
HSelection(32,17)=1;
HSelection(35,18)=1;

Jacobian=RankMatrix*HSelection;
Jacobian=2*Jacobian; % multiplicative scalar 2, see Fanelli 2015 OBES
% Report the rank of the matrix for checking the identification
disp(['The rank of the Jacobian is ' num2str(rank(Jacobian)) ' and the number of parameters is ' num2str(size(SVAR.params,1))])  %full col rank (18)

SVAR_L=SVAR; %storing the results

%% SVAR (proposal)
disp(' ')
disp('Proposal (macro uncertainty exogenous and unidirectionality between macro and financial uncertainty up to covid period):')
initialvalues=[0.0112;-0.1377;0.7540;-0.0259;-0.0037;0.0470;-0.2927;0.0040;0.0535;-0.0318;0.1355;0.0014;-0.0839;-0.0066;0;0;0;0;0]; %initial values (the first 16 are the estimates obtain by ABCF, we set the new parameters initial value equal to 0)

options = optimset('MaxFunEvals',200000,'TolFun',1e-500,'MaxIter',200000,'TolX',1e-100,'Display', 'off'); % optimization settings
fun = @(teta) ll_Proposal_MS(teta,times,sigmas);
[SVAR.params,nLk,~,~,~,SVAR.Hessian] = fminunc(fun,initialvalues,options);

SVAR.Lk=-nLk; %invert the sign of the likelihood
SVAR.SE=diag(pinv(SVAR.Hessian)).^0.5; %standard errors from the Hessian
SVAR.t_stat=SVAR.params./SVAR.SE; %t-stat

PVarl_LRTest = 1-chi2cdf(2*(regimes{1}.Lk+regimes{2}.Lk+regimes{3}.Lk+regimes{4}.Lk-SVAR.Lk),24-size(SVAR.params,1)); % overidentification test
disp(['The p-value of the overidentification restrictions is ' num2str(PVarl_LRTest)]) 
%crucial result above: we cannot reject the model in which macro uncertainty is exogenous, but we have bidirectionality

SVAR.B=[SVAR.params(1),0,0;
        SVAR.params(2),SVAR.params(3),0;
        0,0,SVAR.params(4)];

SVAR.Q2=[SVAR.params(5),0,SVAR.params(8);
         SVAR.params(6),SVAR.params(7),0;
         0,0,SVAR.params(9)];

SVAR.Q3=[0,0,SVAR.params(12);
         SVAR.params(10),SVAR.params(11),SVAR.params(13);
         0,0,SVAR.params(14)];

SVAR.Q4=[SVAR.params(15),0,0;
        SVAR.params(16),SVAR.params(18),SVAR.params(19);
        SVAR.params(17),0,0];  % now (3,1) is a free element

%recall Svar B is for the first period (GI)
SVAR.B2=SVAR.B+SVAR.Q2; %second period (GM)
SVAR.B3=SVAR.B2+SVAR.Q3; %third period (GR)
SVAR.B4=SVAR.B3+SVAR.Q4; %fourth period (COVID)

SVAR.V=pinv(SVAR.Hessian); % var-cov matrix
SVAR.SE_sums=SE_sums(SVAR); % s.e. of B, B2, B3 and B4

%sign normalization:
for i = 1:M
    if SVAR.B(i,i)<0
        SVAR.B(:,i)=-SVAR.B(:,i);
    end
end

for i = 1:M
    if SVAR.B2(i,i)<0
        SVAR.B2(:,i)=-SVAR.B2(:,i);
    end
end

for i = 1:M
    if SVAR.B3(i,i)<0
        SVAR.B3(:,i)=-SVAR.B3(:,i);
    end
end

for i = 1:M
    if SVAR.B4(i,i)<0
        SVAR.B4(:,i)=-SVAR.B4(:,i);
    end
end

%obtain Q2, Q3 and Q4 again after the sign normalisation
SVAR.Q2=SVAR.B2-SVAR.B; 
SVAR.Q3=SVAR.B3-SVAR.B2;
SVAR.Q4=SVAR.B4-SVAR.B3;

%Rank condition check
V1=kron(SVAR.B,eye(M));
V2=kron(SVAR.B+ SVAR.Q2, eye(M));
V3=kron(SVAR.B+ SVAR.Q2+ SVAR.Q3, eye(M));
V4=kron(SVAR.B+ SVAR.Q2+ SVAR.Q3+ SVAR.Q4, eye(M));

% Calculates the matrix for checking the rank condition (full column rank)
RankMatrix=kron(eye(4),mDD)*[V1,   zeros(M^2,M^2),  zeros(M^2,M^2),    zeros(M^2,M^2);
                             V2,    V2,           zeros(M^2,M^2),    zeros(M^2,M^2);
                             V3,    V3,            V3,               zeros(M^2,M^2);   
                             V4,    V4,            V4,                  V4]; 
 
% Selection matrix for extracting the structural parameters                         
HSelection=zeros(4*M^2,size(SVAR.params,1));

HSelection(1,1)=1;
HSelection(2,2)=1;
HSelection(5,3)=1;
HSelection(9,4)=1;

HSelection(10,5)=1;
HSelection(11,6)=1;
HSelection(14,7)=1;
HSelection(16,8)=1;
HSelection(18,9)=1;

HSelection(20,10)=1;
HSelection(23,11)=1;
HSelection(25,12)=1;
HSelection(26,13)=1;
HSelection(27,14)=1;

HSelection(28,15)=1;
HSelection(29,16)=1;
HSelection(30,17)=1;
HSelection(32,18)=1;
HSelection(35,19)=1;

Jacobian=RankMatrix*HSelection;
Jacobian=2*Jacobian; % multiplicative scalar 2, see Fanelli 2015 OBES
% Report the rank of the matrix for checking the identification
disp(['The rank of the Jacobian is ' num2str(rank(Jacobian)) ' and the number of parameters is ' num2str(size(SVAR.params,1))])  %full col rank (20)

SVAR_P=SVAR; %storing the results

%% IRF
IRF=cell(1,4);
R=[eye(M) zeros(M,M*(lags-1))];  % selection matrix R used in IRF computation 
% VAR companion matrix 12x12 is kept in the structure regime estimated (i-th cell for i-th regime)

for k=1:4

    if k==1
        B=SVAR.B;
    elseif k==2
        B=SVAR.B2;
    elseif k==3
        B=SVAR.B3;
    else
        B=SVAR.B4;
    end % this if statment chooses the right B matrix for computing the IRF of the k-th period

    for h = 0 : HorizonIRF
        PHI(:,:,h+1)=R*regimes{k}.Companion_Matrix^h*R'*B;  %R C_i^h R' B_i
    
        for i=1:M
            for j=1:M
                IRF{k}(h+1,M*(i-1)+j)=PHI(i,j,h+1);
            end
        end % creating a matrix with M^2 columns and HorizonIRF rows, for plotting
    end
end


%% LL
% We use the unique VAR estimated above and whose results are stored in the
% structure 'All'. We use as Sigmas the matrices contained in the cell
% 'sigmas_test', that were the ones estimated under the assumption of no
% change in Pi (that we keep in this approach).
   
% SVAR (LL)
initialvalues=[
    0.0112;-0.1377;0;0.0003;0.7540;0;0;0;-0.0259;
    ones(9,1)]; % we choose as initial condition for B the same used for the other estimations (ABCF results)

options = optimset('MaxFunEvals',200000,'TolFun',1e-500,'MaxIter',200000,'TolX',1e-100,'Display', 'off'); % optimization settings
fun = @(teta) ll_LL(teta,times,sigmas_test);
[LL.params,nLk] = fminunc(fun,initialvalues,options);
LL.B=reshape(LL.params(1:9),3,3); % storing B

%sign normalization:
for i = 1:M
    if LL.B(i,i)<0
        LL.B(:,i)=-LL.B(:,i);
    end
end

% LL IRF
IRF_LL=[];
for h = 0 : HorizonIRF
    PHI(:,:,h+1)=R*All.Companion_Matrix^h*R'*LL.B;  %R C^h R' B
    for i=1:M
        for j=1:M
            IRF_LL(h+1,M*(i-1)+j)=PHI(i,j,h+1);
        end
    end
end

%% plots of 4 periods in 1 (+LL)
F0=figure('Position', get(0, 'Screensize'));
for k=1:9
    subplot(3,3,k)
    plot(0:60,IRF{1,1}(:,k),'Color','#0968AB','LineWidth',1.7) % IRF 1-st period
    hold on
    plot(0:60,IRF{1,2}(:,k),'Color','#D46737','LineWidth',1.7) % IRF 2-nd period
    plot(0:60,IRF{1,3}(:,k),'Color','#E9B530' ,'LineWidth',1.7) % IRF 3-rd period
    plot(0:60,IRF{1,4}(:,k),'Color','#479D50','LineWidth',1.7) % IRF 4-th period
    plot(0:60,IRF_LL(:,k),'--black','LineWidth',1.3) % IRF LL
    yline(0,'Color','black','LineWidth',0.5)
    hold off

    if k==1||k==4||k==7
       t=title('$U_{Mt}$ shock','Interpreter','latex');
    elseif k==2||k==5||k==8
       t=title('$Y_{t}$ shock','Interpreter','latex');
    elseif k==3||k==6||k==9
       t=title('$U_{Ft}$ shock','Interpreter','latex');
    end

    if k==1||k==2||k==3
       yl=ylabel('$U_{Mt}$','Interpreter','latex');
    elseif k==4||k==5||k==6
       yl=ylabel('$Y_{t}$','Interpreter','latex');
    elseif k==7||k==8||k==9
       yl=ylabel('$U_{Ft}$','Interpreter','latex');
    end
    t.FontSize=20;
    yl.FontSize=20;
    yl.Rotation=360;
end
saveas(F0,'images/IRF_allR_4regimes.png')

%% bootstrapped IRF
regimes{1}.res_boot=[];
regimes{2}.res_boot=[];
regimes{3}.res_boot=[];
regimes{4}.res_boot=[];
IRF_b=cell(1,4);
for boot = 1 : repetitions
    for i=1:4
        T=size(regimes{i}.res,1);
        TBoot=datasample(1:T,T); % resample from 1 to T
        regimes{i}.res_boot=[regimes{i}.res_boot regimes{i}.res(TBoot,:)];  % bootstrap errors 
    end
end % this loop generates the M residuals series 'repetitions' times, and put all of them in a T x M*repetitions matrix
Residuals_B=[regimes{1}.res_boot;regimes{2}.res_boot;regimes{3}.res_boot;regimes{4}.res_boot];

data_B=zeros(size(Residuals_B,1)+lags,M*repetitions);
data_B(1:lags,:)=kron(ones(1,repetitions),DataSet(1:lags,:)); % setting the first elements equal to the original data

for t=1+lags:TBs(1)
    data_h=zeros(1,size(data_B,2)); % empty row of size 1 x M*repetitions
    for p=1:lags
        data_h=data_h+data_B(t-p,:)*kron(eye(repetitions),regimes{1}.Pi{p}');
    end % loop that computes the predicted obs of t+1 (without the error and the deterministic components)
    data_B(t,:)=data_h+kron(ones(1,repetitions),regimes{1}.Const')+Residuals_B(t-lags,:); % adding the error and the deterministic components
    disp(t)
end % loop that creates the bootstrap M series 'repetitions' times and put them in a T x M*repetitions

for t=TBs(1)+1:TBs(2)
    data_h=zeros(1,size(data_B,2)); % empty row of size 1 x M*repetitions
    for p=1:lags
        data_h=data_h+data_B(t-p,:)*kron(eye(repetitions),regimes{2}.Pi{p}');
    end % loop that computes the predicted obs of t+1 (without the error and the deterministic components)
    data_B(t,:)=data_h+kron(ones(1,repetitions),regimes{2}.Const')+Residuals_B(t-lags,:); % adding the error and the deterministic components
    disp(t)
end 

for t=TBs(2)+1:TBs(3)
    data_h=zeros(1,size(data_B,2)); % empty row of size 1 x M*repetitions
    for p=1:lags
        data_h=data_h+data_B(t-p,:)*kron(eye(repetitions),regimes{3}.Pi{p}');
    end % loop that computes the predicted obs of t+1 (without the error and the deterministic components)
    data_B(t,:)=data_h+kron(ones(1,repetitions),regimes{3}.Const')+Residuals_B(t-lags,:); % adding the error and the deterministic components
    disp(t)
end 

for t=TBs(3)+1:size(DataSet,1)
    data_h=zeros(1,size(data_B,2)); % empty row of size 1 x M*repetitions
    for p=1:lags
        data_h=data_h+data_B(t-p,:)*kron(eye(repetitions),regimes{4}.Pi{p}');
    end % loop that computes the predicted obs of t+1 (without the error and the deterministic components)
    data_B(t,:)=data_h+kron(ones(1,repetitions),regimes{4}.Const')+Residuals_B(t-lags,:); % adding the error and the deterministic components
    disp(t)
end 

data_B_all=cell(1,repetitions); % empty cell
DataSets_b=cell(4,repetitions);
for boot = 1 : repetitions
    data_B_all{boot}=data_B(:,1+(boot-1)*M:M+(boot-1)*M); % cell with repetitions different datasets
    DataSets_b{1,boot}=data_B_all{boot}(1:TBs(1),:); % cell with the 4 different datasets
    DataSets_b{2,boot}=data_B_all{boot}(TBs(1)+1-lags:TBs(2),:);
    DataSets_b{3,boot}=data_B_all{boot}(TBs(2)+1-lags:TBs(3),:);
    DataSets_b{4,boot}=data_B_all{boot}(TBs(3)+1-lags:end,:); % dividing the datasets in the 4 volatility period
end 

Sigma_b=zeros(3,3,4,repetitions);
for boot= 1: repetitions
    disp(boot)
    regimes_b=cell(1,4); %4 volatility regimes
    for i=1:4
        [EstVAR,EstSE,regimes_b{i}.Lk,regimes_b{i}.res] = estimate(VAR,DataSets_b{i,boot});
        regimes_b{i}.T = size(DataSets{i},1)-lags;
        regimes_b{i}.Companion_Matrix=[cell2mat(EstVAR.AR); eye(M*(lags-1),M*(lags-1)) zeros(M*(lags-1),M)]; % companion matrix
        regimes_b{i}.Sigma=EstVAR.Covariance; %var-cov matrix
    end

    sigmas_b={regimes_b{1}.Sigma,regimes_b{2}.Sigma,regimes_b{3}.Sigma,regimes_b{4}.Sigma}; %sigma of each regime

    initialvalues=[0.0112;-0.1377;0.7540;-0.0259;-0.0037;0.0470;-0.2927;0.0040;0.0535;-0.0318;0.1355;0.0014;-0.0839;-0.0066;0;0;0;0;0]; %initial values (the first 16 are the estimates obtain by ABCF, we set the new parameters initial value equal to 0)
    options = optimset('MaxFunEvals',200000,'TolFun',1e-500,'MaxIter',200000,'TolX',1e-100,'Display', 'off'); % optimization settings
    fun = @(teta) ll_Proposal_MS(teta,times,sigmas_b);
    [SVAR_b.params,nLk,~,~,~,SVAR_b.Hessian] = fminunc(fun,initialvalues,options);

    SVAR_b.B=[SVAR_b.params(1),0,0;
            SVAR_b.params(2),SVAR_b.params(3),0;
            0,0,SVAR_b.params(4)];
    
    SVAR_b.Q2=[SVAR_b.params(5),0,SVAR_b.params(8);
             SVAR_b.params(6),SVAR_b.params(7),0;
             0,0,SVAR_b.params(9)];
    
    SVAR_b.Q3=[0,0,SVAR_b.params(12);
             SVAR_b.params(10),SVAR_b.params(11),SVAR_b.params(13);
             0,0,SVAR_b.params(14)];
    
    SVAR_b.Q4=[SVAR_b.params(15),0,0;
            SVAR_b.params(16),SVAR_b.params(18),SVAR_b.params(19);
            SVAR_b.params(17),0,0];

    %recall Svar B is for the first period (GI)
    SVAR_b.B2=SVAR_b.B+SVAR_b.Q2; %second period (GM)
    SVAR_b.B3=SVAR_b.B2+SVAR_b.Q3; %third period (GR)
    SVAR_b.B4=SVAR_b.B3+SVAR_b.Q4; %forth period (COVID)

    %sign normalization:
    for i = 1:M
        if SVAR_b.B(i,i)<0
            SVAR_b.B(:,i)=-SVAR_b.B(:,i);
        end
    end
    
    for i = 1:M
        if SVAR_b.B2(i,i)<0
            SVAR_b.B2(:,i)=-SVAR_b.B2(:,i);
        end
    end
    
    for i = 1:M
        if SVAR_b.B3(i,i)<0
            SVAR_b.B3(:,i)=-SVAR_b.B3(:,i);
        end
    end
    
    for i = 1:M
        if SVAR_b.B4(i,i)<0
            SVAR_b.B4(:,i)=-SVAR_b.B4(:,i);
        end
    end
    
    for k=1:4
        if k==1
            B=SVAR_b.B;
        elseif k==2
            B=SVAR_b.B2;
        elseif k==3
            B=SVAR_b.B3;
        else
            B=SVAR_b.B4;
        end % this if statment chooses the right B matrix for computing the IRF of the k-th period
    
        for h = 0 : HorizonIRF
            PHI(:,:,h+1)=R*regimes_b{k}.Companion_Matrix^h*R'*B;  %R C_i^h R' B_i
        
            for i=1:M
                for j=1:M
                    IRF_b{k,M*(i-1)+j}(h+1,boot)=PHI(i,j,h+1);
                end
            end % creating a matrix with M^2 columns and HorizonIRF rows, for plotting
        end
    end
end

%% plots of the single periods with c.i.
IRF_inf=cell(1,4);
IRF_sup=cell(1,4);
for j=1:4
    for i=1:9
        IRF_inf{j}(:,i)=prctile(IRF_b{j,i},5,2);
        IRF_sup{j}(:,i)=prctile(IRF_b{j,i},95,2);
    end
end

for i=1:4
    F=figure('Position', get(0, 'Screensize'));
    for k=1:9
        subplot(3,3,k)
        plot(0:60,IRF{i}(:,k),'Color','#0968AB','LineWidth',2) % IRF
        hold on
        fill([0:60, fliplr(0:60)], [IRF_inf{i}(:,k)', fliplr(IRF_sup{i}(:,k)')],[0.035,0.408,0.671],'FaceAlpha',0.2,'EdgeColor','none') % C.I. shaded area
        yline(0,'Color','black','LineWidth',0.5)
        hold off
    
        if k==1||k==4||k==7
           t=title('$U_{Mt}$ shock','Interpreter','latex');
        elseif k==2||k==5||k==8
           t=title('$Y_{t}$ shock','Interpreter','latex');
        elseif k==3||k==6||k==9
           t=title('$U_{Ft}$ shock','Interpreter','latex');
        end
    
        if k==1||k==2||k==3
           yl=ylabel('$U_{Mt}$','Interpreter','latex');
        elseif k==4||k==5||k==6
           yl=ylabel('$Y_{t}$','Interpreter','latex');
        elseif k==7||k==8||k==9
           yl=ylabel('$U_{Ft}$','Interpreter','latex');
        end
        t.FontSize=20;
        yl.FontSize=20;
        yl.Rotation=360;
        if i==1
            sgtitle('First volatility regime' ,'Interpreter','latex','FontSize',28)
        elseif i==2
            sgtitle('Second volatility regime','Interpreter','latex','FontSize',28)
        elseif i==3
            sgtitle('Third volatility regime','Interpreter','latex','FontSize',28)
        else
            sgtitle('Fourth volatility regime','Interpreter','latex','FontSize',28)
        end        
    end
    saveas(F,['images/IRF_' num2str(i) 'R_4regimes.png'])
end

%% pvalues

% SVAR (upper)
tstat.SVAR_U.B=SVAR_U.B./SVAR_U.SE_sums{1};
tstat.SVAR_U.B2=SVAR_U.B2./SVAR_U.SE_sums{2};
tstat.SVAR_U.B3=SVAR_U.B3./SVAR_U.SE_sums{3};
tstat.SVAR_U.B4=SVAR_U.B4./SVAR_U.SE_sums{4};
pval.SVAR_U.B=(1-normcdf(abs(tstat.SVAR_U.B)))*2;
pval.SVAR_U.B2=(1-normcdf(abs(tstat.SVAR_U.B2)))*2;
pval.SVAR_U.B3=(1-normcdf(abs(tstat.SVAR_U.B3)))*2;
pval.SVAR_U.B4=(1-normcdf(abs(tstat.SVAR_U.B4)))*2;

% SVAR (lower)
tstat.SVAR_L.B=SVAR_L.B./SVAR_L.SE_sums{1};
tstat.SVAR_L.B2=SVAR_L.B2./SVAR_L.SE_sums{2};
tstat.SVAR_L.B3=SVAR_L.B3./SVAR_L.SE_sums{3};
tstat.SVAR_L.B4=SVAR_L.B4./SVAR_L.SE_sums{4};
pval.SVAR_L.B=(1-normcdf(abs(tstat.SVAR_L.B)))*2;
pval.SVAR_L.B2=(1-normcdf(abs(tstat.SVAR_L.B2)))*2;
pval.SVAR_L.B3=(1-normcdf(abs(tstat.SVAR_L.B3)))*2;
pval.SVAR_L.B4=(1-normcdf(abs(tstat.SVAR_L.B4)))*2;

% SVAR (proposal)
tstat.SVAR_P.B=SVAR_P.B./SVAR_P.SE_sums{1};
tstat.SVAR_P.B2=SVAR_P.B2./SVAR_P.SE_sums{2};
tstat.SVAR_P.B3=SVAR_P.B3./SVAR_P.SE_sums{3};
tstat.SVAR_P.B4=SVAR_P.B4./SVAR_P.SE_sums{4};
pval.SVAR_P.B=(1-normcdf(abs(tstat.SVAR_P.B)))*2;
pval.SVAR_P.B2=(1-normcdf(abs(tstat.SVAR_P.B2)))*2;
pval.SVAR_P.B3=(1-normcdf(abs(tstat.SVAR_P.B3)))*2;
pval.SVAR_P.B4=(1-normcdf(abs(tstat.SVAR_P.B4)))*2;

% sigmas
for i=1:4
    tstat.Sigma{i}=regimes{i}.Sigma_t_stat;
    pval.Sigma{i}=(1-normcdf(abs(tstat.Sigma{i})))*2;
end
tstat.Sigma{5}=All.Sigma_t_stat;
pval.Sigma{5}=(1-normcdf(abs(tstat.Sigma{5})))*2;

% correlation matrices
for i=1:4
    A=((regimes{i}.T - 2)./(1 - regimes{i}.Sigma_Corr_Matrix.^2));
    A(1,1)=Inf;
    A(2,2)=Inf;
    A(3,3)=Inf;
    tstat.Corr{i}=regimes{i}.Sigma_Corr_Matrix.*A.^0.5;
    pval.Corr{i}=(1-normcdf(abs(tstat.Corr{i})))*2;
    pval.Corr{i}(1,1)=NaN;
    pval.Corr{i}(2,2)=NaN;
    pval.Corr{i}(3,3)=NaN;
end
A=((All.T - 2)./(1 - All.Sigma_Corr_Matrix.^2));
A(1,1)=Inf;
A(2,2)=Inf;
A(3,3)=Inf;
tstat.Corr{5}=All.Sigma_Corr_Matrix.*A.^0.5;
pval.Corr{5}=(1-normcdf(abs(tstat.Corr{5})))*2;
pval.Corr{5}(1,1)=NaN;
pval.Corr{5}(2,2)=NaN;
pval.Corr{5}(3,3)=NaN;

%% export
latex.SVAR_P.B=to_latex(SVAR_P.B,pval.SVAR_P.B);
latex.SVAR_P.B2=to_latex(SVAR_P.B2,pval.SVAR_P.B2);
latex.SVAR_P.B3=to_latex(SVAR_P.B3,pval.SVAR_P.B3);
latex.SVAR_P.B4=to_latex(SVAR_P.B4,pval.SVAR_P.B4);

latex.SVAR_U.B=to_latex(SVAR_U.B,pval.SVAR_U.B);
latex.SVAR_U.B2=to_latex(SVAR_U.B2,pval.SVAR_U.B2);
latex.SVAR_U.B3=to_latex(SVAR_U.B3,pval.SVAR_U.B3);
latex.SVAR_U.B4=to_latex(SVAR_U.B4,pval.SVAR_U.B4);

latex.SVAR_L.B=to_latex(SVAR_L.B,pval.SVAR_L.B);
latex.SVAR_L.B2=to_latex(SVAR_L.B2,pval.SVAR_L.B2);
latex.SVAR_L.B3=to_latex(SVAR_L.B3,pval.SVAR_L.B3);
latex.SVAR_L.B4=to_latex(SVAR_L.B4,pval.SVAR_L.B4);

for i=1:4
    latex.Sigma{i}=to_latex_sym(regimes{i}.Sigma,pval.Sigma{i});
end
latex.Sigma{5}=to_latex_sym(All.Sigma,pval.Sigma{5});

for i=1:4
    latex.Corr{i}=to_latex_sym(regimes{i}.Sigma_Corr_Matrix,pval.Corr{i});
end
latex.Corr{5}=to_latex_sym(All.Sigma_Corr_Matrix,pval.Corr{5});

for i=1:5
    disp(latex.Sigma{i})
end


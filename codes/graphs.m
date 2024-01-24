clear
clc
cd '/Users/samueleborsini/Library/Mobile Documents/com~apple~CloudDocs/Università/Economics and econometrics/II anno/Structural Macroeconometrics/Paper/final'
data=readtable("Data_Set.csv");

%% time series plots
F=figure(1);
plot(table2array(data(:,'date')),table2array(data(:,'IP')),'Color','#0968AB','LineWidth',1.5) % IP plot
saveas(F,'images/IP.png')

F=figure(1);
plot(table2array(data(:,'date')),table2array(data(:,'UM1')),'Color','#D46737','LineWidth',1.5) %UM1 plot
hold on
plot(table2array(data(:,'date')),table2array(data(:,'UF1')),'Color','#479D50','LineWidth',1.5) %UF1 plot
legend('UM1','UF1','')
hold off
saveas(F,'images/UF1_UM1.png')

%% estimation
lags=4; % number of lags of the VAR

X_all=[table2array(data(:,"UM1")),table2array(data(:,"IP")),table2array(data(:,"UF1"))]; % selected variables of interest

for i = 1:lags
    VAR_PI{i} = nan(3);
end
VAR_c=nan(3,1);
VAR = varm('Constant',VAR_c,'AR',VAR_PI); % VAR setup

sigma_plot_10y=[];
for t=1:(size(data,1)-120)
    X=X_all(t:(t+119),:);
    [EstVAR,EstSE,logLikVAR,Residuals]=estimate(VAR,X);
    T=size(X,1)-lags;
    Sigma=Residuals'*Residuals/T;
    sigma_plot_10y=[sigma_plot_10y;Sigma(1,:),Sigma(2,:),Sigma(3,:)];
    disp(t)
end % 10 years rolling window estimation of the innovations' covariance matrix

sigma_plot_15y=[];
for t=1:(size(data,1)-180)
    X=X_all(t:(t+179),:);
    [EstVAR,EstSE,logLikVAR,Residuals]=estimate(VAR,X);
    T=size(X,1)-lags;
    Sigma=Residuals'*Residuals/T;
    sigma_plot_15y=[sigma_plot_15y;Sigma(1,:),Sigma(2,:),Sigma(3,:)];
    disp(t)
end % 15 years rolling window estimation of the innovations' covariance matrix

sigma_plot_R=[];
for t=1:(size(data,1)-120)
    X=X_all(1:(t+119),:);
    [EstVAR,EstSE,logLikVAR,Residuals]=estimate(VAR,X);
    T=size(X,1)-lags;
    Sigma=Residuals'*Residuals/T;
    sigma_plot_R=[sigma_plot_R;Sigma(1,:),Sigma(2,:),Sigma(3,:)];
    disp(t)
end % recursive estimation of the innovations' covariance matrix

%% plots
date10y=table2array(data(121:end,"date")); % dates
date15y=table2array(data(181:end,"date"));
titles=cell(3,3); % cell with the titles
titles{1,1}='Var$(U_{Mt})$';
titles{2,1}='Cov$(U_{Mt},Y_{t})$';
titles{3,1}='Cov$(U_{Mt},U_{Ft})$';
titles{2,2}='Var$(Y_{t})$';
titles{3,2}='Cov$(Y_{t},U_{Ft})$';
titles{3,3}='Var$(U_{Ft})$';

F=figure('DefaultAxesFontSize',6);
for k=[1,2,3,5,6,9]
    subplot(3,3,k)
    plot(date10y,sigma_plot_10y(:,k),'Color','#D46737','LineWidth',1) % 10 years rolling window estimates
    hold on
    plot(date15y,sigma_plot_15y(:,k),'Color','#EDB120','LineWidth',1) % 15 years rolling window estimates
    plot(date10y,sigma_plot_R(:,k),'Color','#0968AB','LineWidth',1) % recursive estimates
    xline(datetime('1984-03-01','InputFormat','yyyy-MM-dd'),'--black','LineWidth',0.8) % 1-st break date: 03:1984
    xline(datetime('2007-12-01','InputFormat','yyyy-MM-dd'),'--black','LineWidth',0.8) % 2-nd break date: 12:2007
    xline(datetime('2020-02-01','InputFormat','yyyy-MM-dd'),'--black','LineWidth',0.8) % 3-rd break date: 03:2020
    hold off
    title(titles{k},'Interpreter','latex','FontSize',12)
end
saveas(F,'images/cov_matrix_est.png')
%==========================================================================
%                      DSGE MODEL: FORECASTS
%
%
% Author: Luigi Bocola         lbocola@sas.upenn.edu
% Date  : 06/14/2013
%==========================================================================

%=========================================================================
%                              HOUSEKEEPING
%=========================================================================

clc
clear all
close all
delete *.asv

tic

l = path;

path('Mfiles',path);
path('Optimization Routines',path);
path('LRE',path);
path('Matfiles',path);


disp('                                                                  ');
disp('   BAYESIAN ESTIMATION OF DSGE MODEL: FORECASTS      ');
disp('                                                                  ');

%=========================================================================
%    LOAD POSTERIOR DRAWS AND COMPUTE FORECASTS FOR OUTPUT AND HOURS
%=========================================================================

load mhdraws            % Load MH draws

[Nsim,Npam]     = size(Thetasim);

N_hor           = 30;   % Number of quarters for Variance Decomposition

nend            = 2;    % Number of Endogenous Variables

forecast_model  = zeros(N_hor,Nsim,nend);

counter         = 0;

for i=1:Nsim
    
    para      = Thetasim(i,:);
    
    para(3)   = exp((1/4)*para(3));

    load data

    data(:,1) = data(:,1)-log(para(3))*linspace(1,size(data,1),size(data,1))';
    
[T1, TC, T0, TETA, RC, retcode] = model_solution(para);

    [A,B,H,R,Se,Phi] = sysmat(T1,T0,para);

    [liki,measurepredi,statepredi,varstatepredi] = kalman(A,B,H,R,Se,Phi,data);

    [U,S,V] = svd(squeeze(varstatepredi(end,:,:)));
    
    sigmult = U*sqrt(S);
    
    Sigma   = sigmult*sigmult';
    
    state   = statepredi(end,:)' + U*(mvnrnd(zeros(1,size(S,1)),Sigma))';
    
 for m=1:N_hor
    
    state                 = Phi*state;
        
    forecast_model(m,i,:) = A + [log(para(3));0]*(204+m) +B*state;

 end
    
if counter==100
disp('                                                                  ');
disp(['                               DRAW NUMBER:', num2str(i)]         );
counter = 0;
end

    counter   = counter+1;
    
end  

%=========================================================================
%              PLOT FORECAST (Median and 90% credible set)
%=========================================================================

load data

time   = (1955:.25:2008.75)';

pnames = strvcat('Output','Hours worked');

figure('Position',[20,20,900,600],'Color','w')

subplot(1,2,1),plot(time(181:end-12),data(181:end,1),'LineWidth',3), hold on
plot(time(end-11:end),median(squeeze(forecast_model(1:12,:,1)),2),'--r','LineWidth',3), hold on
jbfill(time(end-12:end)',[data(end,1);prctile(squeeze(forecast_model(1:12,:,1)),5,2)]',[data(end,1);prctile(squeeze(forecast_model(1:12,:,1)),95,2)]', rgb('grey'),rgb('grey'),0.4,0.5), hold on
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(1,:),'FontSize',18,'FontWeight','bold');
box off   

subplot(1,2,2),plot(time(181:end-12),data(181:end,2),'LineWidth',3), hold on
plot(time(end-11:end),median(squeeze(forecast_model(1:12,:,2)),2),'--r','LineWidth',3), hold on
jbfill(time(end-12:end)',[data(end,2);prctile(squeeze(forecast_model(1:12,:,2)),5,2)]',[data(end,2);prctile(squeeze(forecast_model(1:12,:,2)),95,2)]', rgb('grey'),rgb('grey'),0.4,0.5), hold on
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(2,:),'FontSize',18,'FontWeight','bold');
box off   

    
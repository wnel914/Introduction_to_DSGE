
%==========================================================================
%                      DSGE MODEL: FORECASTS
%
%
% Author: Luigi Bocola         lbocola@sas.upenn.edu
% Date  : 06/14/2013
%
% Updated: Jacob Warren        jacobwar@sas.upenn.edu
% Date   : 04/04/2016
%==========================================================================


%=========================================================================
%                              HOUSEKEEPING
%=========================================================================

clc; clear all; close all

l = path;

path('Mfiles',path);
path('Optimization Routines',path);
path('LRE',path);
path('Matfiles',path);
path('Extra',path);
picPath = [pwd, '\Figures\'];


disp('                                                                  ');
disp('   BAYESIAN ESTIMATION OF DSGE MODEL: FORECASTS      ');
disp('                                                                  ');

%=========================================================================
%    LOAD POSTERIOR DRAWS AND COMPUTE FORECASTS FOR OUTPUT AND HOURS
%=========================================================================

load mhdraws            % Load MH draws

[Nsim,Npam]     = size(Thetasim);

N_hor           = 30;   % Number of quarters for Variance Decomposition

nend            = 3;    % Number of Endogenous Variables

forecast_model  = zeros(N_hor,Nsim,nend);

counter         = 0;

for i=1:Nsim
    
    para      = Thetasim(i,:);
    
    [T1, TC, T0, TETA, RC, retcode] = model_solution(para);
    
    [A,B,H,R,Se,Phi] = sysmat(T1,T0,para);
    
    load us.txt
    
    data = us;
    
    [liki,measurepredi,statepredi,varstatepredi] = kalman(A,B,H,R,Se,Phi,data);
    
    [U,S,V] = svd(squeeze(varstatepredi(end,:,:)));
    
    sigmult = U*sqrt(S);
    
    Sigma   = sigmult*sigmult';
    
    state   = statepredi(end,:)' + U*(mvnrnd(zeros(1,size(S,1)),Sigma))';
    
    for m=1:N_hor
        
        state                 = Phi*state;
        
        forecast_model(m,i,:) = A  +B*state; 
        
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


time = (1983:.25:2005.75)';

pnames = strvcat('Output Growth','Inflation', 'Fed Funds');

figure('Position',[20,20,1200,600],'Color','w')

subplot(1,3,1),plot(time(length(data)-10:end-12),data(length(data)-10:end,1),'LineWidth',3), hold on
plot(time(end-11:end),median(squeeze(forecast_model(1:12,:,1)),2),'--r','LineWidth',3), hold on
jbfill(time(end-12:end)',[data(end,1);prctile(squeeze(forecast_model(1:12,:,1)),5,2)]',[data(end,1);prctile(squeeze(forecast_model(1:12,:,1)),95,2)]', rgb('grey'),rgb('grey'),0.4,0.5), hold on
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(1,:),'FontSize',18,'FontWeight','bold');
box off

subplot(1,3,2),plot(time(length(data)-10:end-12),data(length(data)-10:end,2),'LineWidth',3), hold on
plot(time(end-11:end),median(squeeze(forecast_model(1:12,:,2)),2),'--r','LineWidth',3), hold on
jbfill(time(end-12:end)',[data(end,2);prctile(squeeze(forecast_model(1:12,:,2)),5,2)]',[data(end,2);prctile(squeeze(forecast_model(1:12,:,2)),95,2)]', rgb('grey'),rgb('grey'),0.4,0.5), hold on
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(2,:),'FontSize',18,'FontWeight','bold');
box off

subplot(1,3,3),plot(time(length(data)-10:end-12),data(length(data)-10:end,3),'LineWidth',3), hold on
plot(time(end-11:end),median(squeeze(forecast_model(1:12,:,3)),2),'--r','LineWidth',3), hold on
jbfill(time(end-12:end)',[data(end,3);prctile(squeeze(forecast_model(1:12,:,3)),5,2)]',[data(end,3);prctile(squeeze(forecast_model(1:12,:,3)),95,2)]', rgb('grey'),rgb('grey'),0.4,0.5), hold on
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(3,:),'FontSize',18,'FontWeight','bold');
box off

print('-dpng', [picPath, 'Forecasts'])

%==========================================================================
%                  DSGE MODEL: POSTERIOR PREDICTIVE CHECKS
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
picPath = [pwd, '\Figures\'];


disp('                                                                  ');
disp('  BAYESIAN ESTIMATION OF DSGE MODEL: POSTERIOR PREDICTIVE CHECKS  ');
disp('                                                                  ');

%=========================================================================
% LOAD POSTERIOR DRAWS AND COMPUTE MODEL IMPLIED MOMENTS BY SIMULATION
%=========================================================================


load mhdraws

[Nsim,Npam] = size(Thetasim);

M           = 100;  % Number of periods for simulation

nend        = 3;    % Number of Endogenous Variables

nex         = 2;    % Number of Shocks

data_model  = zeros(M,nend,Nsim);

counter     = 1;

for i=1:Nsim
    
    para      = Thetasim(i,:);
    
    
    [T1, TC, T0, TETA, RC, retcode] = model_solution(para);
    
    [A,B,H,R,Se,Phi] = sysmat(T1,T0,para);
    
    e            = (Se.^(1/2))*randn(nend,M);
    
    state        = zeros(size(Phi,1),1);
    
    for m=1:M
        
        state              = Phi*state + R*e(:,m);
        
        data_model(m,:,i)  = A + B*state;
        
    end
    
    if counter==100
        disp('                                                                  ');
        disp(['                               DRAW NUMBER:', num2str(i)]         );
        counter = 0;
    end
    
    counter   = counter+1;
    
end

%=========================================================================
%                  PLOT POSTERIOR PREDICTIVE CHECKS
%=========================================================================

data = load('us.txt');

pnames = strvcat('Mean (GDP Growth)','Mean (Inflation)', 'Mean (Fed Funds)',...
                 'StDev (GDP Growth)' ,'StDev (Inflation)', 'Stdev (Fed Funds)');
             
figure('Position',[20,20,900,600],'Color','w')

%mean of all three series

[F,x] = ksdensity(squeeze(mean(data_model(:,1,:),1)));
subplot(2,3,1), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(mean(data(:,1))*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(1,:),'FontSize',16,'FontWeight','bold');
box off   

[F,x] = ksdensity(squeeze(mean(data_model(:,2,:),1)));
subplot(2,3,2), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(mean(data(:,2))*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(2,:),'FontSize',16,'FontWeight','bold');
box off   

[F,x] = ksdensity(squeeze(mean(data_model(:,3,:),1)));
subplot(2,3,3), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(mean(data(:,3))*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(3,:),'FontSize',16,'FontWeight','bold');
box off   


%StDev
[F,x] = ksdensity(squeeze(sqrt(var(data_model(:,1,:),1))) );
subplot(2,3,4), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(mean(data(:,1))*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(4,:),'FontSize',16,'FontWeight','bold');
box off   

[F,x] = ksdensity(squeeze(sqrt(var(data_model(:,2,:),1))) );
subplot(2,3,5), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(mean(data(:,2))*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(5,:),'FontSize',16,'FontWeight','bold');
box off   

[F,x] = ksdensity(squeeze(sqrt(var(data_model(:,3,:),1))) );
subplot(2,3,6), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(mean(data(:,3))*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(6,:),'FontSize',16,'FontWeight','bold');
box off   

print('-dpng', [picPath, 'PosteriorCheckMeanSD'])


%pairwise correlations:
corr_model = zeros(Nsim, 3);    %3 pair-wise correlations

for i = 1:Nsim
    
    data_i          = squeeze( data_model(:,:,i) );
    
    cormat          = corr(data_i);
    
    upperTri        = triu(cormat, 1);
    
    correlations    = upperTri(upperTri ~= 0);
    
    corr_model(i,:) = correlations;
    
end

pnames = strvcat('Corr(GDP Growth, Inflation)', 'Corr(GDP Growth, Inflation',...
                'Corr(Inflation, Fed Funds)'); 

         
figure('Position',[20,20,900,600],'Color','w')


[F,x] = ksdensity(corr_model(:,1));
subplot(1,3,1), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(corr(data(:,1), data(:,2))*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(1,:),'FontSize',16,'FontWeight','bold');
box off   


[F,x] = ksdensity(corr_model(:,2));
subplot(1,3,2), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(corr(data(:,1), data(:,3))*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(2,:),'FontSize',16,'FontWeight','bold');
box off   


[F,x] = ksdensity(corr_model(:,3));
subplot(1,3,3), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(corr(data(:,3), data(:,2))*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(3,:),'FontSize',16,'FontWeight','bold');
box off   

print('-dpng', [picPath, 'PosteriorCheckcorr'])

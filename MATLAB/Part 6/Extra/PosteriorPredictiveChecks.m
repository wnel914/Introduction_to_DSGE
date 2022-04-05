%==========================================================================
%                  DSGE MODEL: POSTERIOR PREDICTIVE CHECKS
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
disp('  BAYESIAN ESTIMATION OF DSGE MODEL: POSTERIOR PREDICTIVE CHECKS  ');
disp('                                                                  ');

%=========================================================================
% LOAD POSTERIOR DRAWS AND COMPUTE MODEL IMPLIED MOMENTS BY SIMULATION
%=========================================================================

load MH

load mhdraws

[Nsim,Npam] = size(Thetasim);

M           = 100; % Number of periods for simulation 

nend        = 2;    % Number of Endogenous Variables

nex         = 2;    % Number of Shocks

data_model  = zeros(M,nend,Nsim);

counter     = 1;

for i=1:Nsim
        
    para      = Thetasim(i,:);
    
    para(3)   = exp((1/4)*para(3));

[T1, TC, T0, TETA, RC, retcode] = model_solution(para);

    [A,B,H,R,Se,Phi] = sysmat(T1,T0,para);

    e            = (Se.^(1/2))*randn(2,M);
    
    state        = zeros(size(Phi,1),1);
    
    for m=1:M
    
    state              = Phi*state + R*e(:,m);
        
    data_model(m,:,i)  = A + [log(para(3));0]*m + B*state;

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

load data

pnames = strvcat('Mean (GDP Growth)','StDev (Gdp Growth)', 'StDev (Hours)','Corr(GDP Growth, Hours)');

figure('Position',[20,20,900,600],'Color','w')

[F,x] = ksdensity(400*squeeze(mean(data_model(2:end,1,:)-data_model(1:end-1,1,:),1)));
subplot(2,2,1), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(400*mean(data(2:end,1)-data(1:end-1,1))*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(1,:),'FontSize',16,'FontWeight','bold');
box off   

[F,x] = ksdensity(400*squeeze(var(data_model(2:end,1,:)-data_model(1:end-1,1,:),1).^(1/2)));
subplot(2,2,2), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(400*var(data(2:end,1)-data(1:end-1,1)).^(1/2)*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(2,:),'FontSize',16,'FontWeight','bold');
box off   

[F,x] = ksdensity(squeeze(var(data_model(1:end,2,:),1).^(1/2)));
subplot(2,2,3), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(var(data(1:end,2)).^(1/2)*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(3,:),'FontSize',16,'FontWeight','bold');
box off   

corr_model        = zeros(Nsim,1);

for i=1:Nsim

    gdp_growth    = data_model(2:end,1,i)-data_model(1:end-1,1,i);
    
    hours         = data_model(2:end,2,i);
    
    corr_model(i) = corr(gdp_growth,hours);
    
end
    
[F,x] = ksdensity(corr_model);
subplot(2,2,4), plot(x,F,'-b','LineWidth',4), hold on
a = max(F);
plot(corr(data(2:end,1)-data(1:end-1,1),data(2:end,2))*ones(100,1),linspace(0,a,100),'--r','LineWidth',4);
set(gca,'FontSize',18,'FontWeight','bold');
title(pnames(4,:),'FontSize',16,'FontWeight','bold');
box off   



    
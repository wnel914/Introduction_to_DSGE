%==========================================================================
%                          DSGE MODEL IRFs
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
disp('  BAYESIAN ESTIMATION OF DSGE MODEL: IMPULSE RESPONSE FUNCTIONS   ');
disp('                                                                  ');

%=========================================================================
%              LOAD POSTERIOR DRAWS AND COMPUTE IRFS
%=========================================================================

load mhdraws

[Nsim,Npam] = size(Thetasim);

Nirf        = 30; % Number of Quarters for IRFs

nend        = 3;  % Number of Endogenous Variables

IRF_tech    = zeros(Nirf,Nsim,nend);

IRF_pref    = zeros(Nirf,Nsim,nend);

counter     = 0;

for i=1:Nsim
    
    para         = Thetasim(i,:);
    
    para(3)      = exp((1/4)*para(3));
    
    [T1, TC, T0, TETA, RC, retcode] = model_solution(para);
    
    [A,B1,H,R,Se,Phi] = sysmat(T1,T0,para);
    
%    B            = [B1;[0,0,1,zeros(1,9),1];[zeros(1,6),1,zeros(1,5),1];... 
        % Select Output, Hours, Consumption, Investment, Wages
 %       [0,1,zeros(1,10),1];[1,zeros(1,12)]];


 
    e            = zeros(3,Nirf);
    
    e(1,1)       = 1;
    
    e            = (Se.^(1/2))*e;
    
    state        = zeros(size(Phi,1),1);
    
    for m=1:Nirf
        
        state              = Phi*state + R*e(:,m);
        
        IRF_tech(m,i,:)    = B1*state;
        
    end
    
    e            = zeros(3,Nirf);
    
    e(2,1)       = 1;
    
    e            = (Se.^(1/2))*e;
    
    state        = zeros(size(Phi,1),1);
    
    for m=1:Nirf
        
        state              = Phi*state + R*e(:,m);
        
        IRF_pref(m,i,:)    = B1*state;
        
    end
    
    if counter==100
        disp('                                                                  ');
        disp(['                               DRAW NUMBER:', num2str(i)]         );
        counter = 0;
    end
    
    counter   = counter+1;
    
end

%=========================================================================
%                  PLOT IRFS (Median and 90% credible set)
%=========================================================================

pnames = strvcat('Output','Hours', 'Consumption','Investment','Wages',...
    'Interest Rate');

pnames = strvcat('Output','Inflation', 'Interest Rates');

figure('Position',[20,20,900,600],'Color','w')

for j=1:3
    
    subplot(3,1,j), plot(median(IRF_tech(:,:,j),2)*100,'LineWidth',3), hold on
    plot(prctile(IRF_tech(:,:,j),5,2)*100,'--r','LineWidth',1), hold on
    plot(prctile(IRF_tech(:,:,j),95,2)*100,'--r','LineWidth',1), hold on
    set(gca,'FontSize',18,'FontWeight','bold');
    title(pnames(j,:),'FontSize',16,'FontWeight','bold');
    box off
    
end

[ax1,h3]=suplabel('IRFs to Technology Shocks','t');
set(h3,'FontSize',18,'FontWeight','bold')


figure('Position',[20,20,900,600],'Color','w')

for j=1:3
    
    subplot(3,1,j), plot(median(IRF_pref(:,:,j),2)*100,'LineWidth',3), hold on
    plot(prctile(IRF_pref(:,:,j),5,2)*100,'--r','LineWidth',1), hold on
    plot(prctile(IRF_pref(:,:,j),95,2)*100,'--r','LineWidth',1), hold on
    set(gca,'FontSize',18,'FontWeight','bold');
    title(pnames(j,:),'FontSize',16,'FontWeight','bold');
    box off
    
end

[ax1,h3]=suplabel('IRFs to Preference Shocks','t');
set(h3,'FontSize',18,'FontWeight','bold')

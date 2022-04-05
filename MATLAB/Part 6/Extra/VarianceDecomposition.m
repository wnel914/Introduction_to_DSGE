%==========================================================================
%                  DSGE MODEL: VARIANCE DECOMPOSITION
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
disp('   BAYESIAN ESTIMATION OF DSGE MODEL: VARIANCE DECOMPOSITION      ');
disp('                                                                  ');

%=========================================================================
%        LOAD POSTERIOR DRAWS AND COMPUTE VARIANCE DECOMPOSITION
%=========================================================================

load mhdraws    % Load MH draws

[Nsim,Npam]     = size(Thetasim);

N_hor           = 30;   % Number of quarters for Variance Decomposition

nend            = 2;    % Number of Endogenous Variables

variance        = zeros(nend,N_hor,Nsim);

variance_nopref = zeros(nend,N_hor,Nsim);

variance_notech = zeros(nend,N_hor,Nsim);

counter         = 0;

for i=1:Nsim
    
    para      = Thetasim(i,:);
    
    para(3)   = exp((1/4)*para(3));

[T1, TC, T0, TETA, RC, retcode] = model_solution(para);

    [A,B,H,R,Se,Phi]    = sysmat(T1,T0,para);

    Se_nopref           = Se;
    
    Se_nopref(2,2)      = 0;

for k=1:N_hor    
       
 variance(:,k,i)         = diag(B*(Phi^(k-1))*(R*Se*R')*(Phi^(k-1))'*B');
 
 variance_nopref(:,k,i)  = diag(B*(Phi^(k-1))*(R*Se_nopref*R')*(Phi^(k-1))'*B');

 variance_notech(:,k,i)  = variance(:,k,i)-variance_nopref(:,k,i);
  
end

if counter==100
disp('                                                                  ');
disp(['                               DRAW NUMBER:', num2str(i)]         );
counter = 0;
end

    counter   = counter+1;
    
end  

%=========================================================================
% PLOT VARIANCE DECOMPOSITION (Median and 90% credible set)
%=========================================================================
    
pnames = strvcat('Output (Technology)','Hours worked (Technology)',...
                 'Output (Preference)','Hours worked (Preference)');

figure('Position',[20,20,900,600],'Color','w')
    

subplot(2,2,1), plot(100*median((squeeze(variance_nopref(1,:,:))./squeeze(variance(1,:,:))),2),'b','LineWidth',3), hold on
                plot(100*prctile((squeeze(variance_nopref(1,:,:))./squeeze(variance(1,:,:))),5,2),'--r','LineWidth',1), hold on
                plot(100*prctile((squeeze(variance_nopref(1,:,:))./squeeze(variance(1,:,:))),95,2),'--r','LineWidth',1), hold on                
                set(gca,'FontSize',18,'FontWeight','bold');
                xlabel('Horizon','FontSize',14,'FontWeight','bold');
                title(pnames(1,:),'FontSize',16,'FontWeight','bold');
                box off   
                
subplot(2,2,2), plot(100*median((squeeze(variance_nopref(2,:,:))./squeeze(variance(2,:,:))),2),'b','LineWidth',3), hold on
                plot(100*prctile((squeeze(variance_nopref(2,:,:))./squeeze(variance(2,:,:))),5,2),'--r','LineWidth',1), hold on
                plot(100*prctile((squeeze(variance_nopref(2,:,:))./squeeze(variance(2,:,:))),95,2),'--r','LineWidth',1), hold on
                xlabel('Horizon','FontSize',14,'FontWeight','bold');               
                set(gca,'FontSize',18,'FontWeight','bold');
                title(pnames(2,:),'FontSize',16,'FontWeight','bold');
                box off   
                
subplot(2,2,3), plot(100*median((squeeze(variance_notech(1,:,:))./squeeze(variance(1,:,:))),2),'b','LineWidth',3), hold on
                plot(100*prctile((squeeze(variance_notech(1,:,:))./squeeze(variance(1,:,:))),5,2),'--r','LineWidth',1), hold on
                plot(100*prctile((squeeze(variance_notech(1,:,:))./squeeze(variance(1,:,:))),95,2),'--r','LineWidth',1), hold on  
                xlabel('Horizon','FontSize',14,'FontWeight','bold');              
                set(gca,'FontSize',18,'FontWeight','bold');
                title(pnames(3,:),'FontSize',16,'FontWeight','bold');
                box off   
                
subplot(2,2,4), plot(100*median((squeeze(variance_notech(2,:,:))./squeeze(variance(2,:,:))),2),'b','LineWidth',3), hold on
                plot(100*prctile((squeeze(variance_notech(2,:,:))./squeeze(variance(2,:,:))),5,2),'--r','LineWidth',1), hold on
                plot(100*prctile((squeeze(variance_notech(2,:,:))./squeeze(variance(2,:,:))),95,2),'--r','LineWidth',1), hold on     
                xlabel('Horizon','FontSize',14,'FontWeight','bold');           
                set(gca,'FontSize',18,'FontWeight','bold');
                title(pnames(4,:),'FontSize',16,'FontWeight','bold');
                box off   
    
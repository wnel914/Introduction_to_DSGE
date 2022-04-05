function [draw] = draw_prior

% The priors considered are:
%  1.) 100*(1/Beta-1) is GAMMA with mean 0.5 and std 0.5
%  2.) 100*log(Gamma) is NORMAL with mean .75 and std 0.5
%  3.) Lambda is GAMMA with mean .2 and std 0.2
%  4.) 100*log(PiStar) is GAMMA with mean 1 and std 0.5
%  5.) zetaP is BETA with mean 0.7 and std 0.15
%  6.) 1/(1+nu) is GAMMA with mean 1.50 and std 0.75
%  7.) rhoPhi is UNIFORM with mean 0 and std 1
%  8.) rhoLambda is UNIFORM with mean 0 and std 1
%  9.) rhoZ is UNIFORM with mean 0 and std 1
%  10.) 100*sigmaPhi is InvGamma with mean 2 and std 4.0
%  11.) 100*sigmaLambda is InvGamma with mean 0.5 and std 4.0
%  12.) 100*sigmaZ is InvGamma with mean 2 and std 4.0
%  13.) 100*sigmaR is InvGamma with mean 0.5 and std 4

% % % % NOTE: mean refers to para1, and std refers to para2

draw = ones(13,1);

P = ones(13,1);

% % GAMMA pdf
para1 = [0.5,  0.2, 1, 1.50];     %mean of 100(1/B-1), lambda, 100ln(phi*), 1/(1+vi)
para2 = [0.5,  0.2,  0.5, 0.75];     %var of 
b = para2.^2./para1;
a = para1./b;
draw([1,3,4,6]) = gamrnd(a,b);

% % NORMAL pdf 
para1 = 0.75;     %mean of 100log(gamma)
para2 = 0.5;     %var of 
draw(2) = normrnd(para1, para2);

% % BETA pdf
para1 = 0.7;     %mean of zetaP
para2 = 0.15;     %var of 
a = (1-para1).*para1.^2./para2.^2 - para1;
b = a.*(1./para1 - 1);
draw(5) = betarnd(a,b);

% % UNIFORM pdf
draw(7:9) = unifrnd(0,1,3,1);     %rho_{phi, lambda, z} 


%Inverse gamma
para1 = [2, 0.5, 2, 0.5]';     %100sigma(phi, lambda, z, r)
para2 = [4, 4, 4, 4]';
alpha = para2./2;
beta = para2./2 .* para1.^2/2;
draw(10:13) = 1./gamrnd(alpha, 1./beta);


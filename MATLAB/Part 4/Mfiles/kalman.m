function [liki] = kalman(A,B,H,R,Se,Phi,y)

%==========================================================================
%                               Kalman Filter                     
%                                                                 
%   This Function Implement the Kalman filter for the state space model:
%
%                        y = A + Bs +  u         u N(0,H)
%
%                        s = Phi*s(-1) + e   e N(0,Se)
%   
%  In the above system, y is an (l*1) vector of observable variables and 
%  s is an n*1 vector of latent variables. 
%
%  The Input of this function are the state space matrices and a T*l vector 
%  of observations for y.  
%  
%  The Output is:
%   - measurepredi = the one step ahead prediction for y;
%   - liki = is a T*1 vector containing p(yt|Yt-1). The first entry is based  
%            on the prediction of the state vector at its unconditional
%            mean;
%   - statepredi = It is the prediction E[st|Y^(t)];
%   - varstatepredi = It is the MSE of the above prediction;
%
%
%==========================================================================

% Initialize the State Vector at the Stationary Distribution

[T,l]    = size(y);
[n,n]    = size(Phi);
s        = zeros(T+1,n);
P        = zeros(T+1,n,n);
s(1,:)   = zeros(n,1)';

A        = A';
a        = inv(eye(n*n) - kron(Phi,Phi))*reshape(Se,n*n,1);
P(1,:,:) = reshape(a,n,n); 

% Kalman Filter Recursion
sprime             = zeros(n,1);
Pprime             = zeros(n,n);
liki               = zeros(T,1);

for i=1:T
    
  % Forecasting s(t)
  sprime = Phi*s(i,:)';
  Pprime = Phi*squeeze(P(i,:,:))*Phi' + Se;

  % Forecasting y(t)
  yprediction = A + B*sprime;
  v = y(i,:)' - yprediction;
  F = B*Pprime*B' + H;
  
  liki(i) = -0.5*l*log(2*pi) - 0.5*log(det(F)) - 0.5*v'*inv(F)*v;

  % Updating Step  
  kgain  = Pprime*B'*inv(F);
  s(i+1,:) = (sprime + kgain*v)';
  P(i+1,:,:) = Pprime - kgain*B*Pprime;

end

 
end

%==========================================================================
%                  QAR(1,1) MODEL: MARGINAL DATA DENSITY                   
%
%
%
% Author: Luigi Bocola     lbocola@sas.upenn.edu
% Date  : 10/09/2011
%==========================================================================

load Matfiles/mhdraws

%==========================================================================
% default parameters
%==========================================================================

block1   = 1;		% starting block
blockn   = 1;		% ending block

%load file
parasim = Thetasim;
likesim = logposterior;
index   = likesim==min(likesim);
parasim = parasim(index==0,:);
likesim = likesim(index==0,:);
    
% for data density (modified harmonic mean)
densfac = min(likesim)*0.9;
hmax    = 20;

%==========================================================================
% means and standard deviations
%==========================================================================

% initialize output
ndraws     = 0;
drawmean   = 0;
drawsqmean = 0;
sumdrawsq  = 0;

% loop for block
for jblock=block1:blockn
    
    % number of simulations in each block
	nsimul = size(parasim,1);
    npara  = size(parasim,2);
	% collect simulations
	drawblock = parasim;
	drawdim   = size(drawblock,2);

	% compute sums of x(i) and x(i)^2 to be used to calculate means and s.d.
	drawmean   = drawmean + sum(drawblock,1);
	drawsqmean = drawsqmean + sum(drawblock.^2,1);
	sumdrawsq  = sumdrawsq + parasim'*parasim;
	ndraws     = ndraws + nsimul;

end

% compute means and standard deviations
drawmean   = drawmean/ndraws;
drawsqmean = drawsqmean/ndraws;
drawstdd   = sqrt(drawsqmean - drawmean.^2);
drawsig    = sumdrawsq/ndraws - drawmean'*drawmean;

% check diagonal elements 
for j=1:size(drawsig,1)
	if drawsig(j,j) < 1E-6
		drawsig(j,j) = 1E-6;
	end
end

drawsiginv = inv(drawsig);
drawsiglndet = log(det(drawsig));

%==========================================================================
% marginal data density:  modified harmonic mean by Geweke (1999)
%==========================================================================
disp('                                                                ');
disp('                                                                ');

p = (.1:.1:.9)';
pcrit = chi2inv(p,ones(length(p),1)*npara);

ndraws  = 0;
suminvlike = zeros(length(p),1);
laginvlike = zeros(hmax,length(p));
gaminvlike = zeros(hmax,length(p));

% loop for block
for jblock=block1:blockn

	% number of simulations in each block
	[nsimul,npara] = size(parasim);				

	paradev  = parasim-repmat(drawmean,nsimul,1);
	quadpara = sum((paradev*drawsiginv).*paradev,2);

	% simulation loop
	for j=1:nsimul
		lnfpara = - 0.5*npara*log(2*pi) - 0.5*drawsiglndet - 0.5*quadpara(j) - log(p);
		indpara = (quadpara(j)<pcrit);
		invlike = exp(lnfpara - likesim(j) + densfac).*indpara;

		laginvlike = [invlike' ; laginvlike(1:hmax-1,:)];
		gaminvlike = gaminvlike + laginvlike.*repmat(invlike',hmax,1);
		suminvlike = suminvlike + invlike;
	end	% simulation loop

	ndraws = ndraws + nsimul;

end	% loop for blocks

meaninvlike = suminvlike/ndraws;
gaminvlike  = gaminvlike/ndraws - repmat((meaninvlike.^2)',hmax,1);

% standard error
suminvlikeerror = gaminvlike(1,:);
for k=2:hmax
   suminvlikeerror = suminvlikeerror + 2*gaminvlike(k,:)*(1-(k-1)/hmax);
end

suminvlikeerror = 100*sqrt(suminvlikeerror/ndraws)./meaninvlike' ;

% marginalized data density
mdd = densfac-log(meaninvlike);
%mdd = [p mdd];

mdd = mean(mdd);

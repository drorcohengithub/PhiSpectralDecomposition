function [S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(G,split_mask_A,split_mask_E,K,fres,iter_max,gamma,min_error)

Cov_X = G(:,:,1);
Cov_XY = G(:,:,2:end);
N = size(Cov_X,1);
%%% obtain adjmat and sigma by minimizing GV
% for the split model we use the provided mask. 
% Intialize with
A_initial=zeros(N,N,K);

% Default parameters sets if not provided
if nargin<=4 
    % quit after this many interations
    iter_max = 6000;
end
if nargin<=5
    % step size when incrementaly calculating the reduced model
    gamma    = 0.01;
end

if nargin<=6
    % min_error. When this is reached we stop
    min_error    = 10^-6;
end

% run the optimization
% try mex
if exist('adjmat_argminGV_mask_mex','file')
    [SIG_r,A_r,masked_Delta] = adjmat_argminGV_mask_mex(Cov_X,Cov_XY,K,split_mask_A,split_mask_E,A_initial,iter_max,gamma,min_error);
else
    % look for matfile
    disp('mex file not available. To generate mex file see "Generate_mex_code.m". For now running slower .m script instead')
    [SIG_r,A_r,masked_Delta] = adjmat_argminGV_mask(Cov_X,Cov_XY,K,split_mask_A,split_mask_E,A_initial,iter_max,gamma,min_error);
end

%%% phi_VAR estimation in frequency domain
[S_r,~] = var_to_cpsd(A_r,SIG_r,fres);

% logdet(S)
det_S_r = zeros(size(S_r,3),1);
trace_S_r = det_S_r;
prod_diag_S_r = det_S_r;

% Potetially useful quantities
for freqIndx = 1:size(S_r,3)
    thisFreqS = squeeze(S_r(:,:,freqIndx));
    det_S_r(freqIndx) = det( thisFreqS );   
    trace_S_r(freqIndx) = trace( thisFreqS );     
    prod_diag_S_r(freqIndx) = prod(diag(thisFreqS));

end

 

end


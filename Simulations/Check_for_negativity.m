%%%%% Checking what causes the spectral decomposition to be negative
clear all;close all
%% create a bidirectionally connected system
N = 2; % # components
K = 2; % order

%% Residual covariance
SIG = [1 0; 0 0.7];

% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.2 -0.05;
            0.4 0.2]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.25 0;
            -0.2 0.1]; % connectivity matrix at time t-2

%% Get the autocov 

% Get all the characteristics of the system
[X,info] = var_to_autocov(A,SIG);

% Check that everything is ok
var_info(info,true);

% Obtain the frequency bins (arbitrary for simulated system)
% frequcy resolution
freq_res = 100;
% sampling rate 
samp_rate = 1000;
freqs = sfreqs(freq_res,samp_rate);

% We will need these too
% The cov of X
Cov_X  = X(:,:,1);

% The auto cov without cov X
Cov_XY = X(:,:,2:end);

% Set max lag or reduced model. This is sufficiently generous 
max_order=size(Cov_XY,3);

% The spectral density matrix of the full model
[S] = autocov_to_cpsd(X,freq_res);

%% Check that the model satisfies the decomposition requirements (Eq. 16 and, unbumbered one before it)
disp('Check that the FULL model satisfies the decomposition requirements')
H = var2trfun(A,freq_res);
H_det=[];

% calc det of product at each freq
for freq=1:size(H,3)
    tmp = squeeze(H(:,:,freq))*squeeze(H(:,:,freq))';
    H_det(freq) = det(tmp); 
end

% because the fft and everything so far has been one sided, and since the
% equation requires the two sided sum (-pi to pi)
two_sided_product = [H_det H_det(2:end-1)];
%this will be very small
disp('integral of log(det(H*conj(H))) (shud be zero)')
mean(log(two_sided_product))

%% for the spectral density matrix

S_det=[];

% calc det of product at each freq
for freq=1:size(H,3)
    tmp = squeeze(S(:,:,freq));
    S_det(freq) = det(tmp); 
end

% because the fft and everything so far has been one sided, and since the
% equation requires the two sided sum (-pi to pi)
two_sided_product = [S_det S_det(2:end-1)];
%this will be very small
fprintf('---------------------------\n')
disp('integral of log(det(S))')
mean(log(two_sided_product))
disp('log det of residual covariance matrix, shud be same as above')
log(det(SIG))
disp('the difference between the two')
mean(log(two_sided_product))-log(det(SIG))
fprintf('---------------------------\n')
%% Now move on to reduced model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Here we set the restriction on the autoregressive coefficient matrix
% we do this by setting zeros in split_mask_A for components of the 
% coefficient matrix we want to set to zero to values of the auto

%% These have positive spectral decompositions
% phi
split_mask_A = eye(N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% GC
split_mask_A=ones(N);
split_mask_A(2,1) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% GC other dir
split_mask_A=ones(N);
split_mask_A(1,2) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% These have negative spectral decompositions
% Pred info
split_mask_A=zeros(N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% %complement of phi
split_mask_A = ones(N);
split_mask_A(1,1) = 0;
split_mask_A(2,2) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% %settinet a self connection to zero
split_mask_A = ones(N);
split_mask_A(1,1) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% %or the other self connection
split_mask_A = ones(N);
split_mask_A(2,2) = 0;

%%
disp('the split mask is:')
split_mask_A

%%
% The code also allows for constratins on the value of the covariance of the residuals, using the 
% split_mask_SIG parameter. Though this works in practice we have not
% proved analytically that this is a legit thing to do. For this reason
% best to just make irrelvant for now.
split_mask_SIG=ones(N);

% some parameters for the optimzation
min_error = 1e-14; %when to stop. Smaller is more accurate
gamma = 0.1; %intial step size for optimization. This heuroestically changes 
            %if the optimization stalls
iter_max = 12000; %max num of iterations

% Both time and frequency domain quantities are estimated
fprintf('---------------------------\n')

disp('Computing reduced model parameters')
            
[S_r,... % spectral density matrix of reduced model
det_S_r,... %determinant of the spectral density matrix of reduced model
trace_S_r,... % trace of the spectral density matrix of reduced model
prod_diag_S_r,... % product of diag entriesof the spectral density matrix of reduced model
A_r,... % autoregressive coeff matrix of reduced model
SIG_r,... % covariace of the residuals of the reduced model
masked_Delta]... % errors during optimization covariace of the residuals of the reduced model
= get_reduced_S_from_autoCov(X,...
    split_mask_A,split_mask_SIG,max_order,freq_res,iter_max,gamma,min_error);    
disp('Optimization complete')
fprintf('---------------------------\n') 
%% Check if there is any issue with reduced model (see if there is warninigs from mvgc)
% MVGC will complain here for A = 0, because of some lag condition being
% violated. It is not a reason to panic.

fprintf('\n') 
disp('Any warnings from mvgc?')
[G info]=var_to_autocov(A_r,SIG_r);
var_info(info)

fprintf('---------------------------\n')
%%
disp('Check that the DISCONNECTED model satisfies the decomposition requirements')

H_r = var2trfun(A_r,freq_res);
H_r_det=[];

% calc det of product at each freq
for freq=1:size(H_r,3)
    tmp = squeeze(H_r(:,:,freq))*squeeze(H_r(:,:,freq))';
    H_r_det(freq) = det(tmp); 
end

% because the fft and everything so far has been one sided, and since the
% equation requires the two sided sum (-pi to pi)
two_sided_product = [H_r_det H_r_det(2:end-1)];
%this will be very small
disp('integral of log(det(H*conj(H))), shud be zero')
mean(log(two_sided_product))

%% for the spectral density matrix

S_r_det=[];

% calc det of product at each freq
for freq=1:size(H_r,3)
    tmp = squeeze(S_r(:,:,freq));
    S_r_det(freq) = det(tmp); 
end

% because the fft and everything so far has been one sided, and since the
% equation requires the two sided sum (-pi to pi)
two_sided_product = [S_r_det S_r_det(2:end-1)];
%this will be very small
fprintf('---------------------------\n')
disp('integral of log(det(S))')
mean(log(two_sided_product))
disp('log det of residual covariance matrix, shud be same as above')
log(det(SIG_r))
disp('the difference between the two')
mean(log(two_sided_product))-log(det(SIG_r))
fprintf('---------------------------\n')

%%
disp('the spectral decomposition')

[sdecomp_ratio time_domain_ratio det_S] = ratio_of_dets(S, S_r, SIG, SIG_r);

sdecomp_Phi_G = 1/2 *real(sdecomp_ratio);% note factor of half. may have tiny imaginery vals
Phi_G = 1/2 *time_domain_ratio;

clf
% may not be quite real
subplot(1,2,1)
bar(Phi_G)
ylim([0 2*Phi_G])
title('Phi_G time domain, this always seems fine') 

subplot(1,2,2)
plot(freqs,sdecomp_Phi_G)
title('Spectral decomp, this is sometime negative')
xlabel('Frequency')
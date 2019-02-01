clear all
close all

%% Some very simple autoreg system
%% system
N = 2; % # components
K = 2; % order

%% Uni directionally connected system with no inst connections
SIG = [1 0; 0 0.7];


% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.2 0;
            0.4 0.2]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.25 0;
            -0.2 0.1]; % connectivity matrix at time t-2

%% we need the autocov seq. For that we use MVGC

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

% We will need these 

% The cov of X
Cov_X  = X(:,:,1);

% The auto cov without cov X
Cov_XY = X(:,:,2:end);


% We need to choose the maximum lag we will calculate the reduced model up to.
% For good estimate the maximum lag of the reduced model will be
% much larger than that of the full model. A safe but potentially
% over-generous strategy is to use the max lag of the autocov function
max_order=size(Cov_XY,3);

% The spectral density matrix of the full model
[S] = autocov_to_cpsd(X,freq_res);

% store the results here. 
results = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Here we set the restriction on the autoregressive coefficient matrix
% we do this by setting zeros in split_mask_A for components of the 
% coefficient matrix we want to set to zero to values of the auto
% for example, for the 2D model we have here, we get Phi be using: 
% split_mask_A=[1 0;
%               0 1];
% More generally, the atomic partition is set using          
split_mask_A=eye(N);
  
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

%% log ratio
[sdecomp_Phi_G Phi_G det_S] = ratio_of_dets(S, S_r, SIG, SIG_r);

sdecomp_Phi_G = real(sdecomp_Phi_G);%may have tiny imaginery vals

clf
% may not be quite real
subplot(1,2,1)
bar(Phi_G)
title('Phi_G time domain') 

subplot(1,2,2)
plot(freqs,sdecomp_Phi_G)
title('Spectral decomp')
xlable('Frequency')
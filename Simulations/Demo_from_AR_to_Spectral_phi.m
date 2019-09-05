clear all
close all

%% A simple autoreg system

N = 2; % # components
K = 2; % AR order

%% Uni directionally connected system with no inst connections
% noise covariance matrix
Sigma_eps = [1 0; 0 0.7];


% AR coeff
A = zeros(N,N,K);
A(:,:,1) = [0.2 0;
            0.4 0.2]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.25 0;
            -0.2 0.1]; % connectivity matrix at time t-2

%% Step 1 - The spectral density of the full model

% some arbitraty settings for the simulation
% frequcy resolution
freq_res = 100;
% sampling rate 
samp_rate = 1000;

freqs = sfreqs(freq_res,samp_rate);

%compute the spectral density matrix
S = var_to_cpsd(A,Sigma_eps,freq_res);


%% Step 2 - Get the autocovariance of the full model

% Autocovariance 
[Sigma_X_x,info] = var_to_autocov(A,Sigma_eps);

% Check that everything is ok
var_info(info,true);


%% Step 3 - Get parameters of the disconnected model

% We need to choose the lag of the disconnected model.
% For good estimates, this lag may be much larger than that of the full
% model. A safe strategy is to use the max lag of the autocov function
disc_model_order=size(Sigma_X_x,3)-1;

% Here we specify the constraints of the system. We do this using the 
% split_mask_A parameter. We 'cut' causal connections by setting the 
% respective values to zero. For example, for the 2D model we have here, we 
% get cut both causal connections by setting 
split_mask_A=[1 0;
              0 1];
  
% The code also allows for constratins on the value of the covariance of
% % the residuals, using the by using the split_mask_SIG parameter. Though 
% this works in practice, we have not proved analytically that this is a
% sound thing to do. To ingore this parameters we impose no restrictions:
split_mask_Sigma_eps=ones(N);

% some parameters for the optimzation
min_error = 1e-14; %when to stop. Smaller is more accurate
gamma = 0.1; %intial step size for optimization. This heuroestically changes 
            %if the optimization stalls
iter_max = 12000; %max num of iterations

% Both time and frequency domain quantities are estimated
disp('Computing disconnected model parameters')
    

[ignore1,... % we'll get back to these later
ignore2,... % we'll get back to these later
ignore3,... % we'll get back to these later
ignore4,... % we'll get back to these later
A_prime,... % autoregressive coeff matrix of disconnected model
Sigma_eps_prime,... % covariace of the residuals of the disconnected model
delta]... % the change in the coeff of A_prime during optimization. If optimization works this should be small in the end
= get_reduced_S_from_autoCov(Sigma_X_x,split_mask_A,split_mask_Sigma_eps,disc_model_order,freq_res,iter_max,gamma,min_error);    

%% Step 4 - Get the spectral density matrix of the full model

S_prime = var_to_cpsd(A_prime,Sigma_eps_prime,freq_res);

%% log ratio of the quantities
[ci_f ci_t det_S] = ratio_of_dets(S, S_prime, Sigma_eps, Sigma_eps_prime);

% factor of half
Phi_G_spectral_decomp = 1/2 *real(ci_f);% note factor of half, may have tiny imaginery vals
Phi_G_time = 1/2 *ci_t;

%% plot
clf
subplot(1,2,1)
bar(Phi_G_time)
title('Phi_G time domain') 

subplot(1,2,2)
plot(freqs,Phi_G_spectral_decomp)
title('Phi_G Spectral decomposition')
xlabel('Frequency')

%% Note
if 0
    % The function get_reduced_S_from_autoCov already computes the spectral
    % density matrix of the disconnected model, so there is no need to do this
    % manually afterwards. Specifically:
    % 
    [S_prime,... % spectral density matrix of reduced model
    det_S_prime,...%determinant of the spectral density matrix of reduced model
    trace_S_prime,... % trace of the spectral density matrix of reduced model
    prod_diag_S_prime,... % product of diag entries of the spectral density matrix of reduced model
    A_prime,... % autoregressive coeff matrix of reduced model
    Sigma_eps_prime,... % covariace of the residuals of the reduced model
    masked_Delta]... % errors during optimization covariace of the residuals of the reduced model
    = get_reduced_S_from_autoCov(Sigma_X_x,split_mask_A,split_mask_Sigma_eps,disc_model_order,freq_res,iter_max,gamma,min_error);    

    %% recompute the quantities and plot on same graph as abovelog ratio of the quantities
    [ci_f ci_t det_S] = ratio_of_dets(S, S_prime, Sigma_eps, Sigma_eps_prime);

    % factor of half
    Phi_G_spectral_decomp = 1/2 *real(ci_f);% note factor of half, may have tiny imaginery vals
    Phi_G_time = 1/2 *ci_t;


    subplot(1,2,1)
    hold on
    bar(1,Phi_G_time,0.3,'r')
    title('Phi_G time domain') 

    subplot(1,2,2)
    hold on
    plot(freqs,Phi_G_spectral_decomp,'r--')
    title('Phi_G Spectral decomposition')
    xlabel('Frequency')
end
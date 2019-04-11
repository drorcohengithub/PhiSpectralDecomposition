clear all
close all

%% load systems
systems = Example_systems();

%% we can try the bivariate sys in simulations (1-5)

system_indx = 5;
A = systems(system_indx).A;
SIG = systems(system_indx).SIG;
split_mask_A = ones(size(SIG));
split_mask_A(2,1)=0; %zero connection from 1 to 2
%from y to x
y = 1;
x = 2;

%% or we can try a multivariate system across an arbitrary partition (sys #6)
A = systems(6).A;
 SIG = systems(6).SIG;

%this defines the constraint on A. Here for GC between elements 1,2 (termed y) and 3-5 (termed x) 
% y = 1:2;
% x = 3:5;
% split_mask_A = ones(size(SIG));
% split_mask_A(x,y)=0;



%% Get all the characteristics of the system
[X,info] = var_to_autocov(A,SIG); 

% Check that everything is ok
var_info(info,true);

% Obtain the frequency bins
freq_res = 100;
samp_rate = 1000;
freqs = sfreqs(freq_res,samp_rate);

% We will need these 
% The cov of X
Cov_X  = X(:,:,1);
N = size(Cov_X,1);

% The auto cov without cov X
Cov_XY = X(:,:,2:end);

% We need to choose the maximum lag we will calculate the reduced model up to.
% Note that for good estimate the maximum lag of the reduced model will be
% much larger than that of the full model. A safe but potentially
% over-generous strategy is to use the max lag of the autocov function
max_order=size(Cov_XY,3);

% The spectral density matrix of the full model
[S_f] = autocov_to_cpsd(X,freq_res);

% For GC we place no retrictions on the noise covariance
split_mask_E=ones(N);

%% Optimize the restricted model, 
% GC in first dir
iter_max = 12000;   
gamma = 0.1;
min_error = 1e-14;
    
[reslts_dir1.S_r,reslts_dir1.det_S_r,reslts_dir1.trace_S_r,reslts_dir1.prod_diag_S_r,...
    reslts_dir1.A_r,reslts_dir1.SIG_r,reslts_dir1.masked_Delta] = ...
    get_reduced_S_from_autoCov(X,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    

% other dir, just transpose the mask
 
split_mask_A = split_mask_A';

reslts_dir2 = [];

[reslts_dir2.S_r,reslts_dir2.det_S_r,reslts_dir2.trace_S_r,reslts_dir2.prod_diag_S_r,...
    reslts_dir2.A_r,reslts_dir2.SIG_r,reslts_dir2.masked_Delta] = ...
    get_reduced_S_from_autoCov(X,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    

%% get the log ratios of the generalized variances and spectrl adensity matrices
[dir1_ratio_S dir1_ratio] = ratio_of_dets(S_f, reslts_dir1.S_r, SIG, reslts_dir1.SIG_r);  
[dir2_ratio_S dir2_ratio] = ratio_of_dets(S_f, reslts_dir2.S_r, SIG, reslts_dir2.SIG_r);      

%% Now, compute the quantities using the MVGC toolbox
[dir1_t_GC_mvgc]= autocov_to_mvgc(X,x,y);
[dir2_t_GC_mvgc]= autocov_to_mvgc(X,y,x);

[dir1_s_GC_mvgc]= autocov_to_smvgc(X,x,y,freq_res);
[dir2_s_GC_mvgc]= autocov_to_smvgc(X,y,x,freq_res);

%% plot
clf
subplot(2,2,1)
bar([dir1_ratio dir1_t_GC_mvgc])
set(gca,'xticklabel',{'Our framework','MVGC'})
ylim([0 0.02])
title({'dir1';['diff between measures is ' num2str(dir1_ratio-dir1_t_GC_mvgc)]})

subplot(2,2,2)
plot(dir1_s_GC_mvgc)
hold on
plot(real(dir1_ratio_S),'r--')
% ylim([-0.0001 0.0001])
xlabel('Frequency')
xlim([0 length(dir1_ratio_S) ])
title('dir 1 spectral decomp ')

subplot(2,2,3)
title('spectral')
bar([dir2_ratio dir2_t_GC_mvgc])
set(gca,'xticklabel',{'Our framework','MVGC'})
title({'dir2';['diff between measures is ' num2str(dir2_ratio-dir2_t_GC_mvgc)]})
% ylim([0 0.0001])

subplot(2,2,4)
plot(dir2_s_GC_mvgc)
hold on
plot(real(dir2_ratio_S),'r--')
legend({'MVGC', 'Our framework'})
xlabel('Frequency')
title('dir 2 spectral decomp ')
xlim([0 length(dir2_s_GC_mvgc) ])

%%
% Below we check a bunch of the equalities in the paper
%%
%% Equivalences between the means and covs of the conditionals of the full and reduced model
part1=y;
part2=x;

disp('')
disp('Checking time domain equalities between the full and reduced models')
disp('Conditional mean of full model:')
meanCondFull = SIG(part1,part2)*SIG(part2,part2)^-1;
disp(meanCondFull)

disp('Conditional mean of disconnected model')
meadnCondPart1 = reslts_dir1.SIG_r(part1,part2)*reslts_dir1.SIG_r(part2,part2)^-1;
disp(meadnCondPart1)

disp('Their difference')
disp(meanCondFull - meadnCondPart1);

%%
disp('Conditional variance of full model:')
dir1_full_partcov = SIG(part1,part1) - SIG(part1,part2)*SIG(part2,part2)^-1*SIG(part2,part1);
disp(dir1_full_partcov)

disp('Conditional variance of disconnected model:')
dir1_partcov=reslts_dir1.SIG_r(part1,part1) - reslts_dir1.SIG_r(part1,part2)*reslts_dir1.SIG_r(part2,part2)^-1*reslts_dir1.SIG_r(part2,part1);
disp(dir1_partcov)

disp('Their diff:')
disp(dir1_partcov - dir1_full_partcov)

%% other direction
meanCondFull = SIG(part2,part1)*SIG(part1,part1)^-1;
meadnCondPart2 = reslts_dir2.SIG_r(part2,part1)*reslts_dir2.SIG_r(part1,part1)^-1;
disp('')
disp('Difference betwee cond mean in full model and disconnected model for other direction GC')
disp(meanCondFull - meadnCondPart2);

dir2_partcov=reslts_dir2.SIG_r(part2,part2) - reslts_dir2.SIG_r(part2,part1)*reslts_dir2.SIG_r(part1,part1)^-1*reslts_dir2.SIG_r(part1,part2);
dir2_full_partcov = SIG(part2,part2) - SIG(part2,part1)*SIG(part1,part1)^-1*SIG(part1,part2);
disp('Difference betwee cond Variance in full model and disconnected model for other direction GC')
disp(dir2_partcov - dir2_full_partcov)

%% Check the relatinoshup between the autoregressive parameters of the full and disconnected model
% First we will need to create the quantities as per the paper
% This is the covariance matrix of the system has a Big matrix 
K = size(Cov_XY,3);
Cov_B = zeros(N*K,N*K);
for i=1: K
    for j=i: K
        if i== j
            Cov_B((i-1)*N+1:i*N,(j-1)*N+1:j*N) = Cov_X;
        else
            Cov_B((i-1)*N+1:i*N,(j-1)*N+1:j*N) = Cov_XY(:,:,(j-i));
            Cov_B((j-1)*N+1:j*N,(i-1)*N+1:i*N) = Cov_XY(:,:,(j-i))';
        end

    end
end

Cov_XY_B = zeros(N,N*K);
for j=1: K
    Cov_XY_B(:,(j-1)*N+1:j*N) = Cov_XY(:,:,j);
end

%% Now get the corresponding big matrices for the ar coeffs
% to make sure the lags match the disconnecte system
A_tmp = reslts_dir1.A_r * 0;
A_tmp(:,:,1:2) = A;

% get big A for full system
A_p=A_tmp;
A_f_B = zeros(N,N*K);
for j=1: K
    A_f_B(:,(j-1)*N+1:j*N) = A_p(:,:,j);
end

% get big A for disc system in dir 1
A_p=reslts_dir1.A_r;
A_dir1_B = zeros(N,N*K);
for j=1: K
    A_dir1_B(:,(j-1)*N+1:j*N) = A_p(:,:,j);
end

% get big A for disc system in dir 2
A_p=reslts_dir2.A_r;
A_dir2_B = zeros(N,N*K);
for j=1: K
    A_dir2_B(:,(j-1)*N+1:j*N) = A_p(:,:,j);
end

%% Now, difference between the var in terms of delta A
delta_A_dir_1 = A_f_B - A_dir1_B;
sig_diff_1 = delta_A_dir_1 * Cov_B * delta_A_dir_1';
disp('Diff between the variance in full and disconnected model version 1')
disp(sig_diff_1)

disp('Diff between the variance in full and disconnected model version 1')
sig_diff_2 = reslts_dir1.SIG_r - SIG;
disp(sig_diff_2)

disp('difference between the two is')
disp(sig_diff_1-sig_diff_2)

disp('othe dir')
delta_A_dir_2 = A_f_B - A_dir2_B;
sig_diff_1 = delta_A_dir_2 * Cov_B * delta_A_dir_2';
disp('Diff between the variance in full and disconnected model version 1')
disp(sig_diff_1)

disp('Diff between the variance in full and disconnected model version 1')
sig_diff_2 = reslts_dir2.SIG_r - SIG;
disp(sig_diff_2)

disp('difference between the two is')
disp(sig_diff_1-sig_diff_2)
%% Now together with transformation by P
P_dir_1 = eye(size(SIG));
P_dir_1(part1,part2) = SIG(part1,part2)*SIG(part2,part2)^-1;
P_dir_1_I = P_dir_1^-1;
sig_diff_1=P_dir_1_I*(reslts_dir1.SIG_r - SIG)*P_dir_1_I';
sig_diff_2 = P_dir_1_I*delta_A_dir_1 * Cov_B * delta_A_dir_1'*P_dir_1_I';
disp('Difference between two expressions')
disp(sig_diff_1-sig_diff_2)

% confirm these are zero
tmp=P_dir_1_I*delta_A_dir_1;
disp('Transformed A is zero:')
disp(tmp(part1,part1))


%% other dir
P_dir2 = eye(size(SIG));
P_dir2(part2,part1) = SIG(part2,part1)*SIG(part1,part1)^-1;
P_dir2_I = P_dir2^-1;

sig_diff_1 = P_dir2_I*delta_A_dir_2 * Cov_B * delta_A_dir_2'*P_dir2_I';
sig_diff_2 = P_dir2_I*(reslts_dir2.SIG_r - SIG)*P_dir2_I';

disp('Other dir')
disp('Difference between two expressions')
disp(sig_diff_1-sig_diff_2)

tmp=P_dir2_I*delta_A_dir_2;
disp('Transformed A is zero:')
disp(tmp(part2,part2))

%% now, frequency domain
H = var2trfun(A,freq_res);
H_dir1 = var2trfun(reslts_dir1.A_r,freq_res);
H_dir2 = var2trfun(reslts_dir2.A_r,freq_res);

chk1 = 0 * H;
chk2 = 0 * H;

for freq_indx = 1:size(H,3)
          chk1(:,:,freq_indx) =  squeeze(H(:,:,freq_indx)*P_dir_1)^(-1) - squeeze(H_dir1(:,:,freq_indx)*P_dir_1)^(-1);         
          chk2(:,:,freq_indx) =  squeeze(H(:,:,freq_indx)*P_dir2)^(-1) - squeeze(H_dir2(:,:,freq_indx)*P_dir2)^(-1);         
end

tmp = chk1(part1,part1,:);
disp('Overall magnitute of of dir 1')
disp(sum(abs(tmp(:))));

tmp = chk2(part2,part2,:);
disp('Overall magnitute of of dir 2')
disp(sum(abs(tmp(:))));

%% Finally

H_tilde = 0 * H;
Hs11_tilde_converted = 0 * H;
Hs22_tilde_converted = 0 * H;

for freq = 1:size(H,3)
    
     H_tilde(:,:,freq) = squeeze(H(:,:,freq))*P_dir_1;
     Hs11_tilde_converted(part1,part1,freq) = H_tilde(part1,part1,freq) -...
          H_tilde(part1,part2,freq) * squeeze(H_tilde(part2,part2,freq))^-1 * H_tilde(part2,part1,freq); 
   
      H_tilde(:,:,freq) = squeeze(H(:,:,freq))*P_dir2;
      Hs22_tilde_converted(part2,part2,freq) = H_tilde(part2,part2,freq) -...
           H_tilde(part2,part1,freq) * squeeze(H_tilde(part1,part1,freq))^-1 * H_tilde(part1,part2,freq); 
   
end

%%
disp('Checking transfer function equality')

diff_chk = Hs11_tilde_converted(part1,part1,:) - H_dir1(part1,part1,:); 
disp(sum(abs(diff_chk(:))));

diff_chk = Hs22_tilde_converted(part2,part2,:) - H_dir2(part2,part2,:); 
disp(sum(abs(diff_chk(:))));

%If these are not multivariate we can plot them
% clf
% subplot(1,2,1)
% plot(squeeze(Hs11_tilde_converted(1,1,:)))
% plot(squeeze(H_dir1(1,1,:)))
% 
% subplot(1,2,2)
% plot(squeeze(Hs22_tilde_converted(2,2,:)))
% plot(squeeze(H_dir2(2,2,:)))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Predictive information 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% To obtain predictive information all lagged ineractions set to zero. 

% no lagged influences
split_mask_A=zeros(N);

% no restriction on the noise covariance 
split_mask_E=ones(N);

% calc reduced model
[S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(X,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    


%% log ratios
[ratio_S ratio det_S_f] = ratio_of_dets(S_f, S_r, SIG, SIG_r);
   
freq_tm_diff = abs(ratio - mean([ratio_S;ratio_S(2:end-1)]));

% keep results here 
t_pred_info = ratio;
t_pred_info_spct = ratio_S;

% compare against analytic result ( MI past present )
t_pred_info2=log( det(Cov_X)/det(SIG) );
t_pred_info_spct2=log( det(Cov_X)./real(det_S_f) );
    
clf
subplot(1,2,1)
bar([t_pred_info t_pred_info2])
title({['tm domain. Diff between analltic and numerical opt is: ' num2str(t_pred_info-t_pred_info2)];...
    ['diff between tm and mean over freqs:' num2str(freq_tm_diff)]})

subplot(1,2,2)
ph = plot(freqs,real(t_pred_info_spct),freqs, t_pred_info_spct2,'--r'); 
set(ph, 'linewidth', 3)
legend({'Numerical optimization','analytic based'})
ylabel('pred. info.')
xlabel('Freq')      
title('spct decomp')
   



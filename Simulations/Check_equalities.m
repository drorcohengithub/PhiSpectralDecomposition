clear all
close all

%% load systems
systems = Example_systems();

A = systems(6).A;
SIG = systems(6).SIG;
split_mask_A = ones(size(SIG));
split_mask_A(3:5,1:2)=0;
% 
% %For MVGC from y to x
y = 1:2;
x = 3:5;

% A = systems(5).A;
% SIG = systems(5).SIG;
% split_mask_A = ones(size(SIG));
% split_mask_A(2,1)=0; %zero connection from 1 to 2
% %from y to x
% y = 1;
% x = 2;

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

%% Optimize
iter_max = 12000;   
gamma = 0.1;
min_error = 1e-14;
    
[reslts_dir1.S_r,reslts_dir1.det_S_r,reslts_dir1.trace_S_r,reslts_dir1.prod_diag_S_r,...
    reslts_dir1.A_r,reslts_dir1.SIG_r,reslts_dir1.masked_Delta] = ...
    get_reduced_S_from_autoCov(X,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    

%other dir

split_mask_A = split_mask_A';

reslts_dir2 = [];

[reslts_dir2.S_r,reslts_dir2.det_S_r,reslts_dir2.trace_S_r,reslts_dir2.prod_diag_S_r,...
    reslts_dir2.A_r,reslts_dir2.SIG_r,reslts_dir2.masked_Delta] = ...
    get_reduced_S_from_autoCov(X,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    
%% compe the ratios
[dir1_ratio_S dir1_ratio] = ratio_of_dets(S_f, reslts_dir1.S_r, SIG, reslts_dir1.SIG_r);  
[dir2_ratio_S dir2_ratio] = ratio_of_dets(S_f, reslts_dir2.S_r, SIG, reslts_dir2.SIG_r);      
%% Compare with MVGC - (factor of half)

[dir1_t_GC_mvgc]= autocov_to_mvgc(X,x,y);
[dir2_t_GC_mvgc]= autocov_to_mvgc(X,y,x);

[dir1_s_GC_mvgc]= autocov_to_smvgc(X,x,y,freq_res);
[dir2_s_GC_mvgc]= autocov_to_smvgc(X,y,x,freq_res);

%% plot
clf
subplot(2,2,1)
bar([dir1_ratio dir1_t_GC_mvgc])
set(gca,'xticklabel',{'Framework','MVGC'})
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
set(gca,'xticklabel',{'Framework','MVGC'})
title({'dir2';['diff between measures is ' num2str(dir2_ratio-dir2_t_GC_mvgc)]})
% ylim([0 0.0001])

subplot(2,2,4)
plot(dir2_s_GC_mvgc)
hold on
plot(real(dir2_ratio_S),'r--')
legend({'MVGC', 'Framework'})
xlabel('Frequency')
title('dir 2 spectral decomp ')
xlim([0 length(dir2_s_GC_mvgc) ])


%% Equivalences between the means and covs of the conditionals of the full and reduced model
part1=1:2;
part2=3:5;

disp('Checking time domain equality')

meanCondFull = SIG(part1,part2)*SIG(part2,part2)^-1;
meadnCondPart1 = reslts_dir1.SIG_r(part1,part2)*reslts_dir1.SIG_r(part2,part2)^-1;
disp(meanCondFull - meadnCondPart1);

dir1_partcov=reslts_dir1.SIG_r(part1,part1) - reslts_dir1.SIG_r(part1,part2)*reslts_dir1.SIG_r(part2,part2)^-1*reslts_dir1.SIG_r(part2,part1);
dir1_full_partcov = SIG(part1,part1) - SIG(part1,part2)*SIG(part2,part2)^-1*SIG(part2,part1);
disp(dir1_partcov - dir1_full_partcov)

meanCondFull = SIG(part2,part1)*SIG(part1,part1)^-1;
meadnCondPart2 = reslts_dir2.SIG_r(part2,part1)*reslts_dir2.SIG_r(part1,part1)^-1;
disp(meanCondFull - meadnCondPart2);

dir2_partcov=reslts_dir2.SIG_r(part2,part2) - reslts_dir2.SIG_r(part2,part1)*reslts_dir2.SIG_r(part1,part1)^-1*reslts_dir2.SIG_r(part1,part2);
dir2_full_partcov = SIG(part2,part2) - SIG(part2,part1)*SIG(part1,part1)^-1*SIG(part1,part2);
disp(dir2_partcov - dir2_full_partcov)

%% numerator 
H=var2trfun(A,freq_res);
dir1_Hprime=var2trfun(reslts_dir1.A_r,freq_res);
dir2_Hprime=var2trfun(reslts_dir2.A_r,freq_res);

%% Check freq domain equivalence


P = eye(size(SIG));
P(part1,part2) = SIG(part1,part2)*SIG(part2,part2)^-1;


H_tilde = 0 * H;
Hs11_tilde_converted = 0 * H;
for freq = 1:size(H,3)
    
     H_tilde(:,:,freq) = squeeze(H(:,:,freq))*P;
     Hs11_tilde_converted(part1,part1,freq) = H_tilde(part1,part1,freq) -...
          H_tilde(part1,part2,freq) * squeeze(H_tilde(part2,part2,freq))^-1 * H_tilde(part2,part1,freq); 
   
    
end

%%
disp('Checking transfer function relation 1')

clf
diff_chk = Hs11_tilde_converted(part1,part1,:) - dir1_Hprime(part1,part1,:); 
plot(abs(diff_chk(:)));
ylim([-0.0001 0.0001])

%%
P = eye(size(SIG));
P(part2,part1) = SIG(part2,part1)*SIG(part1,part1)^-1;


H_tilde = 0 * H;
Hs22_tilde_converted = 0 * H;
for freq = 1:size(H,3)
    
     H_tilde(:,:,freq) = squeeze(H(:,:,freq))*P;
     Hs22_tilde_converted(part2,part2,freq) = H_tilde(part2,part2,freq) -...
          H_tilde(part2,part1,freq) * squeeze(H_tilde(part1,part1,freq))^-1 * H_tilde(part1,part2,freq); 
%     
%  Hs22_tilde_converted(part2,part2,freq) = H(part2,part2,freq) -...
%          H_tilde(part2,part1,freq) * squeeze(H_tilde(part1,part1,freq))^-1 * H(part1,part2,freq); 
%     
    
end

%%
disp('Checking transfer function relation 2')
clf
diff_chk = Hs22_tilde_converted(part2,part2,:) - dir2_Hprime(part2,part2,:); 
plot(abs(diff_chk(:)));
ylim([-0.0001 0.0001])

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
    
% keep results here 
t_pred_info = ratio;
t_pred_info_spct = ratio_S;

% compare against analytic result ( MI past present )
t_pred_info2=log( det(Cov_X)/det(SIG) );
t_pred_info_spct2=log( det(Cov_X)./real(det_S_f) );
    
clf
subplot(1,2,1)
bar([t_pred_info t_pred_info2])
title(['tm domain. Diff between analltic and numerical opt is: ' num2str(t_pred_info-t_pred_info2)])

subplot(1,2,2)
ph = plot(freqs,real(t_pred_info_spct),freqs, t_pred_info_spct2,'--r'); 
set(ph, 'linewidth', 3)
legend({'Numerical optimization','analytic based'})
ylabel('pred. info.')
xlabel('Freq')      
title('spct decomp')
   



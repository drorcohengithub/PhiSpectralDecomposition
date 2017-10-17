clear all
close all

%% load systems
systems = Example_systems();
%% choose system you want
systemIndx = 5;
A = systems(systemIndx).A;
SIG_f = systems(systemIndx).SIG_f; 

%% Get the autocov of the full system
[G,info] = var_to_autocov(A,SIG_f);
% Check that everything is ok
var_info(info,true);

%% input for analysis
% arbitrarily set the freq resolution and samp rate
freq_res  = 100;
samp_rate = 1000;
% obtan the frequency bins
freqs = sfreqs(freq_res,samp_rate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the time domain GC 
[t_GC]= autocov_to_pwcgc(G);
GC_x1_to_x2 = t_GC(2,1);
GC_x2_to_x1 = t_GC(1,2);

%% Phi G calculation
%Covriance and autcov function for phi calc
Cov_X  = G(:,:,1);
Cov_XY = G(:,:,2:end);
% The maximum lag we will calculate the reduced model up to.
% Note that for good estimate the maximum lag of the reduced model will be
% much larger than that of the full model. A safe but potentially
% over-generous strategy is to use the max lag of the autocov function
max_order=size(Cov_XY,3);

%% Chose form of reduced model. 
% The split_mask determines which connections to cut and which to keep. By
% cutting both lagged connections we will obtain Phi_G
N = size(G,1);
split_mask= eye(N);

%% Get the parameters of the reduced model (denoted with subscript r)
% Both time and frequency domain quantities are estimated
[S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(G,split_mask,max_order,freq_res);

%% Finally, from the defintion of Phi G
Phi_G = log (det(SIG_r) / det (SIG_f));

%% Multi information
% The current code does not support restrictions on the noise covariance
% matrix (will add in next version). However we can obtain it in a a number
% of equivalent ways. The first options is to fit the models seperately
univar_autocov1=G(1,1,:);
[univar_A1 univar_SIG1] = autocov_to_var(univar_autocov1);
univar_autocov2=G(2,2,:);
[univar_A2 univar_SIG2] = autocov_to_var(univar_autocov2);
multi_info = log ( univar_SIG1*univar_SIG2  / det (SIG_f) );

% The second option is to get this from the existing GC estimates
log_E22=GC_x1_to_x2+log(SIG_f(2,2));
log_E11=GC_x2_to_x1+log(SIG_f(1,1));
multi_info_2 = log_E22+log_E11 - log( det( SIG_f ) );

% the difference between the two is very small
disp(multi_info-multi_info_2)

%% plot time domain
time_data = [multi_info Phi_G GC_x1_to_x2 GC_x2_to_x1 GC_x1_to_x2+GC_x2_to_x1];   
if systemIndx<3 % for systems in which we are not interested in the sum of GC
    max_iter = length(time_data)-1;
    offset = 0;
else % for systems in which we'd like to see the sum of GC
    max_iter = length(time_data);
    offset = 0.005;
end

clf
subplot(2,1,1);
set(gca,'fontsize',20)
hold on   
for i = 1:max_iter
    h=bar(i,time_data(i));
    text(i,time_data(i)+offset,num2str(time_data(i)),'fontsize',30)    
    if i == 1
        set(h,'FaceColor','g');
    elseif i == 2
        set(h,'FaceColor','k');
    elseif i == 3 || i == 4
        set(h,'FaceColor','r');
    elseif i == 5
        set(h,'FaceColor','b');
    end
end

set(gca,'xtick',1:max_iter)
set(gca,'xticklabel',{'multi info','PhiG','GC1to2','GC2to1','sum of GCs'})
maxy= max(ylim);
ylim([0 maxy])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spectral decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get the spectral decomposition of GC 
[spct_GC]= autocov_to_spwcgc(G,freq_res);
sdecomp_GC_x1_to_x2 = squeeze(spct_GC(2,1,:));
sdecomp_GC_x2_to_x1 = squeeze(spct_GC(2,1,:));

%% get the spectral decomposition multi info
% (directly from the spectral density matrix) (See Eq. in paper) 
%The spectral density matrix
[S] = autocov_to_cpsd(G,freq_res);

% Calculate the determinant
det_S=[];
for freqIndx=1:size(S,3)
    det_S(freqIndx)=det(squeeze(S(:,:,freqIndx)));
end

%% The defintion of the decompistion of multi information
sdecomp_multi_info = real(log(squeeze(S(1,1,:).*S(2,2,:))./det_S'));

%% Spectral decomposition of Phi G
sdecomp_Phi_G  = log ( real(det_S_r./det_S') );

%% plot frequency domain
subplot(2,1,2);

if systemIndx<3
    hold on
    hl=plot(freqs,squeeze(spct_GC(1,2,:)),'r',freqs,squeeze(spct_GC(2,1,:)),'r--');
    set(hl,'linewidth',4);
    % hl=plot(freqs,squeeze(spct_GC(1,2,:))+squeeze(spct_GC(2,1,:)),'b-.');
    % set(hl,'linewidth',2);
    plot(freqs,sdecomp_multi_info,'g');
    plot(freqs,sdecomp_Phi_G,'k:','linewidth',2)
    legend({'GC x_2 to x_1','GC x_1 to x_2','multi. info.','Phi G'},'fontsize',20,'linewidth',2)
    xlabel('Frequency')
    set(gca,'fontsize',20)
else
    hold on
    hl=plot(freqs,squeeze(spct_GC(1,2,:)),'r',freqs,squeeze(spct_GC(2,1,:)),'r--');
    set(hl,'linewidth',4);
    hl=plot(freqs,squeeze(spct_GC(1,2,:))+squeeze(spct_GC(2,1,:)),'b-.');
    set(hl,'linewidth',2);
    plot(freqs,sdecomp_multi_info,'g');
    plot(freqs,sdecomp_Phi_G,'k:','linewidth',2)
    legend({'GC x_2 to x_1','GC x_1 to x_2','sum of GCs','multi. info.','Phi G'},'fontsize',20,'linewidth',2)
    xlabel('Frequency')
    set(gca,'fontsize',20)
%     plot(freqs,squeeze(real(H(2,1,:))));

end


%% save
svFldr = ['/Users/dror/Google Drive/Sasai/FIGURES/Simulations/' systems(systemIndx).nm];
saveFigAs(svFldr,0,'halfScreen')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The code below verfies some analytical relationships presented in the paper. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check: Relationshup between spectral decomp of multi info and coherence
C=squeeze((abs(S(1,2,:)).^2)./(S(1,1,:).*S(2,2,:)));
%mutual info rate
sdecomp_multi_info_2 = -log(1-C);
clf
plot(freqs,sdecomp_multi_info-sdecomp_multi_info_2); ylim(1e-10*[-1 1]);%  note the scale

%% Check: Get GC through through reducde model framework
% do this by enforce off-diagonal to zero for GC
split_mask= ones(N);
split_mask(2,1)=0;
    
%this takes longer to converge
iter_max = 10000;
gamma = 0.01;
min_error = 1e-6;
[S_r_GC,det_S_r_GC,trace_S_r_GC,prod_diag_S_r_GC,A_r_GC,SIG_r_GC,masked_Delta_GC] =...
    get_reduced_S_from_autoCov(G,split_mask,max_order,freq_res,iter_max,gamma,min_error);
    
%compare with MVGC
GC_x1_to_x2_2=log( det (SIG_r_GC) / det (SIG_f) );
%time domain
disp(GC_x1_to_x2-GC_x1_to_x2_2);

%freq domain
sdecomp_GC_x1_to_x2_2=log(det_S_r_GC./det_S');
plot(freqs,real(sdecomp_GC_x1_to_x2_2)-sdecomp_GC_x1_to_x2) % this can be reduced by decreasing the min_error arg

%check summation across freqs
GC_x1_to_x2_3 = mean([sdecomp_GC_x1_to_x2_2;sdecomp_GC_x1_to_x2_2(2:freq_res)]);
disp(real(GC_x1_to_x2_3 - GC_x1_to_x2_2)) % this can be reduced by decreasing the min_error arg
%% other dirn also works, just replace:
%split_mask= ones(N);
%split_mask(2,1)=0;
% With
%split_mask= ones(N);
%split_mask(1,2)=0;
   
%% Check: Using power spectrums to estimate GC.
% Typically, to estimate GC we would have to fit a reduced model. However
% it turns out that we can simply use the power spectrum to estimate the 
% variance 

% variance of the individual models based on the power spct (two sided)
two_sided_S = cat(3,S,S(:,:,2:freq_res));
% We don't care about these here
two_sided_S(1,2,:) = 0;
two_sided_S(2,1,:) = 0;

% sum of log 
tmpS = squeeze(mean(log(two_sided_S),3));
log_E22_1 = tmpS(2,2); 
log_E11_1 = tmpS(1,1);
% compare with log_E22 and log_E11 as estimated above
disp(log_E22_1 - log_E22);
disp(log_E11 - log_E11);

%% Check: correspondance between time and freq domain versions of multi info and Phi G
% Info geom def of multi info
sbased_multi_info = mean([sdecomp_multi_info; sdecomp_multi_info(2:freq_res)]);
disp(sbased_multi_info - multi_info)

sbased_Phi_G = mean([sdecomp_Phi_G; sdecomp_Phi_G(2:freq_res)]);
disp(sbased_Phi_G - Phi_G)



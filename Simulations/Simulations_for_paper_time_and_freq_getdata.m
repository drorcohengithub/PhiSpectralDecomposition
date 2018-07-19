clear all
close all

%% load systems
systems = Example_systems();

%% choose system you want to investigate
system_indx = 5;
A = systems(system_indx).A;
SIG_f = systems(system_indx).SIG_f; 

% Get all the characteristics of the system
[G,info] = var_to_autocov(A,SIG_f);

% Check that everything is ok
var_info(info,true);

% Arbitrarily set the freq resolution and samp rate
freq_res  = 100;
samp_rate = 1000;

% Obtain the frequency bins
freqs = sfreqs(freq_res,samp_rate);

% We will need these 
% The cov of X
Cov_X  = G(:,:,1);
% The auto cov without cov X
Cov_XY = G(:,:,2:end);

% We need to choose the maximum lag we will calculate the reduced model up to.
% Note that for good estimate the maximum lag of the reduced model will be
% much larger than that of the full model. A safe but potentially
% over-generous strategy is to use the max lag of the autocov function
max_order=size(Cov_XY,3);

% The spectral density matrix of the full model
[S_f] = autocov_to_cpsd(G,freq_res);

% What to do
% Compute GC (both directions)
JOB.gc=1;

% Compute stochastic interaction
JOB.stochasticint=1;

% Compute predictive information
JOB.predinfo=1;

% Compute predictive information
JOB.instcaus=1;

% Compute phi
JOB.phig=1;

% Some parameters for the optimization. These seem to produced reasonable
% results for all metrics and simulations, but you can probably get better results 
% by adjusting on a case by case basis 
%max iter size
iter_max = 12000;
%learning rate, this is adjusted automatically if convergence stagnates
gamma = 0.1;
%min error, stop optimizing if reached
min_error = 1e-10;
%dimension of system
N = size(G,1);


%% Start the work
if JOB.gc
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Granger causality
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %The split_mask_A variable determines which of the autoregress coefficients are set to zero
    %The split_mask_E variable determines which of the noise covariance terms are set to zero 
    % To obtain GC, we cut a lagged connection from on var to the other
    split_mask_A=ones(N);
    
    % GC from first to second variable
    split_mask_A(2,1)=0;
    
    % For GC we place no retrictions on the noise covariance
    split_mask_E=ones(N);

    % calc reduced model
    if ~exist(['./' systems(system_indx).nm '_GC2to1.mat'],'file')
        % Both time and frequency domain quantities are estimated
        disp('Computing reduced model parameters')

        [S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(G,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    

        % save the quantities so that we don't have to run this again
        save(['./' systems(system_indx).nm '_GC2to1.mat'],'S_r','det_S_r','trace_S_r','prod_diag_S_r','A_r','SIG_r','masked_Delta',...
            'iter_max','gamma','gamma','split_mask_A','split_mask_E','max_order')    

    else

        disp('Loading reduced model parameters')
        load(['./' systems(system_indx).nm '_GC2to1.mat'])    
    end
    
    
    % Evaluate GC using main equation, time and freq
    [ratio_S ratio] = ratio_of_dets(S_f, S_r, SIG_f, SIG_r);
    
    % keep GC results here 
    t_GC = nan(N);
    
    % first direction gc
    t_GC(2,1) = ratio;

    % spectral decomp
    t_GC_spct = nan(N,N,length(ratio_S));
    t_GC_spct(2,1,:) = ratio_S;
    
    % for checks later
    pcUs1=SIG_r(1,1)-SIG_r(1,2)/SIG_r(2,2)*SIG_r(2,1);
    pcUs2=SIG_r(2,2)-SIG_r(1,2)/SIG_r(1,1)*SIG_r(2,1);
    

    %% now get GC in the other direction
    % to do this we cut connection in the other direction
    split_mask_A=ones(N);
    
    % note the flip below in order to get the opposite direction
    split_mask_A(1,2)=0;
    
    % again, not restriction on noise cov
    split_mask_E=ones(N);

    % calc reduced model
    if ~exist(['./' systems(system_indx).nm '_GC1to2.mat'],'file')
        disp('Computing reduced model parameters')

        [S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(G,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    
        
        save(['./' systems(system_indx).nm '_GC1to2.mat'],'S_r','det_S_r','trace_S_r','prod_diag_S_r','A_r','SIG_r','masked_Delta',...
            'iter_max','gamma','gamma','split_mask_A','split_mask_E','max_order')    

    else

        disp('Loading reduced model parameters')
        load(['./' systems(system_indx).nm '_GC1to2.mat'])    
    end

    % log ratios
    [ratio_S ratio] = ratio_of_dets(S_f, S_r, SIG_f, SIG_r);
    t_GC(1,2) = ratio;
    t_GC_spct(1,2,:) = ratio_S;

    %% get GC from mvgc to compare
    [t_GC_mvgc]= autocov_to_pwcgc(G);
    GC_x1_to_x2 = t_GC_mvgc(2,1);
    GC_x2_to_x1 = t_GC_mvgc(1,2);

    [spct_GC]= autocov_to_spwcgc(G,freq_res);
    sdecomp_GC_x1_to_x2 = squeeze(spct_GC(2,1,:));
    sdecomp_GC_x2_to_x1 = squeeze(spct_GC(1,2,:));
    
    % check partial covariances
    pc1=SIG_f(1,1)-SIG_f(1,2)/SIG_f(2,2)*SIG_f(2,1);
    pc2=SIG_f(2,2)-SIG_f(1,2)/SIG_f(1,1)*SIG_f(2,1);
    
    pcUs3=SIG_r(1,1)-SIG_r(1,2)/SIG_r(2,2)*SIG_r(2,1);
    pcUs4=SIG_r(2,2)-SIG_r(1,2)/SIG_r(1,1)*SIG_r(2,1);
    
    pc1
    pc2
    pcUs1
    pcUs2
    pcUs3
    pcUs4
    
    %% plot
    clf
    % first direction
    subplot(1,2,1)
    
    % freq domain plot
    ph1to2 = plot(freqs,sdecomp_GC_x1_to_x2,freqs,real(squeeze(t_GC_spct(2,1,:))),'--r');
    set(ph1to2, 'linewidth', 3)
    
    %time domain quantities
    title(sprintf('GC x1 to x2. Ours: %d, MVGC: %d',t_GC(2,1), GC_x1_to_x2));
    
    ylim([0 0.2])
    ylabel('spct gc 1 to 2')
    xlabel('Freq')

    %other direction
    subplot(1,2,2)
    ph2to1 = plot(freqs,sdecomp_GC_x2_to_x1,freqs,real(squeeze(t_GC_spct(1,2,:))),'--r');
    title(sprintf('GC x2 to x1. Ours: %d, MVGC: %d',t_GC(1,2), GC_x2_to_x1));

    legend({'Ours','MVGC'})
    set(ph2to1, 'linewidth', 3)
    ylim([0 0.2])
    ylabel('spct gc 2 to 1')
    xlabel('Freq')


    fignm = ['./' systems(system_indx).nm '_GCcompareMVGC'];
    savefig(gcf, fignm)
    % saveas(gcf, fignm, 'bmp')
    print(gcf, fignm, '-dpng')

end

if JOB.stochasticint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stochastic interaction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% To obtain stochastic interaction we set all lagged ineractions to zero
% We also do not allow contemporous correlation (force E to be diagonal) 
% no lagged influences
    split_mask_A=eye(N);

% no correlation in the noise covariance 
    split_mask_E=eye(N);

% calc reduced model
    if ~exist(['./' systems(system_indx).nm '_stochinfo.mat'],'file')
        % Both time and frequency domain quantities are estimated
        disp('Computing reduced model parameters')
        
        [S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(G,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    

        % save the quantities 
        save(['./' systems(system_indx).nm '_stochinfo.mat'],'S_r','det_S_r','trace_S_r','prod_diag_S_r','A_r','SIG_r','masked_Delta',...
            'iter_max','gamma','gamma','split_mask_A','split_mask_E','max_order')    

    else

        disp('Loading reduced model parameters')
        load(['./' systems(system_indx).nm '_stochinfo.mat'])    
    end

    % log ratio
    [ratio_S ratio] = ratio_of_dets(S_f, S_r, SIG_f, SIG_r);
    
    % keep results here 
    t_stoch_info = ratio;

% compare against the residuals of the indiv processes, alternative way to
% estimate it
    univar_autocov1=G(1,1,:);
    [univar_A1 univar_SIG1] = autocov_to_var(univar_autocov1);
    univar_autocov2=G(2,2,:);
    [univar_A2 univar_SIG2] = autocov_to_var(univar_autocov2);
    t_stoch_info_2 = log ( univar_SIG1*univar_SIG2  / det (SIG_f) );

%we can also estimate these by summing the spectrums, yet another way to
%estimate
    two_sided_S = cat(3,S_f,S_f(:,:,2:freq_res));
    % We don't care about these here
    two_sided_S(1,2,:) = 0;
    two_sided_S(2,1,:) = 0;

    % sum of log 
    tmpS = squeeze(mean(log(two_sided_S),3));
    log_E22_1 = tmpS(2,2); 
    log_E11_1 = tmpS(1,1);
    t_stoch_info_3=log_E22_1+log_E11_1-log(det (SIG_f));

    % spectral decomp using our method
    t_stoch_info_spct = ratio_S;

    % in the frequency domain we have this equality with  coherence
    C=squeeze((abs(S_f(1,2,:)).^2)./(S_f(1,1,:).*S_f(2,2,:)));
    sdecomp_stoch_info_2 = -log(1-C);


    %plot
    clf
    ph = plot(freqs,real(t_stoch_info_spct),freqs, sdecomp_stoch_info_2,'--r'); 
    legend({'Ours','Coh based'})
    set(ph, 'linewidth', 3)
    legend({'Ours','Coh based'})
    ylabel('spct stoch. int.')
    xlabel('Freq')      
    %time domain quantities
    title(sprintf('Time domain equalities 1) %d 2) %d 3) %d \n',t_stoch_info,t_stoch_info_2,t_stoch_info_3))

    ylim([0 2*t_stoch_info])
    
    fignm = ['./' systems(system_indx).nm '_StochIntCompareCoh'];
    savefig(gcf, fignm)
    % saveas(gcf, fignm, 'bmp')
    print(gcf, fignm, '-dpng')

end

if JOB.predinfo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Predictive information 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% To obtain predictive information all lagged ineractions to zero. We place
% no retrictions on the noise cov

    N = size(G,1);
% no lagged influences
    split_mask_A=zeros(N);

% no correlation in the noise covariance 
    split_mask_E=ones(N);

% calc reduced model
    if ~exist(['./' systems(system_indx).nm '_predinfo.mat'],'file')
        % Both time and frequency domain quantities are estimated
        disp('Computing reduced model parameters')

        [S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(G,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    

        % save the quantities so that we don't have to run this again
        save(['./' systems(system_indx).nm '_predinfo.mat'],'S_r','det_S_r','trace_S_r','prod_diag_S_r','A_r','SIG_r','masked_Delta',...
            'iter_max','gamma','gamma','split_mask_A','split_mask_E','max_order')    

    else

        disp('Loading reduced model parameters')
        load(['./' systems(system_indx).nm '_predinfo.mat'])    
    end

%log ratios
    [ratio_S ratio det_S_f] = ratio_of_dets(S_f, S_r, SIG_f, SIG_r);
    
% keep results here 
    t_pred_info = ratio;
    t_pred_info_spct = ratio_S;

% compare against analytic result ( MI past present )
    t_pred_info2=log( det(Cov_X)/det(SIG_f) );
    t_pred_info_spct2=log( det(Cov_X)./real(det_S_f) );
    
    clf
    ph = plot(freqs,real(t_pred_info_spct),freqs, t_pred_info_spct2,'--r'); 
    set(ph, 'linewidth', 3)
    legend({'Ours','analytic based'})
    ylabel('pred. info.')
    xlabel('Freq')      
    title(sprintf('Time domain equalities 1) %d 2) %d \n',t_pred_info,t_pred_info2))
   
    
    fignm = ['./' systems(system_indx).nm '_predinfoCompareAnalytic'];
    savefig(gcf, fignm)
    % saveas(gcf, fignm, 'bmp')
    print(gcf, fignm, '-dpng')

end  

if JOB.instcaus
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inst caus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To obtain inst information we place no restrictions on the auto regression,
% but force no contemporous correlation

% no restrictions
    split_mask_A=ones(N);

% no correlation in the noise covariance 
    split_mask_E=eye(N);

    % calc reduced model
    if ~exist(['./' systems(system_indx).nm '_instinfo.mat'],'file')
        % Both time and frequency domain quantities are estimated
        disp('Computing reduced model parameters')

        [S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(G,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    

        % save the quantities so that we don't have to run this again
        save(['./' systems(system_indx).nm '_instinfo.mat'],'S_r','det_S_r','trace_S_r','prod_diag_S_r','A_r','SIG_r','masked_Delta',...
            'iter_max','gamma','gamma','split_mask_A','split_mask_E','max_order')    

    else

        disp('Loading reduced model parameters')
        load(['./' systems(system_indx).nm '_instinfo.mat'])    
    end

    %log ratio s
    [ratio_S ratio det_S_f] = ratio_of_dets(S_f, S_r, SIG_f, SIG_r);
    
% keep results here 
    t_inst_info = ratio;
    t_inst_info_spct = ratio_S;

% compare against instantatneou caus, time domain
    t_inst_info2=log(SIG_f(1,1).*SIG_f(2,2)/det(SIG_f));
    
    
     
    % insta caus, ding
    [spct_GC]= autocov_to_spwcgc(G,freq_res);
    sdecomp_GC_x1_to_x2 = squeeze(spct_GC(2,1,:));
    sdecomp_GC_x2_to_x1 = squeeze(spct_GC(1,2,:));
    t_inst_info_spct2=log(squeeze(S_f(1,1,:).*S_f(2,2,:))./det_S_f)-sdecomp_GC_x1_to_x2-sdecomp_GC_x2_to_x1;
    
    
    % can also get inst caus using tf
    H = var2trfun(A,freq_res);
    H_tilde11=squeeze(H(1,1,:)+SIG_f(1,2)/SIG_f(2,2)*H(1,2,:));
    H_tilde22=squeeze(H(2,2,:)+SIG_f(1,2)/SIG_f(1,1)*H(2,1,:));
    t_inst_info_spct23 = log( H_tilde11.*conj(H_tilde11).*SIG_f(1,1) .* SIG_f(2,2).*H_tilde22.*conj(H_tilde22) ./ det_S_f );
    
    
    clf
    ph = plot(freqs,real(t_inst_info_spct),freqs, t_inst_info_spct2,'--r',freqs, t_inst_info_spct2,':k'); 
    legend({'Ours','Inst caus based','Inst caus based2'})
    set(ph, 'linewidth', 3)
    ylabel('inst. info.')
    xlabel('Freq')      
    title(sprintf('Time domain equalities 1) %d 2) %d \n',t_inst_info,t_inst_info2))

    %set reasonable limits in case quantities are very small
    ylim([0 abs(2*t_inst_info)])

    fignm = ['./' systems(system_indx).nm '_instinfoCompareInstCaus'];
    savefig(gcf, fignm)
    % saveas(gcf, fignm, 'bmp')
    print(gcf, fignm, '-dpng')

end  

if JOB.phig
%% Phi G calculation

% in phi g we cut both lagged influences, but place no restriction on the
% noise cov matrix

% no lagged influences
    split_mask_A=eye(N);

% correlation in the noise covariance 
    split_mask_E=ones(N);

% calc reduced model
    if ~exist(['./' systems(system_indx).nm '_phig.mat'],'file')
        % Both time and frequency domain quantities are estimated
        disp('Computing reduced model parameters')
        [S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(G,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    

        % save the quantities so that we don't have to run this again
        save(['./' systems(system_indx).nm '_phig.mat'],'S_r','det_S_r','trace_S_r','prod_diag_S_r','A_r','SIG_r','masked_Delta',...
            'iter_max','gamma','gamma','split_mask_A','split_mask_E','max_order')    

    else

        disp('Loading reduced model parameters')
        load(['./' systems(system_indx).nm '_phig.mat'])    
    end

%log ratio
    [ratio_S ratio det_S_f] = ratio_of_dets(S_f, S_r, SIG_f, SIG_r);
    
% keep results here 
    t_phig = ratio;
    t_phig_spct = ratio_S;
  
    clf
    ph = plot(freqs,real(t_phig_spct)); 
    legend({'Ours'})
    set(ph, 'linewidth', 3)
    ylabel('inst. info.')
    xlabel('Freq')      
    title(sprintf('Time domain 1) %d \n',t_phig))

   
    fignm = ['./' systems(system_indx).nm '_phig'];
    savefig(gcf, fignm)
    % saveas(gcf, fignm, 'bmp')
    print(gcf, fignm, '-dpng')
    
end


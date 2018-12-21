clear all
close all
%%
save_fig_folder = '/Volumes/nfs/drorcohen/Extrinsic_Intrinsic_Simulations/';


systems = [];

%% system 1
num_ext = 2;
num_int = 2;
N = num_ext + num_int;
order = 2;

%% As
A = zeros(N,N,order); 

A_int = zeros(num_int,num_int,order);

A_ext = zeros(num_ext,num_ext,order);

%% driving noise
SIG = eye(N);

%% system 0 - intr/extr
%internal conn
A_int(:,:,1) = [0.9 0;
                0.16 0.8]; % connectivity matrix at time t-1 (Aij = from j to i

A_int(:,:,2) = [-0.5 0;
                -0.2 -0.5]; % connectivity matrix at time t-2

%external conn
A_ext(:,:,1) = [0.2 0;
                0 0.3]; % connectivity matrix at time t-1 (Aij = from j to i

A_ext(:,:,2) = [-0.1 0;
                0 -0.2]; % connectivity matrix at time t-2

A(1:num_int,1:num_int,:) = A_int;
A(end-num_ext+1:end,end-num_ext+1:end,:) = A_ext;

%% computation

[G,info] = var_to_autocov(A,SIG);

GC_graph = zeros(N);
for causies = 1:N
    for causals = 1:N
        if causies ~= causals;
            F = autocov_to_mvgc(G, causies, causals);
            GC_graph(causies,causals) = F;
        end
    end
end

GC_graph(GC_graph<1e-12) = 0;
GC_graph_obj = digraph(GC_graph);
GC_graph_obj.Edges
nodelabels = {'internal 1','internal 2','external 1','external 2'};
h = plot(GC_graph_obj,'EdgeLabel',GC_graph_obj.Edges.Weight,'NodeLabel',nodelabels,...
    'ArrowSize',30, 'LineWidth', 50*GC_graph_obj.Edges.Weight,'MarkerSize',20);

nl = h.NodeLabel;
h.NodeLabel = '';
xd = get(h, 'XData');
yd = get(h, 'YData');
text(xd, yd, nl, 'FontSize',20, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')


%



%% system 2 - uni directional no inst;
SIG_f = [1 0; 0 0.7];

autocov_to_mvgc
% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.9 0;
    0.16 0.8]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.5 0;
    -0.2 -0.5]; % connectivity matrix at time t-2

systems(2).A = A;
systems(2).SIG_f = SIG_f;
systems(2).nm='unidir_no_inst';



%% system 3 - unidir with inst
SIG_f = [1 0.1; 0.1 0.7];
normalizer = sqrt(det(SIG_f));
SIG_f = SIG_f/normalizer;


% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.9 0;
    0.16 0.8]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.5 0;
    -0.2 -0.5]; % connectivity matrix at time t-2


systems(3).A = A;
systems(3).SIG_f = SIG_f;
systems(3).nm='unidir_with_inst';

%% system 4 - bidir no inst
SIG_f = [1 0; 0 0.7];

% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.9 0.2;
    0.16 0.8]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.5 0;
    -0.2 -0.5]; % connectivity matrix at time t-2


systems(4).A = A;
systems(4).SIG_f = SIG_f;
systems(4).nm='bidir_no_inst';

%% choose system you want to investigate
system_indx = 4;
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
JOB.stochasticint=0;

% Compute predictive information
JOB.predinfo=0;

% Compute predictive information
JOB.instcaus=0;

% Compute phi
JOB.phig=0;

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
    t_GC_spct(2,1,:) = ratio_S; %x1 to x2
    
    % for checks later
    pcUs1=SIG_r(1,1)-SIG_r(1,2)/SIG_r(2,2)*SIG_r(2,1);
    pcUs2=SIG_r(2,2)-SIG_r(1,2)/SIG_r(1,1)*SIG_r(2,1);
    S_r_dir1 = S_r 
    

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
    t_GC_spct(1,2,:) = ratio_S; %x2 to x1

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
    
    
    %% Check freq domain equivalence
    %%Check Sxx in full and reduced model 
    S_r_dir2 = S_r; %(x2 to x1)
    
    clf
    subplot(2,1,1)
    plot(squeeze(S_f(1,1,:)))
    hold on
    plot(squeeze(real(S_r_dir2(1,1,:))),'--')
    
    
    
    subplot(2,1,2)
    plot(squeeze(S_f(2,2,:)))
    hold on
    plot(squeeze(real(S_r_dir1(2,2,:))),'--')


     
    %% now check the partial spectrums
    
     
    Ps2_f = S_f(2,2,:) - S_f(1,2,:) ./ S_f(1,1,:) .* S_f(2,1,:);
    Ps1_f = S_f(1,1,:) - S_f(1,2,:) ./ S_f(2,2,:) .* S_f(2,1,:);
    
    % The first way to express is to factorize with the power leading the
    % terms, since the power vanishes, GC is a function of the partial
    % spcts
    
    Ps2_r = S_r_dir2(2,2,:) - S_r_dir2(1,2,:) ./ S_r_dir2(1,1,:) .* S_r_dir2(2,1,:);      
    Ps1_r = S_r_dir1(1,1,:) - S_r_dir1(1,2,:) ./ S_r_dir1(2,2,:) .* S_r_dir1(2,1,:);
    
   %we have GC in %(x2 to x1)
    log(S_f(1,1,:) ./ S_r_dir2(1,1,:)) + log(real(Ps2_r./Ps2_f))
    %but first term is zero, so
    GC_dir2 = log(real(Ps2_r./Ps2_f));
    %same for other dir
    GC_dir1 = log(real(Ps1_r./Ps1_f));
    
    
    %% a second way is to factorize is in the opposite way
    Ps2_r2 = S_r_dir2(1,1,:) - S_r_dir2(1,2,:) ./ S_r_dir2(2,2,:) .* S_r_dir2(2,1,:);
    Ps1_r2 = S_r_dir1(2,2,:) - S_r_dir1(1,2,:) ./ S_r_dir1(1,1,:) .* S_r_dir1(2,1,:);
    
    %we have GC in %(x2 to x1)
    GC_dir2_2 = log(S_r_dir2(2,2,:) ./ S_f(2,2,:)) + log(real(Ps2_r2./Ps1_f));
        
    %we have GC in %(x1 to x2)
    GC_dir1_1 = log(S_r_dir1(1,1,:) ./ S_f(1,1,:)) + log(real(Ps1_r2./Ps2_f));

    
    clf
    subplot(2,1,1)
    plot(squeeze(GC_dir1))
    hold on
    plot(sdecomp_GC_x1_to_x2, '--')
    plot(real(squeeze(GC_dir1_1)), '--')

    subplot(2,1,2)
    plot(squeeze(GC_dir2))
    hold on
    plot(sdecomp_GC_x2_to_x1, '--')
    plot(real(squeeze(GC_dir2_2)), '--')
    
    %% check the ration of other spectrums  - these are rubbish   
% %     subplot(2,1,1)
% %     plot(log(squeeze(S_f(1,1,:))./squeeze(real(S_r_dir1(1,1,:))) ))
% %     hold on
% %     
% %        
% %     subplot(2,1,2)
% %     plot(log(squeeze(S_f(2,2,:))./squeeze(real(S_r_dir2(2,2,:))) ))
    
    %% finally, check the transer func equiv    
    H = var2trfun(A,freq_res);
    H_tilde11=squeeze(H(1,1,:)+SIG_f(1,2)/SIG_f(1,1)*H(1,2,:));
    H_tilde22=squeeze(H(2,2,:)+SIG_f(1,2)/SIG_f(2,2)*H(2,1,:));
    
    Sxx_1 = H_tilde11 .* conj(H_tilde11) .* SIG_f(1,1) + squeeze(H(1,2,:) .* SIG_f(2,2) .* conj(H(1,2,:)));
    
    clf
    plot(real(Sxx_1))
    hold on
    plot(squeeze(S_f(1,1,:)))
    
    clf
    plot(H_tilde11 .* conj(H_tilde11) .* SIG_f(1,1) )
    hold on
    plot(real(squeeze(squeeze(H(1,2,:) .* SIG_f(2,2) .* conj(H(1,2,:))))));
    hold on
    plot(real((squeeze(Ps1_f))));
    plot(real(Sxx_1)-real((squeeze(Ps1_f))));

    
    GC_x2_to_x1_tf =  squeeze(S_f(1,1,:)) ./ (H_tilde11 .* conj(H_tilde11) .* SIG_f(1,1));
    GC_x1_to_x2_tf =  squeeze(S_f(2,2,:)) ./ (H_tilde22 .* conj(H_tilde22) .* SIG_f(2,2));
    
    
    % check num equivalence
    %x2 to x1
    plot(real(squeeze(S_f(1,1,:) - Ps2_r)))
    % check denom equiv
    plot(real(squeeze(H_tilde11 .* conj(H_tilde11) .* SIG_f(1,1)) - squeeze(Ps2_f)))
    
    %check partialization equivalence
    clf
    plot(real(squeeze(H_tilde11 .* conj(H_tilde11) .* SIG_f(1,1))) - real(squeeze(Ps1_f)))
%     plot(real(squeeze(H_tilde11 .* conj(H_tilde11) .* SIG_f(1,1))) - real(squeeze(Ps1_f)))
%     plot(real(squeeze(H_tilde11 .* conj(H_tilde11) .* SIG_f(1,1))) - real(squeeze(Ps1_f)))

    
    %check ratio
    clf
    plot(real(GC_x2_to_x1_tf) - real(squeeze(Ps2_r./Ps2_f)))

    
    
    %
    clf
    subplot(2,1,1)
    plot(sdecomp_GC_x1_to_x2, '--')
    hold on
    plot(log(GC_x1_to_x2_tf),'s')
    
    subplot(2,1,2)
    plot(sdecomp_GC_x2_to_x1, '--')
    hold on
    plot(log(GC_x2_to_x1_tf),'s')
    

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
    H_tilde11=squeeze(H(1,1,:)+SIG_f(1,2)/SIG_f(1,1)*H(1,2,:));
    H_tilde22=squeeze(H(2,2,:)+SIG_f(1,2)/SIG_f(2,2)*H(2,1,:));
    
    t_inst_info_spct23 = log( H_tilde11.*conj(H_tilde11).*SIG_f(1,1) .* SIG_f(2,2).*H_tilde22.*conj(H_tilde22) ./ det_S_f );
    mean(t_inst_info_spct23)
    
    clf
    ph = plot(freqs,real(t_inst_info_spct),freqs, real(t_inst_info_spct2),'--r',freqs, real(t_inst_info_spct2),':k',freqs,real(t_inst_info_spct23),'s'); 
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


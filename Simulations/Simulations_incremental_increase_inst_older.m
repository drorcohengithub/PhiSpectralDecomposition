%% you need to seriously improve this for before making the github repo

clear all
close all

%% load systems
systems = Example_systems();


%% choose system you want to investigate
system_indx = 2;
A = systems(system_indx).A;
SIG = systems(system_indx).SIG; 

% Arbitrarily set the freq resolution and samp rate
freq_res  = 100;
samp_rate = 1000;

 % Obtain the frequency bins
freqs = sfreqs(freq_res,samp_rate);

%% we will calculate the quantities for increasing levels of residual noise correlation, which increases instantaneous interaction
% minimum level
min_inst = 0;
% maximum level, makes sure we remain positive definint
max_inst = min(diag(SIG)) - 0.05;
% numver of different levels
num_levels = 10;
% the actual levels
inst_levels = linspace(min_inst,max_inst,num_levels);


%% keep results here
incremental_inst_increase_results = [];

%start looping
for level_indx = 1:num_levels
   
    %Get the original noise correlations
    this_SIG = SIG;
   
    %this inst level
    this_inst = inst_levels(level_indx);
    
    % recreaste the residual cov matrix
    this_SIG(1,2) = this_inst;
    this_SIG(2,1) = this_inst;
    
    %store it here for later
    incremental_inst_increase_results(level_indx).this_SIG = this_SIG;
    
    % get the corresponing autocov function 
    [X,info] = var_to_autocov(A,this_SIG);
    
    % store it, and any info
    incremental_inst_increase_results(level_indx).X = X;
    incremental_inst_increase_results(level_indx).info = info;
        
    % Check that everything is ok
    var_info(info,true);  

    % We will need these 
    % The cov of X
    Cov_X  = X(:,:,1);
    
    % The autocov without cov X
    Cov_XY = X(:,:,2:end);
    
    %max order of the disconnected model
    max_order=size(Cov_XY,3);
    
    % store it
    incremental_inst_increase_results(level_indx).max_order = max_order;

    % The spectral density matrix of the full model
    [S] = autocov_to_cpsd(X,freq_res);
    %store it
    incremental_inst_increase_results(level_indx).S_f = S;
    
    %get spectral GC
    [spct_GC]= autocov_to_spwcgc(X,freq_res);
    spct_GC_x1_to_x2 = squeeze(spct_GC(2,1,:));
    spct_GC_x2_to_x1 = squeeze(spct_GC(1,2,:));
    
    % time domain
    [GC]= autocov_to_pwcgc(X);
    GC_x1_to_x2 = squeeze(GC(2,1,:));
    GC_x2_to_x1 = squeeze(GC(1,2,:));
    
    % store them
    incremental_inst_increase_results(level_indx).spct_GC_x1_to_x2 = spct_GC_x1_to_x2;
    incremental_inst_increase_results(level_indx).spct_GC_x2_to_x1 = spct_GC_x2_to_x1;

    incremental_inst_increase_results(level_indx).GC_x1_to_x2 = GC_x1_to_x2;
    incremental_inst_increase_results(level_indx).GC_x2_to_x1 = GC_x2_to_x1;

    
    % stoch int
    C=squeeze((abs(S(1,2,:)).^2)./(S(1,1,:).*S(2,2,:)));
    spct_stoch_int = -log(1-C);
    incremental_inst_increase_results(level_indx).spct_stoch_int = spct_stoch_int;
    
    univar_autocov1=X(1,1,:);
    [univar_A1 univar_SIG1] = autocov_to_var(univar_autocov1);
    univar_autocov2=X(2,2,:);
    [univar_A2 univar_SIG2] = autocov_to_var(univar_autocov2);
    stoch_int = log ( univar_SIG1*univar_SIG2  / det (this_SIG) );
    incremental_inst_increase_results(level_indx).stoch_int = stoch_int;
    
    %pred info
    det_S_f=[];
    for freqIndx=1:size(S,3)
        det_S_f(freqIndx)=det(squeeze(S(:,:,freqIndx)));
    end
    spct_pred_info=log( det(Cov_X)./real(det_S_f) );
    incremental_inst_increase_results(level_indx).spct_pred_info = spct_pred_info;

    pred_info=log( det(Cov_X)/det(this_SIG) );
    incremental_inst_increase_results(level_indx).pred_info = pred_info;
    
    %inst
    spct_inst_intr=log(this_SIG(1,1).*this_SIG(2,2)/det(this_SIG));
    incremental_inst_increase_results(level_indx).spct_inst_intr=spct_inst_intr;
    
          
%         t_inst_info_spct2=log(squeeze(S_f(1,1,:).*S_f(2,2,:))./det_S_f)-sdecomp_GC_x1_to_x2-sdecomp_GC_x2_to_x1;

    % no lagged influences
    split_mask_A=eye(N);

    % correlation in the noise covariance 
    split_mask_E=ones(N);

    % Both time and frequency domain quantities are estimated
    disp('Computing reduced model parameters')
    [S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(X,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    

 
    %log ratio
    [ratio_S ratio det_S_f] = ratio_of_dets(S, S_r, this_SIG, SIG_r);
    
    % keep results here 
    t_phig = ratio;
    incremental_inst_increase_results(level_indx).spct_phi = ratio_S;
    incremental_inst_increase_results(level_indx).phi = t_phig;
    incremental_inst_increase_results(level_indx).phi_S_r = S_r;
    incremental_inst_increase_results(level_indx).phi_A_r = A_r;
    incremental_inst_increase_results(level_indx).phi_SIG_r = SIG_r;

    
   
end

%%
params_desc = ['_min' num2str(min_inst) 'max' num2str(max_inst) 'nsteps' num2str(num_levels)];

save(['./' systems(system_indx).nm '_incremental_inst_increase_res' params_desc '.mat'],'incremental_inst_increase_results','masked_Delta',...
            'iter_max','gamma','split_mask_A','split_mask_E','max_order', 'freqs') 
        
%% also get model params for gc
        
GC_incremental_inst_increase_results = [];
for level_indx = 1:num_levels
   
    % Get all the characteristics of the system
    %get the autocov
    this_inst = inst_levels(level_indx);
    this_SIG = SIG_f;
    
    this_SIG(1,2) = this_inst;
    this_SIG(2,1) = this_inst;
    
    GC_incremental_inst_increase_results(level_indx).this_SIG = this_SIG;
    
    [X,info] = var_to_autocov(A,this_SIG);
    GC_incremental_inst_increase_results(level_indx).G = X;
    GC_incremental_inst_increase_results(level_indx).info = info;
        
    N = size(X,1);

    % Check that everything is ok
    var_info(info,true);  

    % We will need these 
    % The cov of X
    Cov_X  = X(:,:,1);
    % The auto cov without cov X
    Cov_XY = X(:,:,2:end);
    
    %max ord reduced
    max_order=size(Cov_XY,3);
    GC_incremental_inst_increase_results(level_indx).max_order = max_order;

    % The spectral density matrix of the full model
    [S] = autocov_to_cpsd(X,freq_res);
    GC_incremental_inst_increase_results(level_indx).S_f = S;
    
    %spct GC
    
    % no lagged influences
    split_mask_A=ones(N);
    split_mask_A(2,1) = 0;

    % correlation in the noise covariance 
    split_mask_E=ones(N);

    % Both time and frequency domain quantities are estimated
    disp('Computing reduced model parameters')
    [S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(X,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    

 
    GC_incremental_inst_increase_results(level_indx).S_r = S_r;
    GC_incremental_inst_increase_results(level_indx).det_S_r = det_S_r;
    GC_incremental_inst_increase_results(level_indx).trace_S_r = trace_S_r;
    GC_incremental_inst_increase_results(level_indx).prod_diag_S_r = prod_diag_S_r; 
    GC_incremental_inst_increase_results(level_indx).A_r = A_r;
    GC_incremental_inst_increase_results(level_indx).SIG_r = SIG_r;
    %log ratio
    [ratio_S ratio det_S_f] = ratio_of_dets(S, S_r, this_SIG, SIG_r);
    
    % keep results here 
    t_phig = ratio;
    GC_incremental_inst_increase_results(level_indx).spct_phi = ratio_S;
    
   
end

%%
params_desc = ['_min' num2str(min_inst) 'max' num2str(max_inst) 'nsteps' num2str(num_levels)];

save(['./' systems(system_indx).nm 'GC_incremental_inst_increase_res' params_desc '.mat'],'GC_incremental_inst_increase_results','masked_Delta',...
            'iter_max','gamma','split_mask_A','split_mask_E','max_order', 'freqs','det_S_f') 
        
%%
H = var2trfun(A,freq_res);
H_r = var2trfun(A_r, freq_res);

%%
H(1,2,:)
H_r(1,2,:)

        


function [ results ] = calc_measures(A,SIG,measure_names,check,freq_res,samp_rate,iter_max,gamma,min_error)


%Calculate the various measures from the paper
% Only for bivariate systems
% A - connectivity matrix of full system
% Sig - covariance of residuals of full system
% measure_names - the names of the measures you want to calc; 'phi','gc',
% 'stoch_int','inst_int','pred_info'. Use 'all' for calculating all
% check - boolean, whether or not to check if the various ways of calculating the quantities are identical. Takes longer
% freq_res - frequency resolution, default 100Hz
% samp_rate- sampling rate, default 1000Hz 

%% only for phi params
% iter_max - maximum iteration number for finding A_p, default 12000
% gamma - initial step size for optimization 
% min_error - error thresh at which to stop optimizing 

%% check inputs
if strcmp(lower(measure_names),'all')
    measure_names = {'phi','gc','stoch_int','inst_int','pred_info'};
end

if nargin < 4 || isempty(check)
    
    check = 1;
    
end

if nargin < 5 || isempty(freq_res)
    freq_res = 100;
end
    

if nargin < 6 || isempty(samp_rate)
    samp_rate = 1000;
end

if nargin < 7 || isempty(iter_max)
    
    iter_max = 12000;
    
end

if nargin <= 8 || isempty(gamma)
    
    gamma = 0.1;
    
end

if nargin <= 9 || isempty(min_error)
    
    min_error = 1e-14;
    
end

%% some basics
% Get all the characteristics of the system
[X,info] = var_to_autocov(A,SIG);

% Check that everything is ok
var_info(info,true);

% Obtain the frequency bins
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
[S] = autocov_to_cpsd(X,freq_res);


% store the results here. 
results = [];
results.Cov_XY = Cov_XY;
results.Cov_X = Cov_X;
results.A = A;
results.SIG = SIG;
results.S = S;
results.freq_res = freq_res;
results.samp_rate = samp_rate;
results.freqs = freqs;

for measure_indx = 1 : length (measure_names)
    
    switch measure_names{measure_indx}
        
        case 'phi'
            
            % no lagged influences
            split_mask_A=eye(N);
            results.Phi_model.split_mask_A = split_mask_A;

            % correlation in the noise covariance 
            split_mask_E=ones(N);
            results.Phi_model.split_mask_E = split_mask_E;
            % Both time and frequency domain quantities are estimated
            disp('Computing reduced model parameters')
            
            [results.Phi_model.S_r,...
                results.Phi_model.det_S_r,...
                results.Phi_model.trace_S_r,...
                results.Phi_model.prod_diag_S_r,...
                results.Phi_model.A_r,...
                results.Phi_model.SIG_r,...
                results.Phi_model.masked_Delta] = get_reduced_S_from_autoCov(X,...
                split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    

             %log ratio
            [sdecomp_phi tmp results.det_S] = ratio_of_dets(S, results.Phi_model.S_r, results.SIG, results.Phi_model.SIG_r);
            results.phi =1/2 * tmp;
            % may not be quite real
            results.sdecomp_phi = 1/2 * real(sdecomp_phi);
        
        case 'gc'
            
            [t_GC_mvgc]= autocov_to_pwcgc(X);
            %note factor of half
            results.GC_x1_to_x2 = 1/2 *t_GC_mvgc(2,1);
            results.GC_x2_to_x1 = 1/2 *t_GC_mvgc(1,2);

            [spct_GC]= autocov_to_spwcgc(X,freq_res);
            results.sdecomp_GC_x1_to_x2 =1/2 * squeeze(spct_GC(2,1,:));
            results.sdecomp_GC_x2_to_x1 = 1/2 *squeeze(spct_GC(1,2,:));
            
            
            
            %% we need these to get the details of the model
            %first dir
            split_mask_A=ones(N);
            split_mask_A(2,1) = 0;
            results.GC_x2_to_x1_model.split_mask_A = split_mask_A;

            % correlation in the noise covariance 
            split_mask_E=ones(N);
            results.GC_x2_to_x1_model.split_mask_E = split_mask_E;
            % Both time and frequency domain quantities are estimated
            disp('Computing reduced model parameters')
            
            [results.GC_x2_to_x1_model.S_r,...
                results.GC_x2_to_x1_model.det_S_r,...
                results.GC_x2_to_x1_model.trace_S_r,...
                results.GC_x2_to_x1_model.prod_diag_S_r,...
                results.GC_x2_to_x1_model.A_r,...
                results.GC_x2_to_x1_model.SIG_r,...
                results.GC_x2_to_x1_model.masked_Delta] = get_reduced_S_from_autoCov(X,...
                split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);
        
            %  
            split_mask_A=ones(N);
            split_mask_A(1,2) = 0;
            results.GC_x1_to_x2_model.split_mask_A = split_mask_A;

            % correlation in the noise covariance 
            split_mask_E=ones(N);
            results.GC_x1_to_x2_model.split_mask_E = split_mask_E;
            % Both time and frequency domain quantities are estimated
            disp('Computing reduced model parameters')
            
            [results.GC_x1_to_x2_model.S_r,...
                results.GC_x1_to_x2_model.det_S_r,...
                results.GC_x1_to_x2_model.trace_S_r,...
                results.GC_x1_to_x2_model.prod_diag_S_r,...
                results.GC_x1_to_x2_model.A_r,...
                results.GC_x1_to_x2_model.SIG_r,...
                results.GC_x1_to_x2_model.masked_Delta] ....
                = get_reduced_S_from_autoCov(X,...
                split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);
            
            
        case 'stoch_int'
            
            
             % stoch int
            C=squeeze((abs(S(1,2,:)).^2)./(S(1,1,:).*S(2,2,:)));
            results.sdecomp_stoch_int = 1/2 * -log(1-C);

            results.univar_autocov1=X(1,1,:);
            [results.univar_A1 results.univar_SIG1] = autocov_to_var(results.univar_autocov1);
            results.univar_autocov2=X(2,2,:);
            [results.univar_A2 results.univar_SIG2] = autocov_to_var(results.univar_autocov2);
            results.stoch_int = 1/2 * log ( results.univar_SIG1*results.univar_SIG2  / det (results.SIG) );
            
            A_r = [];
            A_r(1,1,:) = results.univar_A1;
            A_r(2,2,:) = results.univar_A2;
            
            results.A_r = A_r;
            
            SIG_r = [];
            SIG_r = [results.univar_SIG1,0;0,results.univar_SIG2];
            
            results.SIG_r = SIG_r;

                [t_GC_mvgc]= autocov_to_pwcgc(X);
            results.GC_x1_to_x2 = 1/2 * t_GC_mvgc(2,1);
            results.GC_x2_to_x1 = 1/2 * t_GC_mvgc(1,2);

            [spct_GC]= autocov_to_spwcgc(X,freq_res);
            results.sdecomp_GC_x1_to_x2 = 1/2 * squeeze(spct_GC(2,1,:));
            results.sdecomp_GC_x2_to_x1 = 1/2 * squeeze(spct_GC(1,2,:));
            
            S_r = S;
            S_r(1,2,:) = S_r(1,2,:)*0;
            S_r(2,1,:) = S_r(2,1,:)*0;
            
            
           results.stoch_int_model.S_r = S_r;
           results.stoch_int_model.log_det_S_r = log_of_dets(S_r);
            
                                                         
            
        case 'inst_int'
            
           %inst
            results.inst_int = 1/2 *log(SIG(1,1).*SIG(2,2)/det(SIG));
            results.sdecomp_inst_int = ones(length(freqs),1) * results.inst_int;
            
%             [results.Phi_model.S_r,results.Phi_model.det_S_r,results.Phi_model.trace_S_r,...
%                 results.Phi_model.prod_diag_S_r,results.Phi_model.A_r,results.Phi_model.SIG_r,...
%                 results.Phi_model.masked_Delta]
            
        
        case 'pred_info'

            results.det_S=[];
            for freqIndx=1:size(S,3)
                results.det_S(freqIndx)=det(squeeze(S(:,:,freqIndx)));
            end
            
            results.sdecomp_pred_info=1/2 *log( det(Cov_X)./real(results.det_S) )';
            results.pred_info=1/2 *log( det(Cov_X)/det(SIG) );
            
%             [results.Phi_model.S_r,results.Phi_model.det_S_r,results.Phi_model.trace_S_r,...
%                 results.Phi_model.prod_diag_S_r,results.Phi_model.A_r,results.Phi_model.SIG_r,...
%                 results.Phi_model.masked_Delta]
       
        otherwise
            warning(['unknown measure ' measure_names{measure_indx} ', skipping']) 
    end


end

end


    

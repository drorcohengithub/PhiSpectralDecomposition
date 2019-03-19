clear all
close all

%% load systems
systems = Example_systems();

%% choose system you want to investigate
for system_indx = 1:length(systems)
%     system_indx =  1;
    A = systems(system_indx).A;
    SIG = systems(system_indx).SIG; 


    % Arbitrarily set the freq resolution and samp rate
    freq_res  = 100;
    samp_rate = 1000;

    % 
    % % We will need these 
    % % The cov of X
    % Cov_X  = G(:,:,1);
    % % The auto cov without cov X
    % Cov_XY = G(:,:,2:end);
    % 
    % % We need to choose the maximum lag we will calculate the reduced model up to.
    % % Note that for good estimate the maximum lag of the reduced model will be
    % % much larger than that of the full model. A safe but potentially
    % % over-generous strategy is to use the max lag of the autocov function
    % max_order=size(Cov_XY,3);
    % 


    % What to do
    % Compute GC (both directions)
    JOB.gc=1;

    % Compute stochastic interaction
    JOB.stoch_int=1;

    % Compute predictive information
    JOB.pred_info=1;

    % Compute predictive information
    JOB.inst_int=1;

    % Compute phi
    JOB.phig=1;

    % Some parameters for the optimization. These seem to produced reasonable
    % results for all metrics and simulations, but you can probably get better results 
    % by adjusting on a case by case basis 
    %max iter size
    % iter_max = 12000;
    % %learning rate, this is adjusted automatically if convergence stagnates
    % gamma = 0.1;
    % %min error, stop optimizing if reached
    % min_error = 1e-10;
    % %dimension of system
    % N = size(G,1);


    %% Start the work
    if JOB.gc


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Granger causality
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % obtain in the standard way
        measure_names = {'gc'};
        results = calc_measures(A,SIG,measure_names);

            % save the quantities so that we don't have to run this again
        save(['./' systems(system_indx).nm '_GC.mat'],'results')


    end

    if JOB.stoch_int
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stochastic interaction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % obtain in the standard way
        measure_names = {'stoch_int'};
        results = calc_measures(A,SIG,measure_names);

        % save the quantities so that we don't have to run this again
        save(['./' systems(system_indx).nm '_stoch_int.mat'],'results')


    end

    if JOB.pred_info
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Predictive information 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % obtain in the standard way
        measure_names = {'pred_info'};
        results = calc_measures(A,SIG,measure_names);

        % save the quantities so that we don't have to run this again
        save(['./' systems(system_indx).nm '_pred_info.mat'],'results')

    end  

    if JOB.inst_int
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Instantaneous interaction (aka instantaneous causality)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % obtain in the standard way
        measure_names = {'inst_int'};
        results = calc_measures(A,SIG,measure_names);

        % save the quantities so that we don't have to run this again
        save(['./' systems(system_indx).nm '_inst_int.mat'],'results')
    end  

    if JOB.phig
    %% Phi G calculation
    % obtain in the standard way
        measure_names = {'phi'};
        results = calc_measures(A,SIG,measure_names);

        % save the quantities so that we don't have to run this again
        save(['./' systems(system_indx).nm '_phi.mat'],'results')

    end
end

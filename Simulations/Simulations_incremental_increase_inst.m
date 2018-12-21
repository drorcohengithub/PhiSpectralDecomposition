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

% calculate all the measures
measure_names = 'all';

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
    results = calc_measures(A,this_SIG,measure_names,[],freq_res,samp_rate);
    
    incremental_inst_increase_results(level_indx).results = results;
   
   
   
end

%%
params_desc = ['_min' num2str(min_inst) 'max' num2str(max_inst) 'nsteps' num2str(num_levels)];
save(['./' systems(system_indx).nm '_incremental_inst_increase_res' params_desc '.mat'],'incremental_inst_increase_results') 
        
%% calculate inst interaction and Granger causality in reduced models

incremental_inst_increase_results = load(['./' systems(system_indx).nm '_incremental_inst_increase_res' params_desc '.mat']);
incremental_inst_increase_results = incremental_inst_increase_results.incremental_inst_increase_results;

incremental_inst_phi_reduced_model_results = [];
results_red_model_gc = [];
results_red_model_stoch = [];

%%
for level_indx = 1:num_levels
   
    
    %
    this_result = incremental_inst_increase_results(level_indx).results;

    %Get params of the phi reduced model
    this_SIG = this_result.Phi_model.SIG_r;
    this_A_r = this_result.Phi_model.A_r;

    
    incremental_inst_phi_reduced_model_results(level_indx).results = calc_measures(this_A_r,this_SIG,measure_names,[],freq_res,samp_rate);
    
    
    %Get params of the phi reduced model
    this_SIG = this_result.GC_x1_to_x2_model.SIG_r;
    this_A_r = this_result.GC_x1_to_x2_model.A_r;

    
    incremental_inst_GC_x1_to_x2_model_reduced_model_results(level_indx).results = calc_measures(this_A_r,this_SIG,measure_names,[],freq_res,samp_rate);
    
       
    
    
    this_SIG = this_result.GC_x2_to_x1_model.SIG_r;
    this_A_r = this_result.GC_x2_to_x1_model.A_r;

    
    incremental_inst_GC_x2_to_x1_model_reduced_model_results(level_indx).results = calc_measures(this_A_r,this_SIG,measure_names,[],freq_res,samp_rate);
    
    
     
   
   
end





save(['./' systems(system_indx).nm 'incremental_inst_phi_reduced_model_results' params_desc '.mat'],'incremental_inst_phi_reduced_model_results') 

save(['./' systems(system_indx).nm 'incremental_inst_GC_x1_to_x2_model_reduced_model_results' params_desc '.mat'],'incremental_inst_GC_x1_to_x2_model_reduced_model_results') 

save(['./' systems(system_indx).nm 'incremental_inst_GC_x2_to_x1_model_reduced_model_results' params_desc '.mat'],'incremental_inst_GC_x2_to_x1_model_reduced_model_results') 











        


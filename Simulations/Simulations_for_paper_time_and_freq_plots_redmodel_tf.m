clear all
close all

%% load sys
systems = Example_systems();

%% choose system you want to investigate
system_indx = 3;
system_nm = systems(system_indx).nm;

if strcmp(system_nm,'unidir_with_inst')
    ylims = [-5 1;
            -2 0;
            -2.5 2;
            -0.05 2];
elseif strcmp(system_nm,'bidir_with_inst')
    
    ylims = [-4 2;
        -1 0;
        -2 2;
        -0.05 1];
end

%% load

GC_data = load(['./' system_nm '_GC.mat']); 
GC_x2_to_x1 = GC_data.results.GC_x2_to_x1;
sdecomp_GC_x2_to_x1 = GC_data.results.sdecomp_GC_x2_to_x1;

GC_x1_to_x2 = GC_data.results.GC_x1_to_x2;
sdecomp_GC_x1_to_x2 = GC_data.results.sdecomp_GC_x1_to_x2;


stoch_int_data = load(['./' system_nm '_stoch_int.mat']);
stoch_int = stoch_int_data.results.stoch_int;
sdecomp_stoch_int = stoch_int_data.results.sdecomp_stoch_int;
% 
% 
% predinfo_data = load(['./' system_nm '_pred_info.mat']);
% pred_info = predinfo_data.results.pred_info;
% sdecomp_pred_info = predinfo_data.results.sdecomp_pred_info;
% 
% 
% inst_int_data = load(['./' system_nm '_inst_int.mat']);   
% inst_int = inst_int_data.results.inst_int;
% sdecomp_inst_int = inst_int_data.results.sdecomp_inst_int;


Phig_data = load(['./' system_nm '_phi.mat']);
% phig.results.
Phi_G = Phig_data.results.phi;
sdecomp_Phi_G = Phig_data.results.sdecomp_phi;

% these are the same for all measures
freqs = Phig_data.results.freqs;

%% orig sys

A = systems(system_indx).A;
SIG = systems(system_indx).SIG;
freq_res = GC_data.results.freq_res;
tf = var2trfun(A,freq_res);
S_f = var_to_cpsd(A,SIG,freq_res);
log_det_S_f = log_of_dets(S_f);

plt_transfer_func_decomp( tf,SIG, freqs, log_det_S_f, ylims )

savefig(['./' system_nm '_tfdecomp_orig.fig'])
saveas(gcf, ['./' system_nm '_tfdecomp_orig.svg'],'svg')

%% GCs
A_r = GC_data.results.GC_x2_to_x1_model.A_r;
SIG_r = GC_data.results.GC_x2_to_x1_model.SIG_r;
freq_res = GC_data.results.freq_res;
tf = var2trfun(A_r,freq_res);

plt_transfer_func_decomp( tf,SIG_r, freqs, log_det_S_f, ylims )

savefig(['./' system_nm '_tfdecomp_redmodelGC_x2_to_x1.fig'])
saveas(gcf, ['./' system_nm '_tfdecomp_redmodelGC_x2_to_x1.svg'],'svg')

%%
A_r = GC_data.results.GC_x1_to_x2_model.A_r;
SIG_r = GC_data.results.GC_x1_to_x2_model.SIG_r;
freq_res = GC_data.results.freq_res;
tf = var2trfun(A_r,freq_res);

plt_transfer_func_decomp( tf,SIG_r, freqs, log_det_S_f, ylims )

savefig(['./' system_nm '_tfdecomp_redmodelGC_x1_to_x2.fig'])
saveas(gcf, ['./' system_nm '_tfdecomp_redmodelGC_x1_to_x2.svg'],'svg')

%% stoch int

A_r = stoch_int_data.results.A_r;
SIG_r = stoch_int_data.results.SIG_r;
freq_res = stoch_int_data.results.freq_res;
tf = var2trfun(A_r,freq_res);

plt_transfer_func_decomp( tf,SIG_r, freqs, log_det_S_f, ylims )

savefig(['./' system_nm '_tfdecomp_redmodelStochintPhi.fig'])
saveas(gcf, ['./' system_nm '_tfdecomp_redmodelStochint.svg'],'svg')

%% phi

A_r = Phig_data.results.Phi_model.A_r;
SIG_r = Phig_data.results.Phi_model.SIG_r;
freq_res = GC_data.results.freq_res;
tf = var2trfun(A_r,freq_res);

plt_transfer_func_decomp( tf,SIG_r, freqs, log_det_S_f, ylims )

savefig(['./' system_nm '_tfdecomp_redmodelPhi.fig'])
saveas(gcf, ['./' system_nm '_tfdecomp_redmodelPhi.svg'],'svg')



%% stoch int
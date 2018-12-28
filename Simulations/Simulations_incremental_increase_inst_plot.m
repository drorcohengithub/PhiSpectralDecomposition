
clear all
close all

%% load results
fl_desc = 'unidir_no_inst_incremental_inst_increase_res_min0max0.65nsteps10';
% fl_desc = 'bidir_no_inst_incremental_inst_increase_res_min0max0.65nsteps10';
load(['./' fl_desc '.mat']) 

%%
inst_ints = [];
gcs_x1tox2 = [];
gcs_x2tox1 = [];
pred_infos = [];
stoch_ints = [];
phis = [];

for inst_inc_level = 1:length(incremental_inst_increase_results)
    this_results = incremental_inst_increase_results(inst_inc_level).results;
    inst_ints(inst_inc_level) = this_results.inst_int;
    gcs_x1tox2(inst_inc_level) = this_results.GC_x1_to_x2;
    gcs_x2tox1(inst_inc_level) = this_results.GC_x2_to_x1;
    stoch_ints(inst_inc_level) = this_results.stoch_int;
    pred_infos(inst_inc_level) = this_results.pred_info;
    phis(inst_inc_level) = this_results.phi;
    
    
end

%% time domain 
plot(inst_ints,phis,'-ob',inst_ints,gcs_x1tox2,'-or',inst_ints,gcs_x2tox1,'or--',...
inst_ints,stoch_ints,'-og',inst_ints,pred_infos,'-ok');

legend({'phis','gcx1tox2','gcs_x2tox1','stochint','pred'})

savefig(['./' fl_desc '_td_vsinstintr.fig'])
saveas(gcf, ['./' fl_desc '_td_vsinstintr.svg'],'svg')

%%  reduced models investigations

phi_reduced_model_results = load('unidir_no_instincremental_inst_phi_reduced_model_results_min0max0.65nsteps10.mat');
phi_reduced_model_results = phi_reduced_model_results.incremental_inst_phi_reduced_model_results;

GC_x1_to_x2_model_reduced_model_results = load('unidir_no_instincremental_inst_GC_x1_to_x2_model_reduced_model_results_min0max0.65nsteps10.mat');
GC_x1_to_x2_model_reduced_model_results = GC_x1_to_x2_model_reduced_model_results.incremental_inst_GC_x1_to_x2_model_reduced_model_results;

GC_x2_to_x1_model_reduced_model_results = load('unidir_no_instincremental_inst_GC_x2_to_x1_model_reduced_model_results_min0max0.65nsteps10.mat');
GC_x2_to_x1_model_reduced_model_results = GC_x2_to_x1_model_reduced_model_results.incremental_inst_GC_x2_to_x1_model_reduced_model_results;

%%
phi_reduced_model_coord = [];
GC_x1_to_x2_reduced_model_coord = [];
GC_x2_to_x1_reduced_model_coord = [];

%%
for inst_inc_level = 1:length(incremental_inst_increase_results)
    
    % phi
    phi_reduced_model_coord(inst_inc_level).GC_x1_to_x2 = phi_reduced_model_results(inst_inc_level).results.GC_x1_to_x2;
    phi_reduced_model_coord(inst_inc_level).GC_x2_to_x1 = phi_reduced_model_results(inst_inc_level).results.GC_x2_to_x1;
    phi_reduced_model_coord(inst_inc_level).inst_int = phi_reduced_model_results(inst_inc_level).results.inst_int;
    
    % gc dir 1
    GC_x1_to_x2_reduced_model_coord(inst_inc_level).GC_x1_to_x2 = GC_x1_to_x2_model_reduced_model_results(inst_inc_level).results.GC_x1_to_x2; % this must be zero for this model
    GC_x1_to_x2_reduced_model_coord(inst_inc_level).GC_x2_to_x1 = GC_x1_to_x2_model_reduced_model_results(inst_inc_level).results.GC_x2_to_x1; 
    GC_x1_to_x2_reduced_model_coord(inst_inc_level).inst_int = GC_x1_to_x2_model_reduced_model_results(inst_inc_level).results.inst_int;
    
    % gc dir 2, this is zero in full model, so reduced model should be the same as full model
    GC_x2_to_x1_reduced_model_coord(inst_inc_level).GC_x1_to_x2 = GC_x2_to_x1_model_reduced_model_results(inst_inc_level).results.GC_x1_to_x2;
    GC_x2_to_x1_reduced_model_coord(inst_inc_level).GC_x2_to_x1 = GC_x2_to_x1_model_reduced_model_results(inst_inc_level).results.GC_x2_to_x1; % this is must be zero for this model
    GC_x2_to_x1_reduced_model_coord(inst_inc_level).inst_int = GC_x2_to_x1_model_reduced_model_results(inst_inc_level).results.inst_int;
    
end

%% start again with full model

clf
inst_inc_level = 1;
j = 1;
tmp = incremental_inst_increase_results(inst_inc_level).results;
full_results = [tmp.GC_x1_to_x2; tmp.GC_x2_to_x1;tmp.inst_int];
subplot(4,2,j)
bar(full_results')
ylim([-0.1 1.1])

inst_inc_level = length(incremental_inst_increase_results);
j = j +1;
tmp = incremental_inst_increase_results(inst_inc_level).results;
full_results = [tmp.GC_x1_to_x2; tmp.GC_x2_to_x1;tmp.inst_int];
subplot(4,2,j)
bar(full_results')
ylim([-0.1 1.1])

% phi reduced model properties
inst_inc_level = 1;
phi_no_inst_level_results = [phi_reduced_model_coord(inst_inc_level).GC_x1_to_x2 phi_reduced_model_coord(inst_inc_level).GC_x2_to_x1  ...
    phi_reduced_model_coord(inst_inc_level).inst_int];

j = j +1;
subplot(4,2,j)
bar(phi_no_inst_level_results')
ylim([-0.1 1.1])

inst_inc_level = length(incremental_inst_increase_results);
phi_no_inst_level_results = [phi_reduced_model_coord(inst_inc_level).GC_x1_to_x2 phi_reduced_model_coord(inst_inc_level).GC_x2_to_x1  ...
    phi_reduced_model_coord(inst_inc_level).inst_int];

j = j +1;
subplot(4,2,j)
bar(phi_no_inst_level_results')
ylim([-0.1 1.1])
% gc reduced model
inst_inc_level = 1;
phi_no_inst_level_results = [GC_x2_to_x1_reduced_model_coord(inst_inc_level).GC_x1_to_x2 GC_x2_to_x1_reduced_model_coord(inst_inc_level).GC_x2_to_x1  ...
    GC_x2_to_x1_reduced_model_coord(inst_inc_level).inst_int];

j = j +1;

subplot(4,2,j)
bar(phi_no_inst_level_results')
ylim([-0.1 1.1])

inst_inc_level = length(incremental_inst_increase_results);
phi_no_inst_level_results = [GC_x2_to_x1_reduced_model_coord(inst_inc_level).GC_x1_to_x2 GC_x2_to_x1_reduced_model_coord(inst_inc_level).GC_x2_to_x1  ...
    GC_x2_to_x1_reduced_model_coord(inst_inc_level).inst_int];

j = j +1;

subplot(4,2,j)
bar(phi_no_inst_level_results')
ylim([-0.1 1.1])
% delta

inst_inc_level = 1;
phi_no_inst_level_results = [phi_reduced_model_coord(inst_inc_level).GC_x1_to_x2 phi_reduced_model_coord(inst_inc_level).GC_x2_to_x1  ...
    phi_reduced_model_coord(inst_inc_level).inst_int] -[GC_x2_to_x1_reduced_model_coord(inst_inc_level).GC_x1_to_x2 GC_x2_to_x1_reduced_model_coord(inst_inc_level).GC_x2_to_x1  ...
    GC_x2_to_x1_reduced_model_coord(inst_inc_level).inst_int];

j = j +1;

subplot(4,2,j)
bar(phi_no_inst_level_results')
ylim([-0.04 0.04])
set(gca,'xticklabel',{'GC dir 1', 'GC dir 2', 'inst int'})


inst_inc_level = length(incremental_inst_increase_results);
phi_no_inst_level_results = [phi_reduced_model_coord(inst_inc_level).GC_x1_to_x2 phi_reduced_model_coord(inst_inc_level).GC_x2_to_x1  ...
    phi_reduced_model_coord(inst_inc_level).inst_int] - [GC_x2_to_x1_reduced_model_coord(inst_inc_level).GC_x1_to_x2 GC_x2_to_x1_reduced_model_coord(inst_inc_level).GC_x2_to_x1  ...
    GC_x2_to_x1_reduced_model_coord(inst_inc_level).inst_int];


j = j +1;

subplot(4,2,j)
bar(phi_no_inst_level_results')
ylim([-0.04 0.04])
set(gca,'xticklabel',{'GC dir 1', 'GC dir 2', 'inst int'})
title(['Delta is ' num2str(phi_no_inst_level_results)])

savefig(['./' fl_desc '_td_comparing_reduced_model_props_bar.fig'])
saveas(gcf, ['./' fl_desc '_td_comparing_reduced_model_props_bar.svg'],'svg')

%%
clf
hold on
% plot([GC_x2_to_x1_reduced_model_coord(:).GC_x1_to_x2])
plot(inst_ints,[GC_x2_to_x1_reduced_model_coord(:).GC_x2_to_x1],'r')
plot(inst_ints,[GC_x2_to_x1_reduced_model_coord(:).inst_int], '--r')


% clf
% subplot(1,3,2)
% hold on
% plot([GC_x1_to_x2_reduced_model_coord(:).GC_x1_to_x2])
% plot([GC_x1_to_x2_reduced_model_coord(:).GC_x2_to_x1])
% plot(inst_ints,[GC_x1_to_x2_reduced_model_coord(:).inst_int])

% clf
% hold on
% subplot(1,3,3)
% plot(phi_reduced_model_coord(:).GC_x1_to_x2)
% plot(phi_reduced_model_coord(:).GC_x2_to_x1)
plot(inst_ints,[phi_reduced_model_coord(:).inst_int],'b')
legend({'GC_x2_to_x1 red. model. GC x2 to x1', 'GC_x2_to_x1 red. model.', 'Phi red. model. inst int'})

savefig(['./' fl_desc '_td_comparing_reduced_model_props_lplot.fig'])
saveas(gcf, ['./' fl_desc '_td_comparing_reduced_model_props_lplot.svg'],'svg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Frequency domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
quant_names = {'sdecomp_phi','sdecomp_GC_x1_to_x2','sdecomp_stoch_int','sdecomp_pred_info','sdecomp_inst_int'};

for inst_inc_level = 1:length(incremental_inst_increase_results)
    
    this_results = incremental_inst_increase_results(inst_inc_level).results;
    freqs = this_results.freqs;
    for quant_indx = 1:length(quant_names)
        quant_indx
        subplot(length(quant_names),1,quant_indx)  
        title(quant_names{quant_indx})
        hold on
        vals = squeeze(real(this_results.(quant_names{quant_indx})));
        if size(vals,2) == 1 && size(vals,1) == 1
            vals = ones(1,length(freqs)) * vals;
        end
        hl = plot(freqs,vals,'k','color', [1,1,1]/(inst_inc_level),'linewidth', 3/(length(incremental_inst_increase_results)+1-inst_inc_level));
%      hl_GC = plot(freqs,incremental_inst_increase_results(level_indx).spct_GC_x1_to_x2,'k-.');
     
%          if level_indx == 1 || level_indx == length(incremental_inst_increase_results)         
%             set(hl,'linewidth',3)
%          end
    end
     
end

%increase size, looks better
set(gcf,'position',[1 82 1280 1263])

savefig(['./' fl_desc '_spct_grad_all.fig'])
saveas(gcf, ['./' fl_desc '_spct_grad_all.svg'],'svg')


%% now plot only highest level
clf
j=1;
subplot(2,4,j);

quant_names = {'sdecomp_phi','sdecomp_GC_x1_to_x2'};
colors = 'br';

inst_inc_level = length(incremental_inst_increase_results);
    
this_results = incremental_inst_increase_results(inst_inc_level).results;
freqs = this_results.freqs;
for quant_indx = 1:length(quant_names)
    quant_indx
%     subplot(length(quant_names),1,quant_indx)
%     title(quant_names{quant_indx})
    hold on
    vals = squeeze(real(this_results.(quant_names{quant_indx})));
    if size(vals,2) == 1 && size(vals,1) == 1
        vals = ones(1,length(freqs)) * vals;
    end
    hl = plot(freqs,vals,colors(quant_indx),'linewidth', 3/(length(incremental_inst_increase_results)+1-inst_inc_level));
    %      hl_GC = plot(freqs,incremental_inst_increase_results(level_indx).spct_GC_x1_to_x2,'k-.');
    
    %          if level_indx == 1 || level_indx == length(incremental_inst_increase_results)
    %             set(hl,'linewidth',3)
    %          end
end




%% plot diff
j = j + 4;
subplot(2,4,j);

this_results = incremental_inst_increase_results(end).results;

vals = squeeze(real(this_results.sdecomp_phi - this_results.sdecomp_GC_x1_to_x2));

if size(vals,2) == 1 && size(vals,1) == 1
   vals = ones(1,length(freqs)) * vals;
end
hl = plot(freqs,vals,'k','linewidth', 3/(length(incremental_inst_increase_results)+1-inst_inc_level));

%% now the componenets
freq_res = 100;
GC_x1_to_x2_model_H_r =  var2trfun(this_results.GC_x2_to_x1_model.A_r,freq_res); %arg, defined the other way in the scrpit
GC_x1_to_x2_model_H_r_11 = log( squeeze( GC_x1_to_x2_model_H_r(1,1,:) .* conj(GC_x1_to_x2_model_H_r(1,1,:))) );
GC_x1_to_x2_model_H_r_22 = log( squeeze( GC_x1_to_x2_model_H_r(2,2,:) .* conj(GC_x1_to_x2_model_H_r(2,2,:))) );
% 
% clf
% plot(GC_x1_to_x2_model_H_r_11)
% hold on
% plot(GC_x1_to_x2_model_H_r_22)

phi_model_H_r =  var2trfun(this_results.Phi_model.A_r,freq_res);
phi_model_H_r_11 = log( squeeze( phi_model_H_r(1,1,:) .* conj(phi_model_H_r(1,1,:))) );
phi_model_H_r_22 = log( squeeze( phi_model_H_r(2,2,:) .* conj(phi_model_H_r(2,2,:))) );

%%


ylims = [-1.4 1.1]
subplot(2,4,2);
plot(freqs,GC_x1_to_x2_model_H_r_11 ,'r', 'linewidth',3);
hold on 
plot(freqs,phi_model_H_r_11 ,'b', 'linewidth',3);
ylim(ylims)

subplot(2,4,3);
plot(freqs,GC_x1_to_x2_model_H_r_22,'r', 'linewidth',3);
hold on 
plot(freqs,phi_model_H_r_22 ,'b', 'linewidth',3);
ylim(ylims)

subplot(2,4,4);
plot(freqs,ones(1,length(freqs))*log(det(this_results.GC_x1_to_x2_model.SIG_r)) ,'r', 'linewidth',3);
hold on 
plot(freqs,ones(1,length(freqs))*log(det(this_results.Phi_model.SIG_r)) ,'b', 'linewidth',3);
ylim(ylims)

%% diffs

subplot(2,4,6);
plot(freqs,phi_model_H_r_11 - GC_x1_to_x2_model_H_r_11 ,'k', 'linewidth',3);

subplot(2,4,7);
plot(freqs,phi_model_H_r_22 - GC_x1_to_x2_model_H_r_22,'k', 'linewidth',3);


subplot(2,4,8);
plot(freqs,ones(1,length(freqs))*(log(det(this_results.Phi_model.SIG_r)) - log(det(this_results.GC_x1_to_x2_model.SIG_r))) ,'k', 'linewidth',3);


%%

%increase size, looks better
set(gcf,'position',[1 750 2560 595])

savefig(['./' fl_desc '_spct_phi_and_gc_only_eiglike_decomp.fig'])
saveas(gcf, ['./' fl_desc '_spct_phi_and_gc_only_eiglike_decomp.svg'],'svg')


%% examine the highest level of inst gc and phi


%% break down to structure of h

%% 
clf
freq_res = 100;
for level_indx = 1:length(GC_res.GC_incremental_inst_increase_results)
    hold on
    H_r =  var2trfun(GC_res.GC_incremental_inst_increase_results(level_indx).A_r,freq_res);
    hl = plot(freqs,squeeze(abs(H_r(1,2,:)).^2),'K')
    
     if level_indx == 1 || level_indx == length(incremental_inst_increase_results)         
        set(hl,'linewidth',3)
     end
%     plot(freqs,squeeze(abs(H_r(1,2,:)))/max(squeeze(abs(H_r(1,2,:)))))
%       plot(freqs,squeeze(real(incremental_inst_increase_results(level_indx).spct_phi)),'b','linewidth',1)
%      plot(freqs,incremental_inst_increase_results(level_indx).spct_GC_x1_to_x2,'k.')
    
end

title('abs(Hr(1,2)) GX2toX1')

savefig(['./' fl_desc '_spct_grad_Gcinoppdir.fig'])
saveas(gcf, ['./' fl_desc '_spct_grad_Gcinoppdir.svg'],'svg')

%% spectral gc 


%firs confirm that you obtain the same value as
%gc
level_indx_chk = length(GC_res.GC_incremental_inst_increase_results);
H_r =  var2trfun(GC_res.GC_incremental_inst_increase_results(level_indx_chk).A_r,freq_res);
Sig_r = GC_res.GC_incremental_inst_increase_results(level_indx_chk).SIG_r;

logdetSr = zeros(1,size(H_r,3));
logdetSr1 = logdetSr;

Hr_11_sq = abs(H_r(1,1,:)).^2;
Hr_22_sq = abs(H_r(2,2,:)).^2;
tmp1 = Hr_11_sq .* Hr_22_sq * det(Sig_r);

for freqind = 1:size(H_r,3)
        
    tmp =det( squeeze(H_r(:,:,freqind)) * Sig_r * squeeze( H_r(:,:,freqind))' );
  
    det_S_f = det( squeeze(GC_res.GC_incremental_inst_increase_results(level_indx_chk).S_f(:,:,freqind)) );
    det_S_r = det( squeeze(GC_res.GC_incremental_inst_increase_results(level_indx_chk).S_r(:,:,freqind)));
    
    logdetSr(freqind) = real(log ( tmp / det_S_f) );
    logdetSr1(freqind) = real(log ( det_S_r / det_S_f ));
    logdetSr2(freqind) = real(log ( tmp1(freqind) / det_S_f ));

end

%
clf
hold on
plot(real(GC_res.GC_incremental_inst_increase_results(level_indx_chk).spct_phi))
plot(real(incremental_inst_increase_results(level_indx_chk).spct_GC_x1_to_x2),'--')
% plot(real(incremental_inst_increase_results(level_indx_chk).spct_phi),'r--')
plot(logdetSr2,':','linewidth',3)
plot(logdetSr,'.-')
plot(logdetSr1,'--')

%% phi
H_r =  var2trfun(incremental_inst_increase_results(level_indx_chk).phi_A_r,freq_res);
Sig_r = incremental_inst_increase_results(level_indx_chk).phi_SIG_r;

logdetSr = zeros(1,size(H_r,3));
logdetSr1 = logdetSr;

Hr_11_sq = abs(H_r(1,1,:)).^2;
Hr_22_sq = abs(H_r(2,2,:)).^2;
tmp1 = Hr_11_sq .* Hr_22_sq * det(Sig_r);

log_det_S_fs = []

for freqind = 1:size(H_r,3)
        
    tmp =det( squeeze(H_r(:,:,freqind)) * Sig_r * squeeze( H_r(:,:,freqind))' );
  
    det_S_f = det( squeeze(GC_res.GC_incremental_inst_increase_results(level_indx_chk).S_f(:,:,freqind)) );
    det_S_r = det( squeeze(incremental_inst_increase_results(level_indx_chk).phi_S_r(:,:,freqind)));
    
    det_S_fs(freqind) = det_S_f;
    
    logdetSr(freqind) = real(log ( tmp / det_S_f) );
    logdetSr1(freqind) = real(log ( det_S_r / det_S_f ));
    logdetSr2(freqind) = real(log ( tmp1(freqind) / det_S_f ));

end

%%
clf
hold on
diff = real(incremental_inst_increase_results(level_indx_chk).spct_phi)-real(incremental_inst_increase_results(level_indx_chk).spct_GC_x1_to_x2);
plot(freqs,diff)
savefig(['./' fl_desc '_phiMinusGC.fig'])
saveas(gcf, ['./' fl_desc '_phiMinusGC.svg'],'svg')

% plot(logdetSr2,':','linewidth',3)
% plot(logdetSr,'.-')
% plot(logdetSr1,'--')

%% diff

%% now, decompose
clf
for level_indx = length(GC_res.GC_incremental_inst_increase_results)
    
    phi_H_r =  var2trfun(incremental_inst_increase_results(level_indx).phi_A_r,freq_res);
    gc_H_r =  var2trfun(GC_res.GC_incremental_inst_increase_results(level_indx).A_r,freq_res);
    
    phi_Sig_r = incremental_inst_increase_results(level_indx).phi_SIG_r;
    gc_Sig_r = GC_res.GC_incremental_inst_increase_results(level_indx).SIG_r; 
    
    
    phi_Hr_11_sq = squeeze(abs(phi_H_r(1,1,:)).^2);
    phi_Hr_22_sq = squeeze(abs(phi_H_r(2,2,:)).^2);
    
%     phi_Hr_11_sq = squeeze(phi_H_r(1,1,:) .* conj(phi_H_r(1,1,:)));
%     phi_Hr_22sq = squeeze(phi_H_r(2,2,:) .* conj(phi_H_r(2,2,:)));
% %     phi_Hr_22_sq = squeeze(phi_H_r(2,2,:)).^2);
    
    gc_Hr_11_sq = squeeze(abs(gc_H_r(1,1,:)).^2);
    gc_Hr_22_sq = squeeze(abs(gc_H_r(2,2,:)).^2);
        
    hold on
%     hl = plot(freqs,squeeze(gc_Hr_11_sq),'r',...
%         freqs,squeeze(gc_Hr_22_sq),'r--',...
%         freqs,squeeze(phi_Hr_11_sq),'b',...
%         freqs,squeeze(phi_Hr_22_sq),'b--',...
%         freqs,log_det_S_f);
    
   hl = plot(freqs,log(squeeze(gc_Hr_11_sq)),'r',...
        freqs,log(squeeze(gc_Hr_22_sq)),'r--',...
        freqs,log(squeeze(phi_Hr_11_sq)),'b',...
        freqs,log(squeeze(phi_Hr_22_sq)),'b--',...
        freqs,ones(size(freqs)) * log(det(phi_Sig_r)),'b:',...
        freqs,ones(size(freqs)) * log(det(gc_Sig_r)),'r:'); % ,...freqs,log(det_S_fs))
    
     if level_indx == 1 || level_indx == length(incremental_inst_increase_results)         
        set(hl,'linewidth',3)
     end
    
end

legend(hl,{'GCHr11','GCHr22','phiHr11','phiHr22'});
title('H_R^2 for phi and gc')
 
%
savefig(['./' fl_desc '_spct_grad_GcEigendecomp.fig'])
saveas(gcf, ['./' fl_desc '_spct_grad_GcEigendecomp.svg'],'svg')

%%  chk

clf
plot(diff, 'b')
hold on
plot(-log(squeeze(gc_Hr_11_sq))-log(squeeze(gc_Hr_22_sq))...
    +log(squeeze(phi_Hr_11_sq))+log(squeeze(phi_Hr_22_sq)) + log(det(phi_Sig_r)) - log(det(gc_Sig_r)),'r--');

%%
clf
   hl = plot(freqs,log(squeeze(phi_Hr_11_sq)) - log(squeeze(gc_Hr_11_sq)),'r',...
        freqs, log(squeeze(phi_Hr_22_sq)) - log(squeeze(gc_Hr_22_sq)),'b--',...
        freqs,ones(size(freqs)) * (log(det(phi_Sig_r))-log(det(gc_Sig_r))),'b:',...
        freqs,diff,'k'); % ,...freqs,log(det_S_fs))

%%
clf
 hl = plot(freqs,log(squeeze(gc_Hr_11_sq)) + log(squeeze(gc_Hr_22_sq)),'r--',...
        freqs,log(squeeze(phi_Hr_11_sq)) +log(squeeze(phi_Hr_22_sq)),'b--');
    
clf
hl = plot(freqs,log(squeeze(phi_Hr_11_sq)) +log(squeeze(phi_Hr_22_sq)) - ...
    (log(squeeze(gc_Hr_11_sq)) + log(squeeze(gc_Hr_22_sq)) ) ,'k--');

hold on
plot(freqs,real(incremental_inst_increase_results(end).spct_phi) - ...
    real(incremental_inst_increase_results(end).spct_GC_x1_to_x2 ))

clf
hold on
hl = plot(freqs,log(squeeze(gc_Hr_11_sq)) - log(squeeze(phi_Hr_11_sq)),'b--');
hl = plot(freqs,log(squeeze(gc_Hr_22_sq)) - log(squeeze(phi_Hr_22_sq)),'b--');
%%
% clf
% plot3(freqs,phi_Hr_11_sq,phi_Hr_22_sq)
% hold on
% plot3(freqs,gc_Hr_11_sq,gc_Hr_22_sq)
% 
% %% Add decomp for stoch int
% clf
% plot(abs(S_gchack'))
% % plot(abs(S_gchack')-log(abs(GC_res.det_S_f)))
% 
% %% what explains the diff between gc and phi?
% clf
% % plot(log(real(det_S_f)))
% 
% %% get tra fun for phi
% clf
% hold on
% % [H_phi sig_phi] = cpsd_to_var(S_r);
% numerator_phi = det(sig_phi) * squeeze( H_phi(1,1,:) .* H_phi(2,2,:) .* conj(H_phi(1,1,:) .* H_phi(2,2,:)));
% 
% %% gc
% split_mask_A=ones(N);
% split_mask_A(2,1)=0;
% % correlation in the noise covariance 
% split_mask_E=ones(N);
% 
% 
% % Both time and frequency domain quantities are estimated
% disp('Computing reduced model parameters')
% [S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(G,split_mask_A,split_mask_E,max_order,freq_res,iter_max,gamma,min_error);    
% 
%  
% %%
% % [H_gc sig_gc] = cpsd_to_var(S_r);
% H_gc = var2trfun(A_r,freq_res);
% numerator_gc = det(SIG_r) * squeeze( H_gc(1,1,:) .* H_gc(2,2,:) .* conj(H_gc(1,1,:) .* H_gc(2,2,:)));
% 
% %%
% clf
% hold on
% plot(freqs,log(numerator_gc ./ real(det_S_f) ) ,'r.')
% plot(freqs,log(numerator_phi ./real(det_S_f) ) ,'b.')
% plot(freqs,incremental_inst_increase_results(level_indx).spct_GC_x1_to_x2,'k.')
% %       plot(freqs,squeeze(real(incremental_inst_increase_results(level_indx).spct_phi)),'b','linewidth',5)
% 
% diff = log(numerator_gc ./ real(det_S_f) ) - log(numerator_phi ./real(det_S_f) );
% %%
% clf
% hold on
% % plot(freqs,log(det(SIG_r)*squeeze(H_gc(1,1,:).*conj(H_gc(1,1,:))) ./ real(det_S_f) ) ,'r')
% % plot(freqs,log(det(SIG_r)*squeeze(H_phi(1,1,:).*conj(H_phi(1,1,:))) ./real(det_S_f) ) ,'b')
% 
% % plot(freqs,log(squeeze(H_gc(1,1,:).*conj(H_gc(1,1,:)))  ) ,'r')
% plot(freqs,log(squeeze(H_phi(1,1,:).*conj(H_phi(1,1,:))) ) ,'b')
% 
% 
% % plot(freqs,incremental_inst_increase_results(level_indx).spct_GC_x1_to_x2,'k.')
%  plot(freqs,squeeze(real(incremental_inst_increase_results(level_indx).spct_phi)),'b','linewidth',1)
% 
% % plot(freqs,log(squeeze(H_gc(2,2,:).^2) ./ real(det_S_f) ) ,'r--')
% % plot(freqs,log(squeeze(H_phi(2,2,:) .* conj(H_phi(2,2,:))) ./real(det_S_f) ) + log(squeeze(H_phi(1,1,:).^2))+log(det(sig_phi)),'r--')
% % % plot(freqs,log(squeeze(H_phi(2,2,:) .* conj(H_phi(2,2,:))) ./real(det_S_f) ) + log(squeeze(H_phi(1,1,:).^2))+log(det(sig_phi)),'r--')
% % plot(freqs,log(squeeze(H_phi(2,2,:) .* conj(H_phi(1,1,:))) ./real(det_S_f) ) + log(det(sig_phi)),'r--')
% % plot(freqs,log(squeeze(H_phi(2,2,:) .* conj(H_phi(2,2,:))) ./real(det_S_f) ) + log(det(sig_phi)),'r--')
% 
% % plot(freqs,log(squeeze(H_gc(2,2,:).*conj(H_gc(2,2,:)))) ,'r--')
% plot(freqs,log(squeeze(H_phi(2,2,:).*conj(H_phi(2,2,:)))) ,'b--')
% 
% 
% % diff_h1 = log(squeeze(H_gc(1,1,:).^2)./ real(det_S_f)) - log(squeeze(H_phi(1,1,:).^2) ./ real(det_S_f));
% % diff_h2 = log(squeeze(H_gc(2,2,:).^2)./ real(det_S_f)) - log(squeeze(H_phi(2,2,:).^2) ./ real(det_S_f));
% 
% %
% diff_h1 = log( squeeze( H_gc(1,1,:) .* H_gc(2,2,:) .* conj(H_gc(1,1,:) .* H_gc(2,2,:)))) - log(squeeze(H_phi(1,1,:).^2) );
% diff_h2 = log( squeeze( H_gc(1,1,:) .* H_gc(2,2,:) .* conj(H_gc(1,1,:) .* H_gc(2,2,:)))) - log(squeeze(H_phi(2,2,:).^2) );
% 
% %%
% clf
% hold on
% plot(freqs,diff_h1,'r',freqs,diff, 'b')
% plot(freqs,diff_h2,'r--',freqs,diff, 'b--')
% plot(freqs,squeeze(real(incremental_inst_increase_results(level_indx).spct_phi)),'b','linewidth',5)
% 
% 
% %%
% clf
% hold on
% tmp = squeeze(H_gc(1,2,:).* conj(H_gc(1,2,:)));
% % tmp = tmp - min(tmp);
% % tmp = tmp / max(tmp) * max();
% plot(freqs,(tmp))
% plot(freqs,squeeze(real(incremental_inst_increase_results(level_indx).spct_phi)),'b','linewidth',5)

 




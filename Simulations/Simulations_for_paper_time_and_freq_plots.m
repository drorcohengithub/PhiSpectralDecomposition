clear all
close all

%% load sys
systems = Example_systems();

%% choose system you want to investigate
system_indx = 1;
system_nm = systems(system_indx).nm;

%% load

GC_data = load(['./' system_nm '_GC.mat']); 
GC_x2_to_x1 = GC_data.results.GC_x2_to_x1;
sdecomp_GC_x2_to_x1 = GC_data.results.sdecomp_GC_x2_to_x1;

GC_x1_to_x2 = GC_data.results.GC_x1_to_x2;
sdecomp_GC_x1_to_x2 = GC_data.results.sdecomp_GC_x1_to_x2;


stoch_int_data = load(['./' system_nm '_stoch_int.mat']);
stoch_int = stoch_int_data.results.stoch_int;
sdecomp_stoch_int = stoch_int_data.results.sdecomp_stoch_int;


predinfo_data = load(['./' system_nm '_pred_info.mat']);
pred_info = predinfo_data.results.pred_info;
sdecomp_pred_info = predinfo_data.results.sdecomp_pred_info;


inst_int_data = load(['./' system_nm '_inst_int.mat']);   
inst_int = inst_int_data.results.inst_int;
sdecomp_inst_int = inst_int_data.results.sdecomp_inst_int;


Phig_data = load(['./' system_nm '_phi.mat']);
% phig.results.
Phi_G = Phig_data.results.phi;
sdecomp_Phi_G = Phig_data.results.sdecomp_phi;

% these are the same for all measures
freqs = Phig_data.results.freqs;

%% plot time domain
time_data = [stoch_int pred_info Phi_G GC_x1_to_x2 GC_x2_to_x1 inst_int];   
clrs = {'g','k','b','b--','r','y'};
max_iter = length(time_data);
offset = 0.005;

%%

clf
subplot(2,1,1);
set(gca,'fontsize',20)
hold on   
for i = 1:max_iter
    h=bar(i,time_data(i));
    text(i,time_data(i)+offset,num2str(time_data(i)),'fontsize',30)  
    set(h,'FaceColor',clrs{i}(1));
end

set(gca,'xtick',1:max_iter)
% set(gca,'xticklabel',{'multi info','PhiG','GC1to2','GC2to1','sum of GCs'})
set(gca,'xticklabel',{'stochint','predinfo','phig','GC1to2','GC2to1','instinfo'})
maxy= max(time_data) * 1.05;
ylim([0 maxy]);

%% plot spectral decomposition
%% with pred info


subplot(2,1,2);

spect_data = real([sdecomp_stoch_int sdecomp_pred_info sdecomp_Phi_G sdecomp_GC_x1_to_x2 sdecomp_GC_x2_to_x1 sdecomp_inst_int])';   

hold on
for i = 1:size(spect_data,1)
    h=plot(freqs,spect_data(i,:),clrs{i});
   
end

%%
legend({'stochint','predinfo','phig','GC1to2','GC2to1','instinfo'}, 'fontsize',20)

saveas(gcf, ['./' system_nm '_rawlimits_withpredinfo'],'svg')

cf = gcf();
set(cf,'Position', [1281 87 1280 1258])

%%
ylim([0 max(real(sdecomp_stoch_int))*1.05])

saveas(gcf, ['./' system_nm '_withpredinfo'],'svg')

%% without pred info
clf
subplot(2,1,2);

spect_data = real([sdecomp_stoch_int sdecomp_Phi_G sdecomp_GC_x1_to_x2 sdecomp_GC_x2_to_x1 sdecomp_inst_int])';   

hold on
for i = 1:size(spect_data,1)
    h=plot(freqs,spect_data(i,:),clrs{i});
   
end

%%
legend({'stochint','phig','GC1to2','GC2to1','instinfo'}, 'fontsize',20)

saveas(gcf, ['./' system_nm '_rawlimits_withoutpredinf'],'svg')

cf = gcf();
set(cf,'Position', [1281 87 1280 1258])

%%
ylim([0 max(real(sdecomp_stoch_int))*1.05])

saveas(gcf, ['./' system_nm '_withoutpredinfo'],'svg')

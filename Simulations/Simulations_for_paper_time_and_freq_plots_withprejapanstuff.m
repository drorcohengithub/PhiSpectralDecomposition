clear all
close all

%% load systems
systems = Example_systems();
%% choose system you want
system_indx = 5;


%% plot time domain
time_data = [multi_info Phi_G GC_x1_to_x2 GC_x2_to_x1 GC_x1_to_x2+GC_x2_to_x1];   

max_iter = length(time_data);
offset = 0.005;

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
set(gca,'xticklabel',[])

%% save
svFldr = ['/Users/dror/Google Drive/Sasai/FIGURES/Simulations/' systems(system_indx).nm];
saveFigAs(svFldr,0,'halfScreen')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examine eigenvalue decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Breakdown the determinants of the reduced models to their eigenvalues
    
% full system
eigvals_S = [];

% Phi G
eigvals_S_r = [];

% GC
eigvals_S_r_GC = [];

% multi info (eigenvalues are just the power spcts)
eigvals_multiinfo = [];

S_multiinfo = S;

S_multiinfo(1,2,:) = 0;
S_multiinfo(2,1,:) = 0;

%% also look at GC
% to do this by enforce off-diagonal to zero for GC
split_mask= ones(N);
split_mask(2,1)=0;
    
%this takes longer to converge
iter_max = 10000;
gamma = 0.01;
min_error = 1e-6;
[S_r_GC,det_S_r_GC,trace_S_r_GC,prod_diag_S_r_GC,A_r_GC,SIG_r_GC,masked_Delta_GC] =...
    get_reduced_S_from_autoCov(G,split_mask,max_order,freq_res,iter_max,gamma,min_error);
%% get the eigen valus

for freqIndx=1:size(S,3)
    
    [ignr tmp] = eig(squeeze(S(:,:,freqIndx)));
    eigvals_S(freqIndx,:)=diag(tmp);
    
    [ignr tmp] = eig(squeeze(S_r_GC(:,:,freqIndx)));
    eigvals_S_r_GC(freqIndx,:)=diag(tmp); 
    
    [ignr tmp] = eig(squeeze(S_r(:,:,freqIndx)));
    eigvals_S_r(freqIndx,:)=diag(tmp);
    
    [ignr tmp] = eig(squeeze(S_multiinfo(:,:,freqIndx)));
    eigvals_multiinfo(freqIndx,:)=diag(tmp);
    
    
end

%% order the eigen values
eigvals_S = log(real(sort(eigvals_S,2,'descend')));
eigvals_S_r_GC = log(real(sort(eigvals_S_r_GC,2,'descend')));
eigvals_S_r = log(real(sort(eigvals_S_r,2,'descend')));
eigvals_multiinfo = log(real(sort(eigvals_multiinfo,2,'descend')));

%% normalize by det_S
eigvals_S_r_GC_norm = eigvals_S_r_GC - eigvals_S;
eigvals_S_r_norm = eigvals_S_r - eigvals_S;
eigvals_multiinfo_norm = eigvals_multiinfo - eigvals_S;

%store here
eig_cell={};
eig_cell{1,1} = eigvals_S_r_GC;
eig_cell{1,2} = eigvals_S_r_GC_norm;

eig_cell{2,1} = eigvals_S_r;
eig_cell{2,2} = eigvals_S_r_norm;

eig_cell{3,1} = eigvals_multiinfo;
eig_cell{3,2} = eigvals_multiinfo_norm;

reduced_S_nms = {'S_r_GC','S_r_Phi','S_r_multiinfo'};

%% plot properties for reduced sys 

clf
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultLineLinewidth',2)
% set(0,'DefaultAxesXlim','tight')
sys_y_lims = [];
sys_y_lims(3,1,:) = [-4 4];
sys_y_lims(3,2,:) = [-0.8 0.8];
sys_y_lims(4,1,:) = [-5 5];
sys_y_lims(4,2,:) = [-1.1 1.1 ];
sys_y_lims(5,1,:) = [-5 5];
sys_y_lims(5,2,:) = [-1.1 1.1 ];

j=1;
ax(j) = subplot(2,3,j); j = j+1;
hold on
plot(freqs,eigvals_S(:,1),...
    freqs,eigvals_S_r_GC(:,1),'r',...
    freqs,eigvals_S_r(:,1),'b',...
    freqs,eigvals_multiinfo(:,1),'g');

ylim(sys_y_lims(system_indx,1,:))
title('log(Eig 1)')

ax(j) = subplot(2,3,j); j = j+1;
plot(freqs,eigvals_S(:,2),...
    freqs,eigvals_S_r_GC(:,2),'r',...
    freqs,eigvals_S_r(:,2),'b',...
    freqs,eigvals_multiinfo(:,2),'g');

ylim(sys_y_lims(system_indx,1,:))
title('log(Eig 2)')

ax(j) = subplot(2,3,j); j = j+1;
plot(freqs,sum(eigvals_S,2),...
    freqs,sum(eigvals_S_r_GC,2),'r',...
    freqs,sum(eigvals_S_r,2),'b',...
    freqs,sum(eigvals_multiinfo,2),'g');

ylim(sys_y_lims(system_indx,1,:))
title('log(Eig 1) + log(Eig 1)')

% diff
ax(j) = subplot(2,3,j); j = j+1;
hold on
plot(freqs,eigvals_S_r_GC_norm(:,1),'r',...
    freqs,eigvals_S_r_norm(:,1),'b',...
    freqs,eigvals_multiinfo_norm(:,1),'g');

title('\Deltalog(Eig 1)')
xlabel('Frequency');
ylim(sys_y_lims(system_indx,2,:))


ax(j) = subplot(2,3,j); j = j+1;
plot(freqs,eigvals_S_r_GC_norm(:,2),'r',...
    freqs,eigvals_S_r_norm(:,2),'b',...
    freqs,eigvals_multiinfo_norm(:,2),'g');

title('\Deltalog(Eig 2)')
xlabel('Frequency');
ylim(sys_y_lims(system_indx,2,:))


ax(j) = subplot(2,3,j); j = j+1;
plot(freqs,sum(eigvals_S_r_GC_norm,2),'r',...
    freqs,sum(eigvals_S_r_norm,2),'b',...
    freqs,sum(eigvals_multiinfo_norm,2),'g');

title('\Deltalog(Eig 1)+\Deltalog(Eig 2)')
xlabel('Frequency');
ylim(sys_y_lims(system_indx,2,:))
  
% suptitle(reduced_S_nms{reduced_S_indx});

set(ax,'xlim',[0 500],'xticklabel',[])

% ,freqs,eigvals_S_norm(:,2))
% plot(freqs,log(eigvals_S_r(:,1)),freqs,log(eigvals_S_r(:,2)),freqs,sum(log(eigvals_S_r),2),freqs,log(det_S'))

%     set(ax,'Nextplot','add')
%     hold(ax,'on')
% save
svFldr = ['/Users/dror/Google Drive/Sasai/FIGURES/Simulations/' systems(system_indx).nm '_Eigdecomp']; % reduced_S_nms{reduced_S_indx}
saveFigAs(svFldr,0,'halfScreen')

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The code below verfies some analytical relationships presented in the paper. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    %% Check: Relationshup between spectral decomp of multi info and coherence
    C=squeeze((abs(S(1,2,:)).^2)./(S(1,1,:).*S(2,2,:)));
    %mutual info rate
    sdecomp_multi_info_2 = -log(1-C);
    clf
    plot(freqs,sdecomp_multi_info-sdecomp_multi_info_2); ylim(1e-10*[-1 1]);%  note the scale

    %% Check: Get GC through through reducde model framework

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
end


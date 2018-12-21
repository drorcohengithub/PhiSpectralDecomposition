clear all
close all

%% Multi information
% The current code does not support restrictions on the noise covariance
% matrix (will add in next version). However we can obtain it in a a number
% of equivalent ways. The first options is to fit the models seperately
univar_autocov1=G(1,1,:);
[univar_A1 univar_SIG1] = autocov_to_var(univar_autocov1);
univar_autocov2=G(2,2,:);
[univar_A2 univar_SIG2] = autocov_to_var(univar_autocov2);
multi_info = log ( univar_SIG1*univar_SIG2  / det (SIG_f) );

% The second option is to get this from the existing GC estimates
log_E22=GC_x1_to_x2+log(SIG_f(2,2));
log_E11=GC_x2_to_x1+log(SIG_f(1,1));
multi_info_2 = log_E22+log_E11 - log( det( SIG_f ) );

% the difference between the two is very small
disp(multi_info-multi_info_2)

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
% svFldr = ['/Users/dror/Google Drive/Sasai/FIGURES/Simulations/' systems(system_indx).nm];
% saveFigAs(svFldr,0,'halfScreen')





clear all
close all

%% load systems
systems = Example_systems();
%% choose system you want
systemIndx = 3;
A = systems(systemIndx).A;
SIG_f = systems(systemIndx).SIG_f; 

%% Get the autocov of the full system
[G,info] = var_to_autocov(A,SIG_f);
% Check that everything is ok
var_info(info,true);

%% input for analysis
% arbitrarily set the freq resolution and samp rate
freq_res  = 100;
samp_rate = 1000;
% obtan the frequency bins
freqs = sfreqs(freq_res,samp_rate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the time domain GC 
[t_GC]= autocov_to_pwcgc(G);
GC_x1_to_x2 = t_GC(2,1);
GC_x2_to_x1 = t_GC(1,2);

%% Phi G calculation
%Covriance and autcov function for phi calc
Cov_X  = G(:,:,1);
Cov_XY = G(:,:,2:end);
% The maximum lag we will calculate the reduced model up to.
% Note that for good estimate the maximum lag of the reduced model will be
% much larger than that of the full model. A safe but potentially
% over-generous strategy is to use the max lag of the autocov function
max_order=size(Cov_XY,3);

%% Chose form of reduced model. 
% The split_mask determines which connections to cut and which to keep. By
% cutting both lagged connections we will obtain Phi_G
N = size(G,1);
split_mask= eye(N);

%% Get the parameters of the reduced model (denoted with subscript r)
% Both time and frequency domain quantities are estimated
[S_r,det_S_r,trace_S_r,prod_diag_S_r,A_r,SIG_r,masked_Delta] = get_reduced_S_from_autoCov(G,split_mask,max_order,freq_res);

%% Finally, from the defintion of Phi G
Phi_G = log (det(SIG_r) / det (SIG_f));

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
if systemIndx<3 % for systems in which we are not interested in the sum of GC
    max_iter = length(time_data)-1;
    offset = 0;
else % for systems in which we'd like to see the sum of GC
    max_iter = length(time_data);
    offset = 0.005;
end

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

if systemIndx<3
    hold on
    hl=plot(freqs,squeeze(spct_GC(1,2,:)),'r',freqs,squeeze(spct_GC(2,1,:)),'r--');
    set(hl,'linewidth',4);
    % hl=plot(freqs,squeeze(spct_GC(1,2,:))+squeeze(spct_GC(2,1,:)),'b-.');
    % set(hl,'linewidth',2);
    plot(freqs,sdecomp_multi_info,'g');
    plot(freqs,sdecomp_Phi_G,'k:','linewidth',2)
    legend({'GC x_2 to x_1','GC x_1 to x_2','multi. info.','Phi G'},'fontsize',20,'linewidth',2)
    xlabel('Frequency')
    set(gca,'fontsize',20)
else
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
%     plot(freqs,squeeze(real(H(2,1,:))));

end


%% save
svFldr = ['/Users/dror/Google Drive/Sasai/FIGURES/Simulations/' systems(systemIndx).nm];
saveFigAs(svFldr,0,'halfScreen')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The code below verfies some analytical relationships presented in the paper. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check: Relationshup between spectral decomp of multi info and coherence
C=squeeze((abs(S(1,2,:)).^2)./(S(1,1,:).*S(2,2,:)));
%mutual info rate
sdecomp_multi_info_2 = -log(1-C);
plot(freqs,sdecomp_multi_info-sdecomp_multi_info_2); ylim(1e-10*[-1 1]);%  note the scale

%% Check: Get GC through through reducde model framework
% do this by enforce off-diagonal to zero for GC
split_mask= ones(N);
split_mask(2,1)=0;
    
%this takes longer to converge
iter_max = 10000;
[S_r_GC,det_S_r_GC,trace_S_r_GC,prod_diag_S_r_GC,A_r_GC,SIG_r_GC,masked_Delta_GC] = get_reduced_S_from_autoCov(G,split_mask,max_order,freq_res,iter_max);
    
%compare with MVGC
GC_x1_to_x2_2=log( det (SIG_r_GC) / det (SIG_f) );
%time domain
disp(GC_x1_to_x2-GC_x1_to_x2_2);

%freq domain
sdecomp_GC_x1_to_x2_2=log(det_S_r_GC./det_S');
plot(freqs,real(sdecomp_GC_x1_to_x2_2)-sdecomp_GC_x1_to_x2)

%check summation across freqs
GC_x1_to_x2_3 = mean([sdecomp_GC_x1_to_x2_2;sdecomp_GC_x1_to_x2_2(2:freq_res)]);
disp(real(GC_x1_to_x2_3 - GC_x1_to_x2_2))
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


% THE FOLLOWING DEFINTION OF time domain MULTI INFO IS NOT CORRECT
% multi_info_t1 = log ( prod ( diag (SIG_r) ) / det (SIG_f) );
% % this are not similar enough
% disp(multi_info_t1 - multi_info_t)
% 
% % this is because we need to evaluate the univar models individiaully!
% % get univar models
% univar_autocov1=G(1,1,:);
% [univar_A1 univar_SIG1] = autocov_to_var(univar_autocov1);
% univar_autocov2=G(2,2,:);
% [univar_A2 univar_SIG2] = autocov_to_var(univar_autocov2);
% multi_info_t2 = log ( univar_SIG2*univar_SIG2  / det (SIG_f) );
% % this is correct
% disp(multi_info_t2 - multi_info_t1)
% 
% % the power spectrums of the reduced and full model should be identical, 
% % due to the existence of the spectral decomposision
% % first check
% univar_S1 = var_to_cpsd(univar_A1,univar_SIG1,freq_res);
% clf
% plot(log(squeeze(real(univar_S1))));
% hold on
% plot( log(squeeze(S(1,1,:))),'b--');
% 
% univar_S2 = var_to_cpsd(univar_A2,univar_SIG2,freq_res);
% clf
% plot(log(squeeze(real(univar_S2))));
% hold on
% plot( log(squeeze(S(2,2,:))),'b--');
% 
% % now try 
% spctbased_log_univar_SIG2 = mean(log(squeeze(real(cat(3,univar_S2,univar_S2(:,:,2:freq_res)) ))))
% disp(spctbased_log_univar_SIG2 - log(univar_SIG2))
% 
% spctbased_log_univar_SIG1 = mean(log(squeeze(real(cat(3,univar_S1,univar_S1(:,:,2:freq_res)) ))))
% disp(spctbased_log_univar_SIG1 - log(univar_SIG1))

%% plot time domain
subplot(2,1,2);
set(gca,'fontsize',20)
hold on
if systemIndx<4
    tdata = [multi_info_t Phi_G GC_x1_to_x2 GC2to1];
    for i = 1:length(tdata)
        h=bar(i,tdata(i));
        text(i,tdata(i)+0.005,num2str(tdata(i)),'fontsize',30)
        if i == 1  
            set(h,'FaceColor','g');
        elseif i == 2
            set(h,'FaceColor','k');
        elseif i == 3 || i == 4
            set(h,'FaceColor','r');
        end
    end

    set(gca,'xtick',1:4)
    set(gca,'xticklabel',{'multi info','PhiG','GC1to2','GC2to1'})
else
    tdata = [multi_info_t Phi_G GC_x1_to_x2 GC2to1 GC_x1_to_x2+GC2to1];    
    for i = 1:length(tdata)
        h=bar(i,tdata(i));
        text(i-0.5,tdata(i)+0.005,num2str(tdata(i)),'fontsize',30)
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

    set(gca,'xtick',1:5)
    set(gca,'xticklabel',{'multi info','PhiG','GC1to2','GC2to1','sum of GCs'})
end
maxy= max(ylim);
ylim([0 maxy])



% more details stuff
if systemIndx>=3
    
    % enforce off-diagonal to zero for GC
    split_mask= ones(N);
    split_mask(2,1)=0;
    
    %this takes longer to converge
    iter_max = 10000;
    [S_r_GC,det_Sr_GC,trace_Sr_GC,prod_diag_Sr_GC,A_r_GC,SIG_r_GC,masked_Delta_GC] = get_reduced_S_from_autoCov(G,split_mask,max_order,freq_res,iter_max);
    
    %Compare var estimation with mvgc    
    univar_autocov2=G(2,2,:);
    [univar_A2 univar_SIG2] = autocov_to_var(univar_autocov2);
    % this is very similar
    disp(SIG_r_GC(2,2)-univar_SIG2)
    
    % compare time domain estimations of gc
    Phi_based_tm_domain_GC = log(det(SIG_r_GC)/det(SIG_f));
    % this is very good
    disp(Phi_based_tm_domain_GC - t_GC(2,1))
    
    % check freq domain
    Phi_GC_s  = log ( real(det_Sr_GC./det_S') );   
    clf
    hl=plot(freqs,squeeze(spct_GC(1,2,:)),'r',freqs,squeeze(spct_GC(2,1,:)),'r--');
    hold on
    hl=plot(freqs,squeeze(Phi_GC_s),'k-.');
    % very similar
    
    %check freq/time domain correspondance
    % phi based
    spct_based_phi_based_tGC = mean([Phi_GC_s;Phi_GC_s(2:freq_res)]);
    disp(spct_based_phi_based_tGC - Phi_based_tm_domain_GC);
    
    %mvgc
    mvgc_spct_based_tGC = mean(cat(3,spct_GC,spct_GC(:,:,2:freq_res)),3);
    disp(mvgc_spct_based_tGC(2,1) - t_GC(2,1));
    
    
    %% So all the time domain inequalities as well as the time/freq domain equalities hold.
    
    % The next thing is to try and explain the relationship between the measueres  
    % For example, what does it mean that mutual info rate is sometimes
    % smaller and sometime greater than Phi G?
    
    % get other GC 
    split_mask= ones(N);
    split_mask(1,2)=0;
    
    %this takes longer to converge
    iter_max = 10000;
    [S_r_GC1,det_Sr_GC1,trace_Sr_GC1,prod_diag_Sr_GC1,A_r_GC1,SIG_r_GC1,masked_Delta_GC1] = get_reduced_S_from_autoCov(G,split_mask,max_order,freq_res,iter_max);
    

    % lets breakdown the determinants of the reduced models to their eigen
    % values     
    S_eva = []; 
    Sr_eva = [];
    S_GC_eva = [];
    S_GC1_eva = [];

    S_eve = [];
    Sr_eve = [];
    S_GC_eve = [];
    S_GC1_eve = [];

    for freqIndx=1:size(S,3)
        
        [S_eve(:,:,freqIndx) tmp] = eig(squeeze(S(:,:,freqIndx)));
        [Sr_eve(:,:,freqIndx) tmp2] = eig(squeeze(S_r(:,:,freqIndx)));
        
        [S_GC_eve(:,:,freqIndx) tmp3] = eig(squeeze(S_r_GC(:,:,freqIndx)));
        [S_GC1_eve(:,:,freqIndx) tmp4] = eig(squeeze(S_r_GC1(:,:,freqIndx)));


        S_eva(1,freqIndx)=max(tmp(1,1),tmp(2,2));
        S_eva(2,freqIndx)=min(tmp(1,1),tmp(2,2));
        
        Sr_eva(1,freqIndx)=max(tmp2(1,1),tmp2(2,2));
        Sr_eva(2,freqIndx)=min(tmp2(1,1),tmp2(2,2));
              
        S_GC_eva(1,freqIndx)=max(tmp3(1,1),tmp3(2,2));
        S_GC_eva(2,freqIndx)=min(tmp3(1,1),tmp3(2,2));
              
        
                    
        S_GC1_eva(1,freqIndx)=max(tmp4(1,1),tmp4(2,2));
        S_GC1_eva(2,freqIndx)=min(tmp4(1,1),tmp4(2,2));
        
    end
    
    %% plot
    clf
    l = 5;
    w =4;
    cntr = 1;
    ax = [];
    
    set(0,'DefaultAxesFontSize',20) 
    set(0,'DefaultLineLineWidth',2)
 
    %
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,log(S_eva(1,:)))
    
    ylim([-2.5 3]);
    
    ylabel('S')
    title('log(eig 1)')
    
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,log(S_eva(2,:)))
    
    ylim([-2.5 3]);
    
    title('log(eig 2)')
    
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,log(sum(S_eva)))
    
    ylim([-2.5 3]);
    
    title('log(eig 2)+log(eig 1) (=log(det(S))')
    
    % leave blank
     cntr = cntr + 1;
    % GC
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,log(S_eva(1,:)))
    
    hold on
    plot(freqs,log(S_GC_eva(1,:)),'--')
    
    ylim([-2.5 3]);
    ylabel('S__GC X_2 to X_1');
    
    
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,log(S_eva(2,:)))
    hold on
    plot(freqs,log(S_GC_eva(2,:)),'--')
    ylim([-2.5 3]);
    
    
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,sum(log(S_eva)),freqs,sum(log(S_GC_eva)),'r--')
%     title('log(eig 2)+log(eig 1) = log(det(S''))')
    
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,sum(log(S_GC_eva))-sum(log(S_eva)))
    title('GC X_2 to X_1');
    
    
    % GC1
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,log(S_eva(1,:)))
    hold on
    plot(freqs,log(S_GC1_eva(1,:)),'--')
    ylim([-2.5 3]);
    
    ylabel('S__GC X_1 to X_2');
    
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,log(S_eva(2,:)))
    hold on
    plot(freqs,log(S_GC1_eva(2,:)),'--')
    ylim([-2.5 3]);
    
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,sum(log(S_eva)),freqs,sum(log(S_GC1_eva)),'r--')
%     title('log(eig 2)+log(eig 1) = log(det(S''))')
    
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,sum(log(S_GC1_eva))-sum(log(S_eva)))
    title('GC X_1 to X_2');
    
    % phi
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,log(S_eva(1,:)))
    hold on
    plot(freqs,log(Sr_eva(1,:)),'--')
    ylabel('S__r') 
    
    ylim([-2.5 3]);
    line([190 190],ylim)
    
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,log(S_eva(2,:)))
    hold on
    plot(freqs,log(Sr_eva(2,:)),'--')
    ylim([-2.5 3]);
    line([190 190],ylim)
    
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,sum(log(S_eva)),freqs,sum(log(Sr_eva)),'r--')
%     title('log(eig 2)+log(eig 1) = log(det(S''))')
    line([190 190],ylim)

    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,sum(log(Sr_eva))-sum(log(S_eva)))
    title('Spectral Phi');
    line([190 190],ylim)
    
    % coh
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,log(S_eva(1,:)))
    hold on
    plot(freqs,squeeze(log(S(1,1,:))),'--')
    ylim([-2.5 3]);
    line([190 190],ylim)
    ylabel('S__multi')
    
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,log(S_eva(2,:)))
    hold on
    plot(freqs,squeeze(log(S(2,2,:))),'--')
    ylim([-2.5 3]);
    line([190 190],ylim)
    
    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,sum(log(S_eva)),freqs,squeeze(log(S(1,1,:)))+squeeze(log(S(2,2,:))),'r--')
%     title('log(eig 2)+log(eig 1) = log(det(S''))')
    line([190 190],ylim)

    ax(cntr)=subplot(l,w,cntr); cntr = cntr + 1;
    plot(freqs,squeeze(log(S(1,1,:)))+squeeze(log(S(2,2,:)))-sum(log(S_eva))')
    title('Spectral multi info');
    line([190 190],ylim)
    
    set(ax,'XLim',[0 500])
    
    svFldr = ['/Volumes/data/Programs/Sasai/FIGURES/' systems(systemIndx).nm 'eigens'];
    saveFigAs(svFldr,0,'halfScreen')
    
    
    %% Have a detailed look at tf      
    H = var2trfun(A,freq_res);
    
    %ampliture
    clf
    cntr = 1;
    set(0,'DefaultAxesFontSize',20) 
    set(0,'DefaultLineLineWidth',2)
    subplot(2,2,cntr); cntr = cntr+1;
    plot(freqs,squeeze(abs(H(1,1,:))).^2);
    title('From X1 to X1');
    ylabel('Amplitude');
    subplot(2,2,cntr); cntr = cntr+1;
    plot(freqs,squeeze(abs(H(1,2,:))).^2);
    title('From X2 to X1');
    subplot(2,2,cntr); cntr = cntr+1;
    plot(freqs,squeeze(abs(H(2,1,:))).^2);
    title('From X1 to X2');
    xlabel('Frequency');
    ylabel('Amplitude');
    subplot(2,2,cntr); cntr = cntr+1;
    plot(freqs,squeeze(abs(H(2,2,:))).^2);
    title('From X2 to X2');
    xlabel('Frequency');
    suptitle('abs(H)^2')
    
    svFldr = ['/Volumes/data/Programs/Sasai/FIGURES/' systems(systemIndx).nm 'tfAmp'];
    saveFigAs(svFldr,0,'halfScreen')
    
    clf
       
    cntr = 1;
    figure
    subplot(2,2,cntr); cntr = cntr+1;
    plot(freqs,squeeze(unwrap(angle(H(1,1,:)))));
    title('From X1 to X1');
    ylabel('Phase');
    subplot(2,2,cntr); cntr = cntr+1;
    plot(freqs,squeeze(unwrap(angle(H(1,2,:)))));
    title('From X2 to X1');

    subplot(2,2,cntr); cntr = cntr+1;
    plot(freqs,squeeze(unwrap(angle(H(2,1,:)))));
    title('From X1 to X2');
    xlabel('Frequency');
    ylabel('Phase');

    subplot(2,2,cntr); cntr = cntr+1;    
    plot(freqs,squeeze(unwrap(angle(H(2,2,:)))));
    title('From X2 to X2');
    xlabel('Frequency');
    suptitle('unwrap(angle(H))')
    
    svFldr = ['/Volumes/data/Programs/Sasai/FIGURES/' systems(systemIndx).nm 'tfPhase'];
    saveFigAs(svFldr,0,'halfScreen')
    
%     plot(freqs,squeeze(abs(H(2,1,:)))',freqs,squeeze(abs(H1(2,1,:)))')
        
end

    

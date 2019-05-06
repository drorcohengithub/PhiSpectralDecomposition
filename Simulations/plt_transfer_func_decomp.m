function [ fig_h ] = plt_transfer_func_decomp( tf,SIG, freqs, log_det_S_f,ylims )

clf
%%
subplot(4,1,1)
hold on
comp_names = {};
for i = 1:size(tf,1)
    for j =  1:size(tf,2)
        comp_names = {comp_names{:} [num2str(i) ' ' num2str(j)]};
        h =  plot(freqs,log(squeeze(abs(tf(i,j,:)).^2)));
        if i ~= j
            set(h,'linestyle','-');
        end
    end
end
legend(comp_names)
title('log det tf abs^2')
ylim(ylims(1,:))
%%

subplot(4,1,2)
log_det_SIG = ones(length(freqs),1)*log(det(SIG));
plot(freqs,log_det_SIG);
ylim(ylims(2,:))
title('log det SIG')

%%
subplot(4,1,3)
log_det_S_r = log(squeeze(abs(tf(1,1,:)).^2)) + log(squeeze(abs(tf(2,2,:)).^2)) + log_det_SIG;
plot(freqs,log_det_S_r);
hold on
plot(freqs,log_det_S_f,'r');

ylim(ylims(3,:))
title('log_det_S_r')

if size(log_det_S_r)
%%
subplot(4,1,4)
ci_f = (log_det_S_r-log_det_S_f)/2;
plot(freqs,ci_f);
ylim(ylims(4,:))
title('ci_f')

fig_h = gcf;
set(fig_h, 'Position', [1    82   960   903])


end


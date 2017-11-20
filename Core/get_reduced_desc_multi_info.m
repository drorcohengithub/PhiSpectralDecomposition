function [S_r,det_S_r,SIG_r,det_SIG_r, A_r] = get_reduced_desc_multi_info(G,blocks,fres)


% time domain
[SIG_r A_r As SIGs] = autocov_to_blockvar(G,blocks);
det_SIG_r = det(SIG_r);
%freq
S = autocov_to_cpsd(G,fres);
S_r = zeros(size(S));

for blockIndx = 1:length(blocks)
    S_r(blocks{blockIndx},blocks{blockIndx},:) = S(blocks{blockIndx},blocks{blockIndx},:);
end
    
%check spectal identity
for freqIndx = 1:size(S_r,3)
    thisFreqS = squeeze(S_r(:,:,freqIndx));
    det_S_r(freqIndx) = det( thisFreqS );   
end

diff = abs(log(det_SIG_r) - mean(log([det_S_r det_S_r(2:end-1)])));
if  diff > 1e-14
    warning(['spectral - time equivalence may be violated. Diff is: ' num2str(diff)])  
end


end


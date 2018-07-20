function [ratio_S ratio det_S] = ratio_of_dets(S, S_r, SIG, SIG_r)

% logdet(S)
det_S_r = zeros(size(S_r,3),1);
det_S = zeros(size(S_r,3),1);

for freqIndx = 1:size(S_r,3)
    this_freq_S_r = squeeze(S_r(:,:,freqIndx));
    det_S_r(freqIndx) = det( this_freq_S_r );   
  
    this_freq_S = squeeze(S(:,:,freqIndx));
    det_S(freqIndx) = det( this_freq_S );  

end

ratio = log(det(SIG_r)/det(SIG));
ratio_S = log(det_S_r./det_S);

%check spectal identity
diff = abs(ratio - mean([ratio_S;ratio_S(2:end-1)]));
if  diff > 1e-14
    warning(['spectral - time equivalence may be violated. Diff is: ' num2str(diff)])  
end


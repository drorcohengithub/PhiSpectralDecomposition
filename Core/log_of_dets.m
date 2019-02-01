function [log_det_S] = log_of_dets(S)

% logdet(S)
det_S = zeros(size(S,3),1);

for freqIndx = 1:size(S,3)

    this_freq_S = squeeze(S(:,:,freqIndx));
    det_S(freqIndx) = det( this_freq_S );  

end

log_det_S = log(det_S);


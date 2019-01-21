function [ log_det_S ] = log_dets(S)

numfreqs = size(S,3);
log_det_S = nan(size(S,3),1);
for freq = 1:numfreqs
    
    log_det_S(freq) = log( det( squeeze(S(:,:,freq)) ) );
    
end


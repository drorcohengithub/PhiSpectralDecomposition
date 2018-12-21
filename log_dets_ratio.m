function [ ratio ] = log_dets_ratio(numerator,denominator)

numfreqs = size(numerator,3);
ratio = zeros(1,numfreqs);

for freq = 1:numfreqs
    
    ratio(freq) = log( det(squeeze(numerator(:,:,freq))) / det(squeeze(denominator(:,:,freq))) );
    
end


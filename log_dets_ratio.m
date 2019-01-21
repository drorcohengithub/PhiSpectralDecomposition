function [ ratio ] = log_dets_ratio(numerator,denominator)

log_det_numerator = log_dets(numerator);
log_det_denominator = log_dets(denominator);

ratio = log_det_numerator ./ log_det_denominator;


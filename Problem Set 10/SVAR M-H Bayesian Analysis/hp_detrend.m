% HP_DETREND - Return the detrended version of a variable using the HP-filter
%
% usage:
%
% [detrended, trend] = hp_detrend(x, lam, trans)
%
% where
%
% x = the data to be detrended
% lam = the smoothing parameter
% trans = the trans formation to be performed ('pct_dev', 'difference', ...)
%
% detrended = the detrended variable
% trend = the filtered trend.

function [detrend, trend] = hp_detrend(series, lam, varargin)


%Get the trend

if sum(isnan(series))
    warning ('HP filter with NaN; valid entries must be contiguous');
    idx = ~isnan(series); 
    ser_tmp = series(idx);
    trend_tmp = hpfilter(ser_tmp, lam);
    trend = zeros(size(series));
    trend(idx) = trend_tmp;
else
    trend = hpfilter(series, lam);
end
%Do the proper transformation
if length(varargin) > 0
    trans = varargin{1};
else
    trans = 'difference';
end


if strcmp(trans, 'pct_dev')
    detrend = 100*((series - trend)./trend);
elseif strcmp(trans,'difference')
    detrend = series - trend;
end
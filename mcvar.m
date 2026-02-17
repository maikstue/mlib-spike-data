function cv = mcvar(x,mode)
% function cv = mcvar(x,mode)
% 
% computes the coefficient of variation cv for vector x
% 
% INPUT
% x       vector of spike times or interspike intervals (both should be provided in seconds)
%           - if x is a vector of spike times (usually spike times relative to the beginning of the recording until the end), set mode to 'diff'
%             the input vector will be differentiated, so you get cv = std(diff(x)) / mean(diff(x))
%           - if x is a vector of interspike intervals, set mode to 'abs'
%             the computation is then simply cv = std(x) / mean(x)
% 
% If you want to know the CV of a vector of spike counts x (e.g., spikes per s in consecutive 1-s bins), use mcvar(x,'abs').
% The result is called the Fano factor.
% 
% OUTPUT
% cv        coefficient of variation, std(x)/mean(x) OR std(diff(x))/mean(diff(x))
% 
% Maik C. Stüttgen, Rotterdam, May 2013
%% check input and compute cv
if isempty(x)
  cv = nan;
elseif numel(x)~=max(size(x))
  error('vector input please')
else
  switch mode
    case 'diff'
      cv = std(diff(x))/mean(diff(x));
    case 'abs'
      cv = std(x)/mean(x);
  end
end
end
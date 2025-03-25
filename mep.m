function [ep,errb,ci] = mep(fp,tevents,pre,post,varargin)
% function [ep,errb,ci] = mep(fp,tevents,pre,post,varargin)
% computes and plots an evoked = event-related potential from continuous LFP data
% optionally returns tick-wise error bars (SEM) and confidence intervals
% 
% INPUTS
% fp        field potential recording sampled at arbitrary frequency
% tevents   time stamps of events (in ticks)
% pre       number of ticks pre-event
% post      number of ticks post-events
% ciwidth   width of confidence interval (default: 95%); optional
% 
% OUTPUTS
% ep        evoked potential (averaged over trials, ranging from -pre to +post ticks relative to tevents)
% errb      standard error of the mean for each tick
% ci        confidence intervals (lower, upper)
% 
% Maik C. Stüttgen, October 2024
%% check varargin
ciwidth = 95;
if ~isempty(varargin)
  for i=1:2:size(varargin,2)
    switch varargin{i}
      case 'ciwidth'
        ciwidth = varargin{i+1};
      otherwise
        errordlg('unknown input argument')
    end
  end
end

%% compute evoked potential
epmat = nan(numel(tevents),pre+post+1);

for i = 1:numel(tevents)
  try
    epmat(i,:) = fp(tevents(i)-pre:tevents(i)+post,:);
  catch
    epmat(i,:) = nan;
  end
end
ep = nanmean(epmat);

%% compute standard error of the mean
if nargout>1
  errb = std(epmat,0,1)/sqrt(numel(tevents));
end

%% compute confidence intervals
if nargout>2
  ci = prctile(fp,[0.5*(100-ciwidth),ciwidth+0.5*(100-ciwidth)]);
end

end
function nspx = mnspx(spxtimes,trigtimes,pre,post,removenans)
% nspx = mnspx(spxtimes,trigtimes,pre,post,removenans)
% function simply returns the distribution of spike counts within a specific interval relative to a trigger
% IMPORTANT: all timestamp inputs (spxtimes, trigtimes) must be seconds and will be converted to ms in the script!
%
% MANDATORY INPUTS
% spxtimes      vector with timestamps (seconds) of spike events
% trigtimes     vector with timestamps (seconds) of trigger events
% pre           time before trigger to include milliseconds; default 1000 ms)
% post          time after trigger to include (milliseconds; default 1000 ms)
% removenans    remove triggers which are NaN; 'true' (default) or 'false'
% 
% pre and post are inclusive, so pre=200 and post=200 will count spikes from >=-200 to <=200 ms peri-stimulus
% 
% sorry about the two different timing formats... some compatibility issue forces me to do that... feel free to edit below!
%
% EXAMPLE
% get (and plot) the distribution of spike counts relative to the start of event 9 (food reward) and the following 2000 ms
%     nspx = mnspx(spx.timings,spx.eventtimings(spx.eventmarkers==9),0,2000)
%     hist(nspx,0.5:max(nspx)+0.5)
% 
% HISTORY
% 2022-02-12    removal of NaNs in trigtimes is optional
% 2022-01-22    NaNs in spxtimes and trigtimes are removed
% 
% by Maik C. Stuettgen, Summer 2013 @ Erasmus MC Rotterdam, The Netherlands
%% preps
spxtimes(isnan(spxtimes)) = [];
if exist('removenans','var')
  if removenans==1,trigtimes(isnan(trigtimes)) = [];
  end
elseif ~exist('removenans','var'),trigtimes(isnan(trigtimes)) = [];  
end
spxtimes  = spxtimes*1000;
trigtimes = trigtimes*1000;
nspx      = nan(numel(trigtimes),1); % preallocate for speed
%% the works
% for every trigtime, get the number of spikes in the relevant time window
for i = 1:numel(trigtimes)
  if ~isnan(trigtimes(i))
    nspx(i,1) = sum(spxtimes>=(trigtimes(i)-pre) & spxtimes<=(trigtimes(i)+post));
  end
end

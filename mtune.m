function [avFR,semFR,fitobject] = mtune(tspx,tevents,cevents,tpost)
% function [avFR,semFR,fitobject] = mtune(tspx,tevents,cevents,tpost)
% constructs and fits a tuning curve
% 
% INPUTS
% tspx        spike time stamps (in seconds)
% tevents     event time stamps (in seconds)
% cevents     event codes (integers)
% tpost       time post event onset to analyze (in ms)
% 
% OUTPUTS
% avFR        average firing rates for each stimulus
% semFR       standard errors for each stimulus
% fitobject   Gaussian fit to the tuning curve
% 
% Maik C. Stüttgen, November 2024
%% the works
stimIDs = unique(cevents);

for i = 1:numel(stimIDs)
  nspx     = mnspx(tspx,tevents(cevents==stimIDs(i)),0,tpost);
  avFR(i)  = mean(nspx/(tpost/1000));
  semFR(i) = std(nspx/(tpost/1000))/sqrt(numel(nspx));
end
clear nspx

%% generate a plot
figure('name','Tuning curve','units','normalized','position',[.4,.4,.2,.3]),hold on
errorbar(stimIDs,avFR,semFR,'o','CapSize',0)
xlabel('Stimuli'),ylabel('Spikes per s')
axis([stimIDs(1)-0.5,stimIDs(end)+0.5,min(floor(avFR-semFR))-1,max(ceil(avFR+semFR))+1])

% fit with Gaussian curve if the Curve Fitting Toolbox is available
if contains(struct2array(ver),'Curve Fitting Toolbox')
  fitobject = fit(stimIDs,avFR','gauss1');
  x = stimIDs(1):0.01:stimIDs(end);
  y = fitobject.a1*exp(-((x-fitobject.b1)/fitobject.c1).^2);
  plot(x,y)
else
  disp('Curve Fitting Toolbox not available.')
end

end

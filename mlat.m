function [mean1stspike,median1stspike,peaklat,meantspike,mediantspike] = mlat(tspx,tstims,texclude)
% function [mean1stspike,median1stspike,peaklat,meantspike,mediantspike] = mlat(tspx,tstims,texclude)
% 
% Computes different measures of response latency after stimulation.
% 1) Mean first-spike latency.
% 2) Median first-spike latency.
% 3) Time of peak firing rate.
% 4) Mean spike time in time window.
% 5) Median spike time in time window.
% 
% Which measure is the most appropriate? Sakata & Harris, Neuron 2009, Supplemental Figure 11:
% In the main text (Figure 4), we assessed the timing of activity onset in each layer using the median MUA spike time in a 50-ms time window after 
% event onset for each channel. Here we verify that these observations are robust to other definitions of peak latency.
% We note that first-spike timing is not a suitable measure, as it is biased by firing rate (A. Peyrache et al, SFN abstract 690.13, 2008;
% Supp Fig S2 of Luczak et al., 2009). Instead, we use two other measures: mean spike time, and the time of PETH (peri-event time histogram) maxima.
% In both cases, we obtain a similar profile to that seen in Figure 4.
% (A) Here we define latency for each channel and each event as the mean MUA spike time in a 50-ms time window after event onset for each channel.
% (B) Here we defined latency for each channel and each experiment as the time of the PETH peak. For each experiment, we compute the stimulus- or 
% upstate-triggered MUA PETH for all channels (smoothed with a 5-ms Gaussian kernel), and compute for each channel the time at which the PETH attains 
% its peak value in a 50-ms time window after event onset. With this measure, we obtain a single value per experiment, rather than per event.
% 
% INPUT
% tspx        spike timestamps (in seconds)
% tstims      stimulus timestamps (in seconds)
% texclude    maximum time window after stimulation in which spikes are considered (in seconds); set to inf to take all spikes
% 
% OUTPUT
% see text
% 
% Maik C. Stüttgen, University Medical Center Mainz, June 2021
%% get mean and median first-spike latencies
tFirstSpikes = nan(numel(tstims),1);
for i = 1:numel(tstims)
  try tFirstSpikes(i) = tspx(find(tspx>tstims(i),1,'first')) - tstims(i);
  end
end
tFirstSpikes(tFirstSpikes>texclude) = []; % eliminate spike latencies later than texclude; otherwise, crap results
mean1stspike   = mean(tFirstSpikes)*1000;
median1stspike = median(tFirstSpikes)*1000;

%% get PSTH peak time
psth = mpsth(tspx,tstims,'pre',200,'post',max([200,1000*(texclude)+50]),'fr',1);
sdf  = msdf(psth(:,2),'Gauss',5,'notext');
[~,idx] = max(sdf(202:(201+round(texclude*1000)))); % peak index from +1 to +50 is 251
peaklat = psth(idx+201,1);
% if sdf(idx+201)<3*mean(sdf(101:200))
%   disp('WARNING: SDF peak height is less than 3 times the pre-stimulus firing rate')
% end

%% get mean and median spike times in analysis interval
[psth,trialspx] = mpsth(tspx,tstims,'pre',0,'post',1000*texclude,'fr',0);
meantspike   = sum(psth(2:end,1).*psth(2:end,2))/sum(psth(2:end,2));  % wmean(psth(2:end,1),psth(2:end,2)) also works
mediantspike = median(cell2mat(trialspx));

end
function c = mcheck(waveforms,varargin)
% c = mcheck(waveforms,varargin)
% this function performs a quality check on a set of sorted waveforms
%
% After collecting multiple units, use mcheck_overview to plot distributions of most important indices
% 
% MANDATORY INPUT
% waveforms     a matrix whose rows are (putative) spikes and whose columns are samples (or vice versa -
%               mcheck.m assumes that the longer dimension represents separate waveforms)
%               so each row should be n samples from the same spike
%
% OPTIONAL INPUTS
% timevec       vector of timestamps of the unit in the channel; should be of the same size as waveform
%               NOTE THAT TIMESTAMPS SHOULD BE GIVEN IN SECONDS!
% markvec       vector of markers of the units in the channel; should be of the same size as waveforms
% target        which marker should be analyzed
% fs            sampling rate; reciprocal of bin size (?s) of analog sampling (e.g. 20 kHz for 50 ?s)
%               if provided, mcheck will compute spike width (FWHM, full width at half maximum)
% tevents       vector with event times; if provided, an additional plot will be generated that depicts
%               spike times relative to event times +- 200 ms; in addition, spike waveforms in a window +- 20 ms will be plotted to verify that these
%               waveforms are indeed spikes and not some event- or stimulus-related artifact
%               if tevents is provided, timevec is also needed
%               NOTE THAT TIMESTAMPS SHOULD BE GIVEN IN SECONDS!
% plotno        number of waveforms to plot; default is 5,000
% elimWin       window for elimination of spikes in vicinity of a specified event(in ms; default [20 20])
% plotit        if set to '1' (default), a plot is generated; if set to '0', the plot is automatically closed
%
% OUTPUT c      structure with a range of quality indices (mean and sd of waveforms, SNR etc)
%
% HISTORY
% Oct 2024            polished for publication
% Feb 2020            started new version (this one); deleted SNR99 and SNR95; SNRd values now use the 1%-winsorized noise distribution
%                     SNR95 correlated with -SNR_dmin+SNR_dmax at r=0.99 in Jens' data set
%                     refractory period is now limited to 2 ms (should be <1.5%; Reyes-Puerta et al., Cerebral Cortex 2015)
% Jan 2020            if we have >1000 spikes, only every 10th spike amplitude is plotted in subplot(4,4,6)
%                     text output changed for refractory period violations from <4 ms to <3 ms; the former was correct when
%                     considering integer values only
%                     mcheck now returns all parameters from mwave in the field 'waveparms'; the fields FWHM and FWHMax have been deleted
%                     [fwhmMax,fwhmMin,p2tAmp,p2tDur,p2tRat]
%                     I found that SNR95 and SNR99 are correlated by r=0.98 in an example data set!
%                     SNR95 more closely matches visual impression with density plots than SNR_dmin.
% Jul 2019            added measure for unit stability - coefficient of variation of spike counts in 100-s bins
%                     added the same for waveform amplitudes - std(amps)/mean(amps)
% Jan 2017            replaced sum(spxtimes) with numel(spxtimes)
% Oct 2016            added cell2vectors as nested function
% Sep 2016            subplot(445): added value of rank-biserial correlation b/w time and amplitude as well as its significance
%                     added the main diagonal and a horizontal line (mean of spike amplitudes) as visual guides to assess stability
%                     amplitudes in this subplot are now smoothed with a kernel of 9 values (previously: 5)
%                     subplot(446): rescaled left and right axes
%                     rearranged some variables (timevec, waveforms) and code (quality criterion values)
%                     changed name of variable binsz to fs for sampling rate
%                     replaced xlabel and ylabel with set(get(gca,'XLabel'),'String','axisname') to increase function speed
% May 2015            adjustments for Matlab R2015a
%                     speeded up spike elimination
% Oct 2014            added optional argument 'plotit' to close the plot if no plot is desired
% Mar 20, 2014        added total spike count and firing rate to subplot 444
%                     debugged subplot 445 (amplitude over time)
% Feb 6, 2014         renamed checkunit.m to mcheck.m - for consistency
% Feb 5, 2014         adjusted quality criteria to match, assuming Gaussian distributions of signal and noise
% Oct 30, 2013        changed all single-precision arrays to double-precision
% June 25, 2013       added voltage over time in subplot(445)
% July 3, 2012        renamed variable 'time', included estimate of overall firing rate
% June 28, 2012       some polishing work for publication
% June 8, 2012        renamed variable tpecks to tevents
% June 4, 2012        changed output figure size; changed critical values
% May 7, 2012         new outputs: estimated percent false negatives, percent refractory period violations
%                     distill recommendations for classification (take e.g. skew and std)
% April 27, 2012      stepsize for histograms calculated automatically
%                     debugged histogram output for subplot 447
% Aug 16, 2011        changed sd plot to % of sd
% June 16, 2011       new additional outputs:
%                     - no. of waveforms within event range
%                     - zoom in on ISI distribution's first 50 ms
%                     - spike amplitude over session time plus regression line
%                     - firing rate over session time
%                     - Gaussian fit and false negative estimation added optional input: elimWin
% June 16, 2011       cosmetic fixes to subplot 3
% May 24, 2011        added density plot and waveform samples; incorporated cpecks.m code into this function
% May 04, 2011        constrained number of plotted waveforms to 5,000 (default)
% Feb 28, 2011        added cpecks and binsz (optional arguments)
% Feb 16, 2011        added examples
%
% by Maik C. St?ttgen, Feb 2011
%% assign criteria for traffic lights (green, yellow, red)
% criteria are specified such that, if signal and noise distributions are Gaussian and the signal is symmetrical in terms
% of amplitude, they should yield the same recommendations; these assumptions are usually violated but they are at least
% internally consistent
crit_SNRd  = [6,4];            % SNR(dmax) = mean maximum value minus mean noise value divided by 1%-winsorized SD of noise
                                  % so d measures SNR from noise to either maximum or minimum value; sum them up to get the whole range!
crit_SNRdf = [10,8];           % SNR_dmax - SNR_dmin
crit_skew  = [0.5,0.75];
crit_ISI   = [0.015,0.03];
crit_FN    = [0.15,0.30];
crit_CVfr   = [0.4,0.8];
crit_CVa   = [0.1,0.2];
plotno  = 5000;
plotit  = 1;
elimWin = [20 20];
minInt  = 0.000; % minimum time interval between two consecutive events (in seconds!); events following the preceding events within <minInt ms will be discarded
nospx   = 4;     % duration of absolute refractory period in ms
%% get input data and check it
if max(size(waveforms))<64
  disp('less than 64 waveforms detected, aborting...')
  return
end
if size(waveforms,1)<size(waveforms,2)
  waveforms=waveforms';
end
target  = 1;
fs      = [];
if nargin
  for i=1:2:size(varargin,2)
    switch varargin{i}
      case 'markvec'
        markvec = varargin{i+1};
      case 'target'
        target = varargin{i+1};
      case 'timevec'
        spxtimes = varargin{i+1};
      case 'fs'
        fs = varargin{i+1};
      case 'plotno'
        plotno = varargin{i+1};
      case 'tevents'
        tevents  = varargin{i+1};
      case 'elimWin'
        elimWin = varargin{i+1};
      case 'plotit'
        plotit  = varargin{i+1};
      otherwise
        errordlg('unknown argument')
        return
    end
  end
end
clear varargin
%% rearrange input according to markvec and target
if exist('markvec','var')
  if numel(markvec)~=size(waveforms,1)
    errordlg('vectors do not match');return
  end
  wavs     = double(waveforms(markvec==target,:));
  try
    spxtimes = double(spxtimes(markvec==target));
  end
  clear markvec waveforms
elseif ~exist('markvec','var')
  wavs = double(waveforms);
  clear waveforms
end
%% calculate some indices...
c.spx       = size(wavs,1);
c.meanwave  = mean(wavs,1);
c.stdev     = std(wavs,0,1);
c.skew      = skewness(wavs,0,1);
if exist('spxtimes','var')
  c.recdur_spx    = (spxtimes(end)-spxtimes(1))/60;      % recording duration as inferred from the first to the last spike
  c.fr            = 1/mean(diff(spxtimes));
else
  c.recdur_spx    = nan;
  c.fr            = nan;
end
if exist('tevents','var')
  c.recdur_events = (tevents(end)-tevents(1))/60;        % recording duration as inferred from the first to the last event
else
  c.recdur_events = nan;
end
%% plot 'em
figure('units','normalized','position',[0.05 .1 .8 .8],'name',['Unit quality check for marker ',num2str(target),...
  ' with ',num2str(c.spx),' waveforms, overall firing rate is ' num2str(c.fr,'%2.1f') ' Hz'])
%% subplot 1 - raw waveforms
subplot(451)
title(['up to ' num2str(plotno) ' waveforms'],'FontSize',10),hold all
if plotno>size(wavs,1),plotno=size(wavs,1);end
plot(wavs(randsample(1:size(wavs,1),plotno),:)','b')
plot(c.meanwave,'k','LineWidth',2);
plot(c.meanwave+c.stdev,'k');
plot(c.meanwave+2*c.stdev,'k');
plot(c.meanwave-c.stdev,'k');
plot(c.meanwave-2*c.stdev,'k');
lower = min(c.meanwave)-mean(c.stdev)*4;
upper = max(c.meanwave)+mean(c.stdev)*4;
axis([0 size(wavs,2)+1 lower upper])
set(get(gca,'XLabel'),'String','Time (bins)','FontSize',10);set(get(gca,'YLabel'),'String','ADC units','FontSize',10); % much faster than xlabel and ylabel
%% subplot 2 - waveform variability
subplot(452)
title('Waveform variability','FontSize',10),hold all
[ax,h1,h2] = plotyy([5 size(wavs,2)-5],[100 100],[5 size(wavs,2)-5],[mean(c.skew) mean(c.skew)]);
plot(ax(1),1:size(wavs,2),100*c.stdev/mean(c.stdev(1:5)),'b')
plot(ax(1),[5 size(wavs,2)-5],[50 50],'b:')
plot(ax(1),[5 size(wavs,2)-5],[150 150],'b:')
set(gcf,'CurrentAxes',ax(2))
line(1:size(wavs,2),c.skew,'Color','r') % plot seems to destroy parent figure
line([5 size(wavs,2)-5],[-.5 -.5],'Color','r','LineStyle',':')
line([5 size(wavs,2)-5],[.5 .5],'Color','r','LineStyle',':')
set(get(ax(1),'Xlabel'),'String','Time (bins)','FontSize',10)
set(get(ax(1),'Ylabel'),'String','% of noise SD','FontSize',10)
set(ax(1),'Xlim',[0 size(wavs,2)+1],'Ylim',[0 200],'YTick',[0 100 200],'YColor','b')
set(h1,'Color','b')
set(get(ax(2),'Ylabel'),'String','Skewness','FontSize',10)
yskew = [-1 1];
set(ax(2),'Xlim',[0 size(wavs,2)+1],'Ylim',[-1 1],'YTick',[yskew(1) 0 yskew(2)],'YColor','r')
set(h2,'Color','r')
clear ax h1 h2 yskew
%% subplot 3 - distributions of peak and baseline voltages
subplot(453)
% adjust stepsize for histograms according to mean waveform
step = (max(mean(wavs))-min(mean(wavs)))/50;
title('Baseline & peak voltages','FontSize',10),hold on
% baseline voltages
[n,b] = hist(wavs(:,1),min(wavs(:,1)):step:max(wavs(:,1)));
bar(b,n,'b','EdgeColor','b'),hold on
% peak voltages
a1 = max(wavs,[],2);
a3 = min(wavs,[],2);
[n2,b2] = hist(a1,min(a1):step:max(a1));
[n4,b4] = hist(a3,min(a3):step:max(a3));
bar(b2,n2,'r','EdgeColor','r'),bar(b4,n4,'g','EdgeColor','g')
set(get(gca,'XLabel'),'String','ADC units','FontSize',10);
set(gca,'XLim',[lower upper])
clear b n step
%% subplots 4 & 5 - estimating false negatives for minimum and maximum peak distributions
% guess mode - inspired by mode_guesser.m to find the range of the most tightly distributed data use median of the tightest range as the guess of the mode have
% to do for minimum and maximum!
p = 0.5;
x = a3;
num_samples = length(x);
shift = round(num_samples*p);
x = sort(x);
[~,m_spot] = min(x(shift+1:end) - x(1:end-shift));
m = x(round(m_spot + (shift/2)));
% now duplicate sorted data set around the mode
if skewness(x)>0 % skewed towards positive values
  y = vertcat(m+m-x(x>m),x(x>m));
else
  y = vertcat(m+m-x(x<m),x(x<m));
end
Gaussmean = mean(y);
Gaussstd  = std(y);
xval         = linspace(min(y),max(y),1000);
distribution = normpdf(xval,Gaussmean,Gaussstd);
distribution = distribution/(max(distribution));
c.FNmin = 1-normcdf(max(b4),Gaussmean,Gaussstd);

subplot(454)
title('Distribution of minima','FontSize',10),hold all
bar(b4,n4,'g','EdgeColor','g')
plot(xval,distribution*max(n4),'k')
set(get(gca,'XLabel'),'String','ADC units','FontSize',10);
set(gca,'XLim',[lower Gaussmean+4*max(c.stdev)])

x = a1;
num_samples = length(x);
shift = round(num_samples*p);
x = sort(x);
[~,m_spot] = min(x(shift+1:end) - x(1:end-shift));
m = x(round(m_spot + (shift/2)));
% now duplicate sorted data set around the mode
if skewness(x)>0 % skewed towards positive values
  y = vertcat(m+m-x(x>m),x(x>m));
else
  y = vertcat(m+m-x(x<m),x(x<m));
end
Gaussmean = mean(y);
Gaussstd  = std(y);
xval         = linspace(min(y),max(y),1000);
distribution = normpdf(xval,Gaussmean,Gaussstd);
distribution = distribution/(max(distribution));
c.FNmax = 1-normcdf(-min(b2),-Gaussmean,Gaussstd);

subplot(455)
title('Distribution of maxima','FontSize',10),hold all
bar(b2,n2,'r','EdgeColor','r')
hold all
plot(xval,distribution*max(n2),'k')
set(get(gca,'XLabel'),'String','ADC units','FontSize',10);
set(gca,'XLim',[Gaussmean-4*max(c.stdev) upper])

clear distribution Gaussmean Gaussstd m m_spot p b2 n2 b4 n4 shiftxval x y xval num_samples
%% subplot 6 - waveform stability: cumulative spikes over time, amplitude over time
subplot(456)
amps = max(wavs,[],2)-min(wavs,[],2);
if exist('spxtimes','var') && numel(spxtimes)>50 && ~any(isnan(spxtimes))
  try
    title('Stability','FontSize',10),hold on
    mean_amp = nan(floor(max(spxtimes/60)),1);
    sd_amp   = nan(floor(max(spxtimes/60)),1);
    sem_amp  = nan(floor(max(spxtimes/60)),1);
    fr       = nan(floor(max(spxtimes/60)),1);
    for i = 1:floor(max(spxtimes/60))   % mean_amp and fr will be reused for subplot 10
      idx = find(spxtimes>((i-1)*60) & spxtimes<i*60);
      mean_amp(i) = mean(amps(idx));
      sd_amp(i)   = std(amps(idx));
      sem_amp(i)  = sem(amps(idx));
      fr(i)       = numel(idx)/60;
    end
    [n,b] = hist(spxtimes,0:1:max([max(spxtimes) 1]));
    nn = histcounts(spxtimes,0:60:max([max(spxtimes),1]));  % replaced hist by histcounts on Aug 9 2021, minor effect (2nd decimal place)
    % The CV of firing rates was also used by Reyes-Puerta et al. (2015, Cerebral Cortex), they considered units with CV<1.5 as stable.
    % They used 100-s bins; I will use 1-min bins
    c.cvfr  = std(nn)/mean(nn);
    c.cvamps = std(amps)/mean(amps);
    [ax,h1,h2] = plotyy(b/60,cumsum(n),1:numel(mean_amp),mean_amp);
    line([0 max(spxtimes)/60],[0 sum(n)],'Parent',ax(1),'Color','k','LineStyle','-','LineWidth',1)    % plot destroys axes, line does not
    line([0 ceil(max(b/60))],[mean(mean_amp),mean(mean_amp)],'Parent',ax(2),'Color','k','LineStyle','-','LineWidth',1)
    line(1:i,mean_amp+1*sd_amp,'Parent',ax(2),'Color','r','LineStyle','-','LineWidth',1)
    line(1:i,mean_amp-1*sd_amp,'Parent',ax(2),'Color','r','LineStyle','-','LineWidth',1)
    uistack(h1,'top')
    uistack(h2,'top')
    set(get(ax(1),'Ylabel'),'String','Cumsum(spikes) in 1-s bins','FontSize',10)
    set(ax(1),'Xlim',[0 ceil(max(b/60))],'Ylim',[0 sum(n)],'YColor','b')
    set(ax(2),'Xlim',[0 ceil(max(b/60))],'Ylim',[nanmean(mean_amp)-10*nanstd(mean_amp) nanmean(mean_amp)+10*nanstd(mean_amp)],'YColor','r')%,'YTickLabel','')
    set(h1,'LineStyle','-','Color','b','LineWidth',2)
    set(h2,'LineStyle','-','Color','r','LineWidth',2)
    [rho,pval] = corr(spxtimes,amps,'type','Spearman');
    text(ax(1),0.2*spxtimes(end)/60,0.1*numel(amps),['r=',num2str(rho,'%2.2f'),', p=',num2str(pval,'%2.3f')],'Color','r');
    set(get(gca,'XLabel'),'String','Time (min)','FontSize',10);
    clear ax h1 h2 rho pval
  catch
    axis([0 1 0 1])
    text(0,0.5,'This did not work.','FontName','Courier')
    axis off
  end
else
  c.cvamps = std(amps)/mean(amps);
  try
    plot(spxtimes/60,smoothdata(amps,'movmean',5),'r');
    xlim([0 ceil(spxtimes/60)])
    ylim([mean(amps)*0.5 mean(amps)*1.5])
    set(get(gca,'YLabel'),'String','p2p amplitude','FontColor','r','FontSize',10);
    axis off
  catch
    axis([0 1 0 1])
    text(0,0.5,'No timestamps','FontName','Courier')
    axis off
  end
end
%% subplot 7 - waveform stability: spike count histogram over time, amplitudes over time
% w = get(subplot(446),'Position'); % somehow yields different positions
% subplot('Position',[0.34264 0.54832 0.12 0.15554]) % left bottom width height
subplot(457)
if exist('spxtimes','var') && numel(spxtimes)>20 && ~any(isnan(spxtimes))
  title('Stability','FontSize',10),hold on
  [n,b] = hist(spxtimes,0:60:max(spxtimes));
  if numel(spxtimes)<=1000
    [ax,h1,h2] = plotyy(b/60,n,spxtimes/60,amps,'bar','plot');
  elseif numel(spxtimes)>1000 && numel(spxtimes)<=5000
    [ax,h1,h2] = plotyy(b/60,n,spxtimes(1:2:end)/60,amps(1:2:end),'bar','plot');   % plot only every 2nd data point
  elseif numel(spxtimes)>5000 && numel(spxtimes)<25000
    [ax,h1,h2] = plotyy(b/60,n,spxtimes(1:5:end)/60,amps(1:5:end),'bar','plot');   % plot only every 5th data point
  elseif numel(spxtimes)>25000
    [ax,h1,h2] = plotyy(b/60,n,spxtimes(1:10:end)/60,amps(1:10:end),'bar','plot'); % plot only every 10th data point
  end
  set(h1,'FaceColor','k')
  set(gcf,'CurrentAxes',ax(2))
  pars = polyfit(spxtimes/60,amps,1);
  line(spxtimes/60,pars(1)*spxtimes/60+pars(2),'Color','k') % plot seems to destroy parent figure
  set(gcf,'CurrentAxes',ax(1))
  text(0.1*max(b)/60,0.2*max(n),['a=',num2str(pars(1),'%1.2f'),', b=',num2str(pars(2),'%1.2f')],'Color','r','FontSize',12,'FontWeight','bold')
  set(h1,'LineStyle','-')
  set(h2,'LineStyle','none','Marker','.','MarkerSize',3)
  set(get(ax(1),'Xlabel'),'String','Time (min)','FontSize',10)
  set(get(ax(1),'Ylabel'),'String','N spikes in 60-s bins','FontSize',10)
  set(ax(1),'Xlim',[-1 ceil(max(b/60))+1],'Ylim',[0 max(n)*3],'YTick',[0,round(max(n))],'YTickLabel',{'0',round(max(n))})
  set(ax(2),'Xlim',[-1 ceil(max(b/60))+1],'Ylim',[min(amps)-1*range(amps) max(amps)+range(amps)/10],'YColor','r','YTickLabel','')
  set(h2,'Color','r')

  clear b n ax h1 h2 pars
else
  text(0,0.5,'No timestamps','FontName','Courier')
  axis off
end
%% subplot 8 - waveform amplitude vs. firing rate
subplot(4,5,8),hold on
try
  shades = repmat(linspace(0.2,0.8,numel(mean_amp)),3,1)';
  [c.rAmpFr,c.pAmpFr] = corr(mean_amp,fr,'rows','complete');
  for i = 1:numel(mean_amp)
    plot(fr(i),mean_amp(i),'.','Color',shades(i,:))
  end
  pars = polyfit(fr,mean_amp,1);          % will not work with NaNs but mcheck will not crash
  line(fr,pars(1)*fr+pars(2),'Color','k') % plot seems to destroy parent figure
  set(gca,'XLim',[min(fr)-std(fr) max(fr)+std(fr)],'YLim',[min(mean_amp)-nanstd(mean_amp) max(mean_amp)+nanstd(mean_amp)])
  axx = get(gca,'XLim');
  axy = get(gca,'YLim');
  text(axx(1)+0.1*(axx(2)-axx(1)),axy(2),['r=',num2str(c.rAmpFr,'%2.2f'),', p=',num2str(c.pAmpFr,'%2.3f')],'Color','r');
  set(get(gca,'XLabel'),'String','Spikes per s','FontSize',10);set(get(gca,'YLabel'),'String','Amplitude','FontSize',10);
  clear axx axy idx fr pars shades
catch
  axis off
end
%% subplots 9 & 10 - waveform stability
try
  if numel(spxtimes)>50 && ~any(isnan(spxtimes))
    subplot(459),hold on
    mean_amp_5min = nan(floor(max(spxtimes/300)),1);   % 5-min segments
    std_amp_5min  = nan(floor(max(spxtimes/300)),1);
    sem_amp_5min  = nan(floor(max(spxtimes/300)),1);
    varerr_amp_5min  = nan(floor(max(spxtimes/300)),1);
    all_amps_5min = cell(floor(max(spxtimes/300)),1);
    for i = 1:floor(max(spxtimes/300))
      wamps = amps(spxtimes>((i-1)*300) & spxtimes<i*300);
      mean_amp_5min(i)   = mean(wamps);
      std_amp_5min(i)    = std(wamps);
      sem_amp_5min(i)    = sem(wamps);
      varerr_amp_5min(i) = var(wamps)/numel(wamps);  % SEM_diff = sqrt((sd2/na) + (sd2/nb))
      all_amps_5min{i}   = wamps;
      clear wamps
    end
    xval = 2.5:5:5*floor(max(spxtimes/300))-2.5;
    % gray lines plot the overall mean firing rate and its 95% confidence interval
    plot([xval(1),xval(end)],[mean(mean_amp_5min),mean(mean_amp_5min)],'Color',[.7 .7 .7])
    plot([xval(1),xval(end)],[mean(mean_amp_5min)+2*sem(mean_amp_5min),mean(mean_amp_5min)+2*sem(mean_amp_5min)],'Color',[.7 .7 .7])
    plot([xval(1),xval(end)],[mean(mean_amp_5min)-2*sem(mean_amp_5min),mean(mean_amp_5min)-2*sem(mean_amp_5min)],'Color',[.7 .7 .7])
    % red lines plot 5-min average firing rates and their binwise 95% confidence intervals
    plot(xval,mean_amp_5min,'r')
    plot(xval,mean_amp_5min+2*sem_amp_5min,'r')
    plot(xval,mean_amp_5min-2*sem_amp_5min,'r')
    axis([0 max(spxtimes)/60 nanmean(mean_amp)-5*nanstd(mean_amp) nanmean(mean_amp)+5*nanstd(mean_amp)])
    % we consider a unit waveform stable from one 5-min interval to the next if their mean amplitudes differ by less than 2 SE
    isstable = abs(diff(mean_amp_5min'))<2*sqrt(varerr_amp_5min(1:end-1)+varerr_amp_5min(2:end))';  % actually would need SE of the difference
    c.isstable = isstable;
    idx = find(isstable)+1;
    for i = 1:sum(isstable)
      plot([xval(idx(i)-1),xval(idx(i))],[nanmean(mean_amp)+4*nanstd(mean_amp),nanmean(mean_amp)+4*nanstd(mean_amp)],'k.-')
    end
    idx = find(~isstable)+1;
    for i = 1:numel(idx)
      plot([xval(idx(i)-1),xval(idx(i))],[nanmean(mean_amp)+3.5*nanstd(mean_amp),nanmean(mean_amp)+3.5*nanstd(mean_amp)],'m.-')
    end
    x = cell2vectors(all_amps_5min);
    [c.pANOVAamp,c.tableANOVAamp] = anova1(x(:,1),x(:,2),'off');
    c.eta2Amp = c.tableANOVAamp{2,2}/c.tableANOVAamp{4,2};
    set(get(gca,'XLabel'),'String','Time (min)','FontSize',10);
    text(5,nanmean(mean_amp)-4*nanstd(mean_amp),['eta^2=',num2str(c.eta2Amp,'%2.3f'),', p=',num2str(c.pANOVAamp,'%2.4f')],'Color','r')
    
    subplot(4,5,10),hold on
    g = nan(numel(isstable)+1);
    for i = 1:size(g,1)
      g(i,:) = (mean_amp_5min(i,1)-mean_amp_5min) ./ sqrt(nanmean([repmat(std_amp_5min(i)^2,size(g,1),1),std_amp_5min.^2],2));
    end
    g(logical(eye(size(g))))=nan;
    imagesc(xval,xval,g)
    colormap(gca,autumn)
    colorbar
    [row,col] = find(abs(g)>0.5);
    plot(xval(row),xval(col),'k.')
    set(get(gca,'XLabel'),'String','Time (min)','FontSize',10);
    axis([xval(1)-2.5 xval(end)+2.5 xval(1)-2.5 xval(end)+2.5])
    clear mean_amp mean_amp_5min sem_amp_5min std_amp_5min varerr* xval isstable x row col idx1 idx2
  else
    subplot(459),hold on
    axis([0 1 0 1])
    text(0,0.5,'No timestamps','FontName','Courier')
    axis off
    subplot(4,5,10),hold on
    axis([0 1 0 1])
    text(0,0.5,'No timestamps','FontName','Courier')
    axis off
  end
catch
  subplot(459),axis off
  axis([0 1 0 1])
  text(0.9,0.5,'stability analysis failed')
end
%% subplot 11 - interspike intervals, zoomed out
subplot(4,5,11)
if exist('spxtimes','var') && sum(spxtimes)>0
  title('ISI (10-ms bins)','FontSize',10),hold on
  resol = 0.5;
  [n,b]=hist(diff(spxtimes),0.005:0.01:(max(diff(spxtimes))));
  bar(b(1:min([resol*1000 numel(b)])),n(1:min([resol*1000 numel(b)])),'b','EdgeColor','b')
  axis([0 resol 0 ceil(max(n)*1.1)])
  set(get(gca,'XLabel'),'String','Interval (s)','FontSize',10);
  ylabel('Frequency')
  clear b n resol
else
  text(0,0.5,'No timestamps','FontName','Courier')
  axis off
end
%% subplot 12 - interspike intervals, zoomed in
subplot(4,5,12)
if exist('spxtimes','var') && sum(spxtimes)>0
  title('ISI (1-ms bins)','FontSize',10),hold on
  try
    resol = 0.05;
    [n,b]=hist(diff(spxtimes),0.0005:0.001:(max(diff(spxtimes))));
    bar(b(1:resol*1000),n(1:resol*1000),'b','EdgeColor','b')
    axis([0 resol 0 ceil(max(n)*1.1)])
    set(get(gca,'XLabel'),'String','Interval (s)','FontSize',10);
    c.violations = sum(n(1:nospx))/(c.spx-1); % -1 because we are looking at interspike intervals here
    plot([nospx/1000,nospx/1000],[0 ceil(max(n)*1.1)],'k:')
    clear b n resol
  end
else
  text(0,0.5,'No timestamps','FontName','Courier')
  axis off
end
%% subplot 13 - density plot
subplot(4,5,13),hold on
numTicks = size(wavs,2);
clear b n
for i = 1:numTicks % for each tick, check distribution of voltages
  [n(i,:),b(i,:)] = hist(wavs(:,i),linspace(lower,upper,min([numel(unique(wavs))/4 100])));
%   [n(i,:),b(i,:)] = hist(wavs(:,i),-0.1:0.005:0.1); % just if you need this standardized
end
colormap(gca,hot)
% remove extreme outliers
cutoff = 5*std(reshape(n,numel(n),1));
n(n>cutoff) = cutoff;
pcolor(n'),shading interp
axis([1 numTicks 1 size(n,2)])
c.nticks = size(wavs,2);
set(gca,'XTickLabel',[],'YTickLabel',[])

clear cutoff b n i x
%% subplot 14 - spike count distribution over the course of the session
if exist('spxtimes','var') && sum(spxtimes)>0
  subplot(4,5,14),title('Spikes per minute'),hold on
  n = histcounts(spxtimes,60:60:max([max(spxtimes) 1]));
  histogram(n,'FaceColor',[.5,.5,.5])

  % add a Poisson distribution
  ppdf = poisspdf(0:max(n),mean(n));
  [~,F] = mode(n);
  plot(0:max(n),F*ppdf/max(ppdf),'k-','LineWidth',2)
  plot([mean(n),mean(n)],[0,F],'k:','LineWidth',2)
  c.devskew = [skewness(n),(mean(n)-median(n))/std(n)];
  xlabel('N spikes'),ylabel('Frequency')

  clear n ppdf F
else
  text(0,0.5,'No timestamps','FontName','Courier')
  axis off
end
%% subplot 15 - spike shape
if fs
  waveParms = mwaveInMcheck(c.meanwave,10^6/fs);
  % waveParms = [fwhmMax,fwhmMin,p2tAmp,p2tDur,p2tRat]
  c.waveParms = waveParms;
  clear waveParms
end
%% subplots 16&17 - PSTH and eliminated waveforms (if any)
if exist('tevents','var')
  
  % low-pass filter events
  todelete = zeros(numel(tevents),1);     % becomes array to index events for deletion
  i=2;
  while i<=numel(tevents)                 % while loop marks events for deletion
    if tevents(i)-tevents(i-1)<minInt,todelete(i)=1;end,i=i+1;
  end
  tevents2 = tevents(todelete==0);
  
  % eliminate spikes in close proximity to specified events, find spikes close to it and note their index
  elimWin = elimWin/1000;
  dummy   = zeros(numel(spxtimes),1);
  for p = 1:numel(tevents2);
    try
      dummy(spxtimes>=tevents2(p)-elimWin(1) & spxtimes<=tevents2(p)+elimWin(2)) = dummy(spxtimes>=tevents2(p)-elimWin(1) & spxtimes<=tevents2(p)+elimWin(2)) + 1;
    end
  end
  tspikes2 = spxtimes;
  tspikes2(dummy>0) = [];
  
  plotrange = -0.05:0.001:+0.05;
  psth  = zeros(numel(plotrange),1);
  psth2 = zeros(numel(plotrange),1);
  
  % construct psth for original spike array
  for i = 1:numel(tevents2)
    stimes = spxtimes - tevents2(i);
    trialtimes = round(1000*(stimes(stimes>=plotrange(1) & stimes<=plotrange(end))));
    psth(plotrange(1)*-1000+trialtimes+1) = psth(plotrange(1)*-1000+trialtimes+1)+1;
  end
  
  % construct psth for clean spike array
  for i = 1:numel(tevents2)
    stimes = tspikes2 - tevents2(i);
    trialtimes = round(1000*(stimes(stimes>=plotrange(1) & stimes<=plotrange(end))));
    psth2(plotrange(1)*-1000+trialtimes+1) = psth2(plotrange(1)*-1000+trialtimes+1)+1;
  end
  
  subplot(4,5,16)
  title(['PSTH, elimWin=[' num2str(elimWin(1)) ' ' num2str(elimWin(2)) ']'],'FontSize',10),hold on
  plot(plotrange,psth,'m')
  plot(plotrange,psth2,'b')
  set(get(gca,'XLabel'),'String','Peri-event time (s)','FontSize',10);
  
  % waveforms eliminated due to closeness to trigger event
  subplot(4,5,17)
  title([num2str(sum(dummy>0)) ' waveforms in event vicinity'],'FontSize',10), hold on
  plot(wavs(dummy>0,:)')
  set(get(gca,'XLabel'),'String','Time (bins)','FontSize',10);
  axis([0 numTicks+1 lower upper])
  
else
  subplot(4,5,16)
  text(0,0.5,'No timestamps','FontName','Courier')
  axis off
  subplot(4,5,17)
  text(0,0.5,'No timestamps','FontName','Courier')
  axis off
  
end
clear dummy i numTicks p plotrange psth psth2 stimes tevents tevents2 todelete trialtimes tspikes2
%% subplot 18 - recommendations
subplot(4,5,18),axis off,axis([0 1.3 -0.1 1.1]),hold on
c.SNR_dmax = (mean(mwinsor(a1,1)) - mean(mwinsor(wavs(:,1),1))) / std(mwinsor(wavs(:,1),1));  % denominators for SNR used to be pooled variance: sqrt(0.5 * var(wavs(:,1)) + var(a1))
c.SNR_dmin = (mean(mwinsor(a3,1)) - mean(mwinsor(wavs(:,1),1))) / std(mwinsor(wavs(:,1),1));  % SD of the 1%-winsorized noise distribution
c.SNR_d    = c.SNR_dmax - c.SNR_dmin;

% another SNR measure, see Stratton et al. 2012, PLoS ONE: RMS_signal / RMS_noise
% RMS = sqrt(sum(x.^2)/numel(x))
% It is however not quite clear where the signal starts, so this is rather an underestimation because RMS_signal includes some noise data points.
c.SNR_rms  = mean(rms(squeeze(wavs(:,2:end)))) / rms(wavs(:,1));

text(0,1.2,[num2str(c.spx) ' spikes, ' num2str(c.fr,'%2.1f') ' Hz'],'FontName' ,'Courier')
text(0,1.05,'SNR(dmax) ','FontName','Courier'),text(0.7,1.05,num2str(c.SNR_dmax,'%2.1f'),'FontName' ,'Courier')
text(0,0.95,'SNR(dmin) ','FontName','Courier'),text(0.7,0.95,num2str(c.SNR_dmin,'%2.1f'),'FontName' ,'Courier')
text(0,0.85,'SNR(d)    ','FontName','Courier'),text(0.7,0.85,num2str(c.SNR_d,'%2.1f'),'FontName','Courier')
text(0,0.7,'SKEW(max) ','FontName','Courier'),text(0.7,0.7,num2str(max(c.skew(5:end-5)),'%2.1f'),'FontName','Courier')
text(0,0.6,'SKEW(min) ','FontName','Courier'),text(0.7,0.6,num2str(min(c.skew(5:end-5)),'%2.1f'),'FontName','Courier')
if isfield(c,'cvfr'),text(0,0.45,'CV(FR)     ','FontName','Courier'),text(0.7,0.45,num2str(c.cvfr,'%2.2f'),'FontName','Courier'),end
if isfield(c,'cvamps'),text(0,0.35,'CV(amps)   ','FontName','Courier'),text(0.7,0.35,num2str(c.cvamps,'%2.2f'),'FontName','Courier'),end
text(0,0.05,'FNneg','FontName','Courier'),text(0.7,0.05,num2str(c.FNmin,'%2.2f'),'FontName','Courier')
text(0,-0.05,'FNpos','FontName','Courier'),text(0.7,-0.05,num2str(c.FNmax,'%2.2f'),'FontName','Courier')

if c.SNR_dmax>crit_SNRd(1)
  rectangle('Position',[1.2 1.0 0.1 0.1],'FaceColor','g')
elseif c.SNR_dmax>crit_SNRd(2)
  rectangle('Position',[1.2 1.0 0.1 0.1],'FaceColor','y')
else
  rectangle('Position',[1.2 1.0 0.1 0.1],'FaceColor','r')
end
if c.SNR_dmin<-crit_SNRd(1)
  rectangle('Position',[1.2 0.9 0.1 0.1],'FaceColor','g')
elseif c.SNR_dmin<-crit_SNRd(2)
  rectangle('Position',[1.2 0.9 0.1 0.1],'FaceColor','y')
else
  rectangle('Position',[1.2 0.9 0.1 0.1],'FaceColor','r')
end
if c.SNR_d>crit_SNRdf(1)
  rectangle('Position',[1.2 0.8 0.1 0.1],'FaceColor','g')
elseif c.SNR_d>crit_SNRdf(2)
  rectangle('Position',[1.2 0.8 0.1 0.1],'FaceColor','y')
else
  rectangle('Position',[1.2 0.8 0.1 0.1],'FaceColor','r')
end
if max(c.skew(5:end-5))<crit_skew(1)
  rectangle('Position',[1.2 0.65 0.1 0.1],'FaceColor','g')
elseif max(c.skew(5:end-5))<crit_skew(2)
  rectangle('Position',[1.2 0.65 0.1 0.1],'FaceColor','y')
else
  rectangle('Position',[1.2 0.65 0.1 0.1],'FaceColor','r')
end
if min(c.skew(5:end-5))>-crit_skew(1)
  rectangle('Position',[1.2 0.55 0.1 0.1],'FaceColor','g')
elseif min(c.skew(5:end-5))>-crit_skew(2)
  rectangle('Position',[1.2 0.55 0.1 0.1],'FaceColor','y')
else
  rectangle('Position',[1.2 0.55 0.1 0.1],'FaceColor','r')
end
if exist('spxtimes','var') && c.cvfr<crit_CVfr(1)
  rectangle('Position',[1.2 0.4 0.1 0.1],'FaceColor','g')
elseif exist('spxtimes','var') && c.cvfr<crit_CVfr(2)
  rectangle('Position',[1.2 0.4 0.1 0.1],'FaceColor','y')
else
  rectangle('Position',[1.2 0.4 0.1 0.1],'FaceColor','r')
end
if c.cvamps<crit_CVa(1)
  rectangle('Position',[1.2 0.3 0.1 0.1],'FaceColor','g')
elseif c.cvamps<crit_CVa(2)
  rectangle('Position',[1.2 0.3 0.1 0.1],'FaceColor','y')
else
  rectangle('Position',[1.2 0.3 0.1 0.1],'FaceColor','r')
end
if isfield(c,'violations')
  text(0,0.25,['ISI <' (num2str(nospx)) ' ms'],'FontName','Courier'),text(0.7,0.25,num2str(c.violations,'%2.3f'),'FontName','Courier')
  if c.violations<crit_ISI(1)
    rectangle('Position',[1.2 0.2 0.1 0.1],'FaceColor','g')
  elseif c.violations<crit_ISI(2)
    rectangle('Position',[1.2 0.2 0.1 0.1],'FaceColor','y')
  else
    rectangle('Position',[1.2 0.2 0.1 0.1],'FaceColor','r')
  end
else
  rectangle('Position',[1.2 0.2 0.1 0.1],'FaceColor','k')
end
if c.FNmin<crit_FN(1)
  rectangle('Position',[1.2 0 0.1 0.1],'FaceColor','g')
elseif c.FNmin<crit_FN(2)
  rectangle('Position',[1.2 0 0.1 0.1],'FaceColor','y')
else
  rectangle('Position',[1.2 0 0.1 0.1],'FaceColor','r')
end
if c.FNmax<crit_FN(1)
  rectangle('Position',[1.2 -0.1 0.1 0.1],'FaceColor','g')
elseif c.FNmax<crit_FN(2)
  rectangle('Position',[1.2 -0.1 0.1 0.1],'FaceColor','y')
else
  rectangle('Position',[1.2 -0.1 0.1 0.1],'FaceColor','r')
end

axis([0 1.3 -0.1 1.1])
clear a1 a3
%% subplot 19:20 - sample waveforms
subplot(4,3,12)
[numWavs,numTicks] = size(wavs);
toPlot = floor(min([numWavs/3 200]));
if numWavs>200
  plot(1:numTicks,wavs(1:toPlot,:)','b'),hold on
  plot(numTicks+1:numTicks*2,wavs(round(numWavs/2-toPlot/2):round(numWavs/2+toPlot/2),:),'c')
  plot(numTicks*2+1:numTicks*3,wavs(end-ceil(toPlot/2):end,:)','r')
  axis([0 numTicks*3+1 lower upper])
  mult = 16:-4:-16;
  x = mult.*ones(1,numel(mult))*c.stdev(1)+c.meanwave(1);
  for i = 1:numel(x)
    plot([1,numTicks*3+1],[c.meanwave(1)+x(i) c.meanwave(1)+x(i)],'k:')
  end
end
axis off
ypos = double((max(c.meanwave)+2*c.stdev(1)));
text(1,ypos,['first ' num2str(toPlot) ' spikes'],'FontWeight','bold')
text(1+numTicks,ypos,['middle ' num2str(toPlot) ' spikes'],'FontWeight','bold')
text(1+numTicks*2,ypos,['last ' num2str(toPlot) ' spikes'],'FontWeight','bold')
clear i lower upper mult numWavs toPlot ypos
%% close plot if requested
if plotit==0
  close
end
end
%% nested function mwaveInMcheck
function waveParms = mwaveInMcheck(meanwave,si)
% waveParms = mwaveInMcheck(meanwave,si,varargin)
% modified from mwave in September 2016
%
% waveParms = [fwhmMax,fwhmMin,p2tAmp,p2tDur,p2tRat]
%
% ARGUMENTS
% meanwave    required - vector containing the mean spike waveform
% si          required - scalar value indicating sampling interval in ?s
% upsample waveform, get p2tAmp and p2tRat
if ~isa(meanwave,'double')
  meanwave = double(meanwave);
end
meanwave = meanwave-meanwave(1);       % pull start of meanwave to 0
upsamp   = 10;
dummywave    = interp1(meanwave,linspace(1,length(meanwave),length(meanwave)*upsamp),'spline'); % changed to spline from linear on Oct 19 2020
p2tAmp   = max(dummywave)-min(dummywave);      % peak-to-trough amplitude
p2tRat   = abs(max(dummywave)/min(dummywave)); % peak-to-trough ratio
% compute full width at half maximum and peak-to-trough duration
% extract FWHM from upsampled waveform in dummy
[posPeak,posPeakIdx] = max(dummywave);
[negPeak,negPeakIdx] = min(dummywave);

prePosPeakIdx  = find(dummywave(1:posPeakIdx)<posPeak/2,1,'last');
postPosPeakIdx = find(dummywave(posPeakIdx:end)<posPeak/2,1,'first') + posPeakIdx-1;
fwhmMax        = si*abs(postPosPeakIdx - prePosPeakIdx)/upsamp;

preNegPeakIdx  = find(dummywave(1:negPeakIdx)>negPeak/2,1,'last');
postNegPeakIdx = find(dummywave(negPeakIdx:end)>negPeak/2,1,'first') + negPeakIdx-1;
fwhmMin        = si*abs(postNegPeakIdx - preNegPeakIdx)/upsamp;

p2tDur         = si*abs(posPeakIdx - negPeakIdx)/upsamp;

if isempty(fwhmMax),fwhmMax=nan;end
if isempty(fwhmMin),fwhmMin=nan;end
if isempty(p2tAmp),p2tAmp=nan;end
if isempty(p2tDur),p2tDur=nan;end
if isempty(p2tRat),p2tRat=nan;end

waveParms = [fwhmMax,fwhmMin,p2tAmp,p2tDur,p2tRat];

subplot(4,5,15)
plot(dummywave,'Color','b','LineWidth',2),hold on
try
  plot([prePosPeakIdx prePosPeakIdx],[-10^3 10^3],'g:')
  plot([postPosPeakIdx postPosPeakIdx],[-10^3 10^3],'g:')
  plot([0 numel(dummywave)],[posPeak/2 posPeak/2],'g:')
  plot([prePosPeakIdx postPosPeakIdx],[posPeak/2 posPeak/2],'g')
  plot([preNegPeakIdx preNegPeakIdx],[-10^3 10^3],'c:')
  plot([postNegPeakIdx postNegPeakIdx],[-10^3 10^3],'c:')
  plot([0 numel(dummywave)],[negPeak/2 negPeak/2],'c:')
  plot([preNegPeakIdx postNegPeakIdx],[negPeak/2 negPeak/2],'c')
end
axis([0 numel(dummywave)+1 -(max(abs([max(dummywave) min(dummywave)]))) (max(abs([max(dummywave) min(dummywave)])))])
% scale x-axis from upsampled ticks to ms
set(gca,'XTick',0:upsamp*1000/si:numel(dummywave),'XTickLabel',num2str((0:1:floor(numel(dummywave)/(upsamp/si)/1000))'))
set(get(gca,'XLabel'),'String','Time (ms)','FontSize',10);
ylabel('ADC units')
end
%% nested function sem
function y = sem(x,varargin)
% function y = sem(x,varargin)
% function simply returns standard error of the mean
%
% INPUTS:
% x           vector or 2D-matrix containing single or double precision numbers
% varargin    if x is a matrix:
%             scalar specifying the dimension over which the SEM should be computed
%
% taken from mstats on Sep 15 2016
% input check
if numel(size(x))>2
  error('sorry, code only applicable for vectors and 2D matrices')
end
% the works
if isvector(x)
  y = nanstd(x)/sqrt(sum(~isnan(x)));
else
  if nargin<2
    y = nanstd(x,0,1)./sqrt(sum(~isnan(x)));
  else
    y = nanstd(x,0,varargin{1})./sqrt(sum(~isnan(x),varargin{1}));
  end
end
end
%% nested function cell2vectors
function y = cell2vectors(x)
% accepts a one-dimensional cell array containing numbers, concatenates all numbers and
% constructs a second vector indicating the identity of each numbers cell
% function comes in handy when using e.g. anova or kruskalwallis with unequal sample sizes
% per group
% Maik C. St?ttgen, February 2014
% check whether x is one-dimensional
if numel(x)~=max(size(x))
  error('input cell array is not one-dimensional')
end
% determine total number of elements across all cells and check whether they are all numeric
n = 0;
for i = 1:numel(x)
  if ~isnumeric(x{i})
    error(['cell no. ' num2str(i) ' is not numeric'])
  end
  n = n+numel(x{i});
end
% extract values
y = nan(n,2);
c = 1;
for i = 1:numel(x)
  y(c:c+numel(x{i})-1,:) = [x{i} i*ones(numel(x{i}),1)];
  c = c + numel(x{i});
end
end
%% nested function mwinsor
function w = mwinsor(X,percent)
% function w = mwinsor(X,percent) winsorizes a distribution
% 
% Winsorizing involves finding the highest and the lowest values (percent/2) and replacing each of them by the highest
% and lowest untrimmed scores, respectively.
% 
% INPUT
% X         a vector or a matrix which should be winsorized (along columns, by default)
% percent   scalar between 0 and 100; by default, mwinsor rounds to the nearest integer when k (the number of data points 
%           above or below a certain cutoff) is not an integer.

% check input
if ndims(X)>2 %#ok<*ISMAT>
  error('Only matrices with up to two dimensions are supported.')
end
if percent<0 || percent>100
  error('percent should lie between 0 and 100')
end

% generate winsorized vector/matrix
% first, set lowest k data values to nan, where k=n*(percent/100)/2, and where n is the number of elements
k = round(size(X,1)*(percent/100)/2);
w = sort(X);
w(1:k,:) = nan;
minima = min(w);
for i = 1:size(w,2)
  w(1:k,i) = minima(i);
end
w(end-k+1:end,:) = nan;
maxima = max(w);
for i = 1:size(w,2)
  w(end-k+1:end,i) = maxima(i);
end
end
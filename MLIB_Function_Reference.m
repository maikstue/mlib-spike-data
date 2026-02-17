%% MLIB FUNCTION REFERENCE
% This document describes the use of the functions contained in MLIB. Additional information for each function is provided through each function's 
% documentation. The functions are explained with reference to several example data files containing extracellular unit recordings from the brains of rats, 
% mice, and pigeons, taken from the following publications:
% 
% Kasties N, Starosta S, Güntürkün O, Stüttgen MC (2016)
% Neurons in the pigeon caudolateral nidopallium differentiate Pavlovian conditioned stimuli but not their associated reward value.
% Scientific Reports 6: 3549.
% 
% Lengersdorf D, Pusch R, Güntürkün O, Stüttgen MC (2014)
% Neurons in the pigeon nidopallium caudolaterale signal the selection and execution of perceptual decisions.
% European Journal of Neuroscience 40: 3316.
% 
% Stoilova VV, Knauer B, Berg S, Rieber E, Jäkel F, Stüttgen MC (2020)
% Auditory cortex reflects goal-directed movement but is not necessary for behavioral adaptation in sound-cued reward tracking.
% Journal of Neurophysiology 124: 1056.
% 
% van der Bourg A, Yang JW, Reyes-Puerta V, Laurenczy B, Wieckhorst M, Stüttgen MC, Luhmann HJ, Helmchen F (2017)
% Layer-specific refinement of sensory coding in developing mouse barrel cortex.
% Cerebral Cortex 27: 4835.
% 
% Vandevelde JR, Yang JW, Albrecht S, Lam H, Kaufmann P, Luhmann HJ, Stüttgen MC (2023)
% Layer- and cell-type-specific differences in neural activity in mouse barrel cortex during a whisker detection task.
% Cerebral Cortex 33: 1361.
% 
% Yeganeh F, Knauer B, Guimaraes Backhaus R, Yang JW, Stroh A, Luhmann HJ, Stüttgen MC (2022)
% Effects of optogenetic inhibition of a small fraction of parvalbumin-positive interneurons on the representation of sensory stimuli in mouse barrel cortex.
% Scientific Reports 12: 19419.
% 
%% MLIB EXAMPLE FILES
% There are several example files to illustrate the use of the MLIB functions for extracellular spike recordings.
% 
% File #1: unitForMLIBTesting_1_rat.mat
% File #2: unitForMLIBTesting_2_pigeon.mat
% File #3: unitForMLIBTesting_3_pigeon.mat
% File #4: unitForMLIBTesting_4_rat.mat
% File #5: unitForMLIBTesting_5_rat.mat
% File #6: unitForMLIBTesting_6_rat.mat
% File #7: data4mcluster.mat
% File #8: data4mcsd.mat
% File #9: data4mdimreduct.mat
% File #10: data4mep.mat
% File #11: data4mspectan.mat
% 
% The files' contents are described in detail in the accompanying documentation document.
% 
% NO.	FILE NAME	                  |  VARIABLES 	                    |  FUNCTION(S))	  | FIGURE PANEL(S) IN PAPER
% --------------------------------|---------------------------------|-----------------|-------------------------
% 1	  unitForMLIBTesting_1_rat	  |  spike times (tspx)             |  mnspx          |  none
%                                 |  event codes (cevents)          |  mcvar          |
%                                 |  event times (tevents)	        |                 |
% --------------------------------|---------------------------------|-----------------|-------------------------
% 2	  unitForMLIBTesting_2_pigeon	|  spike times (tspx)             |  mcheck         |  Figure 2a-h
%                                 |  event codes (cevents)          |  mpsth          |
%                                 |  event times (tevents)          |  msdf           |
%                                 |  spike waveforms (adc)          |  mwave          |
% --------------------------------|---------------------------------|-----------------|-------------------------
% 3	  unitForMLIBTesting_3_pigeon	|  spike timestamps from two      |  mcc            |  Figure 2g-m
%                                 |  different channels             |  mcheck         |
%                                 |  tspx_ch9 & _ch10)              |                 |
%                                 |  event codes (cevents)          |                 |
%                                 |  event times (tevents)          |                 |
%                                 |  waveforms from each channel    |                 |
%                                 |  (adc_ch9 & _ch10)              |                 |
%                                 |  cluster assignments for both   |                 |
%                                 |  channels| waveforms            |                 |
%                                 |  cluster_ch9 & _ch10)           |                 |
%                                 |  waveform sampling rate (fs)	  |                 |
% --------------------------------|---------------------------------|-----------------|-------------------------
% 4	  unitForMLIBTesting_4_rat	  |  spike times (tspx)             |  mlat           |  Figure 3a
%                                 |  event codes (cevents)          |  mmi            |
%                                 |  event times (tevents)          |  mwa            |
% --------------------------------|---------------------------------|-----------------|-------------------------
% 5	  unitForMLIBTesting_5_rat	  |  spike times (tspx)             |  mwa2           |  none
%                                 |  event codes (cevents)          |                 |
%                                 |  event times (tevents)	        |                 |
% --------------------------------|---------------------------------|-----------------|-------------------------
% 6	  unitForMLIBTesting_6_rat	  |  spike times (tspx)             |  mmultregress   |  Figure 3c
%                                 |  stimulus times (tstims)        |                 |
%                                 |  response times (tresponses)    |                 |
%                                 |  outcome times (toutcomes)      |                 |
%                                 |  stimulus type (stims)          |                 |
%                                 |  response type (responses)      |                 |
%                                 |  outcome type (outcomes)        |                 |
% --------------------------------|---------------------------------|-----------------|-------------------------
% 7	  data4mcluster	              |  concatenated peri-stimulus     |  mcluster       |  none
%                                 |  time histograms for            |                 |
%                                 |  98 single neurons              |                 |
% --------------------------------|---------------------------------|-----------------|-------------------------
% 8	  data4mcsd	                  |  average LFP responses from     |  mcsd           |  Figure 5
%                                 |  20 electrode channels recorded |                 |
%                                 |  over 80 ms                     |                 |
% --------------------------------|---------------------------------|-----------------|-------------------------
% 9	  data4mdimreduct	            |  spike rates of 300 units       |  mdimreduct	    |  Figure 4
%                                 |  relative to presentation of    |                 |
%                                 |  either of two stimuli          |                 |
% --------------------------------|---------------------------------|-----------------|-------------------------
% 10   data4mep	                  |  continuous field potential     |  mep            |  none
%                                 |  trace from a single electrode  |                 |
% --------------------------------|---------------------------------|-----------------|-------------------------
% 11   data4mspectan	            |  continuous EEG recording       |  mspectan	      |  none
%                                 |                                 |                 |
% 
%% MLIB FUNCTIONS IN ALPHABETICAL ORDER
% FUNCTION        BRIEF DESCRIPTION
% mcc             for spike data; visualizes and quantifies cluster isolation quality to aid spike sorting
% mcheck          for spike data; visualizes and quantifies single unit isolation quality to aid spike sorting
% mcluster        for spike data; clusters spike response patterns of several units
% mcsd            for LFP data; performs and visualizes current source density analysis
% mcvar           for spike data; computes coefficient of variation for interspike intervals and Fano factor for spike counts
% mdimreduct      for spike data; performs dimensionality reduction to visualize the population spike response in 2-3 principal components
% mep             for LFP data; computes and plots an evoked potential with error bars and confidence intervals
% mlat            for spike data; computes different latency measures for a single neuron
% mmi             for spike data; computes mutual information between two different spike trains
% mmi_psych       for behavioral data; expresses discrimination performance in a two-choice task as mutual information (for comparison with spike MI)
% mmultregress    for spike data; performs multiple regression analysis to dissect the influence of different factors on the spike response for a single neuron
% mnspx           for spike data; counts the number of spikes in a time window relative to a series of events (e.g., stimulus presentations)
% mpsth           for spike data; computes and plots a peri-stimulus time histogram (PSTH)
% mroc            for spike data; computes the area under the receiver-operating characteristic (ROC) curve for two distributions of spike counts
% msdf            for spike data; constructs a spike-density function from a PSTH using different kernels
% mspectan        for LFP data; performs a spectral analysis and plots the result
% mtune           for spike data; computes 1-D tuning curve and fits a normal curve
% mvecstrength    for spike data; computes vector strength of a spike train relative to a cyclic event
% mwa             for spike data; computes the area under the ROC curve (AUROC) in moving windows of different sizes relative to a fixed reference event
% mwa2            for spike data; similar to mwa, but computes AUROC for two events
% mwave           for spike data; computes measures to characterize the average spike waveform of a single neuron
    
%% mcc.m
% [L_ratio,ID,A,J3,PsF,DBVI] = mcc(waveforms,cid,whichclusters,ts,whichPCs,nocorrs)
% function to assess the quality of spike sorting (post hoc)
% computes standard quality measures (L_ratio, J3,...); see doc mcc for more information
% This function requires the Statistics and Machine Learning Toolbox.
% 
% Recommended reading:
% 1) Schmitzer-Torbert, N., Jackson, J., Henze, D., Harris, K., & Redish, A. D. (2005). Quantitative measures of cluster quality for 
%    use in extracellular recordings. Neuroscience, 131(1), 1–11.
% 2) Nicolelis, M. A. L., Dimitrov, D., Carmena, J. M., Crist, R., Lehew, G., Kralik, J. D., & Wise, S. P. (2003). Chronic, multisite, 
%    multielectrode recordings in macaque monkeys. PNAS 100(19), 11041–11046.

% first, clean up and load example data
clear all,close all,clc
load unitForMLIBTesting_3_pigeon

% This file contains waveforms from two electrode channels. Let's first look at waveforms from channel 9.

% an electrode with two clusters, take the first three principal components
% 'adc' contains all waveforms detected on the respective electrode channel
% 'clusters' contains the results of spike sorting - trash (cluster 0) or putative spike waveforms

whichclusters = [1,2];    % we will look at clusters 1 and 2 (and ignore 0)
whichPCs = 1:3;           % we will use the first three principal components for calculation of quality indices

% Example 1
% run the code; for explanation of the two figures and the subplots, refer to the MLIB documentation
[L_ratio,ID,A,J3,PsF,DBVI] = mcc(adc_ch9,clusters_ch9,whichclusters,tspx_ch9,whichPCs);

% Example 2
% now let's have a look at the waveforms from channel 10; here, three putative spike clusters were isolated
whichclusters = [1,2,3];
[L_ratio,ID,A,J3,PsF,DBVI] = mcc(adc_ch10,clusters_ch10,whichclusters,tspx_ch10,whichPCs,'nocorrs');    % with 'nocorrs', we omit the second figure showing crosscorrelograms
 
return
 
%% mcheck.m
% c = mcheck(waveforms,varargin)
% function for assessing spike sorting quality
% returns a range of quality matrices, see doc mcheck for explanation
% This function requires the Statistics and Machine Learning Toolbox as well as the Curve Fitting Toolbox.
% 
% Recommended reading:
% 1) Starosta, S., Stüttgen, M. C., & Güntürkün, O. (2014). Recording single neurons’ action potentials in freely moving pigeons 
%    across three stages of learning. Journal of Visualized Experiments, 88, e51283.
% 2) Stoilova, V. v, Knauer, B., Berg, S., Rieber, E., Jäkel, F., & Stüttgen, M. C. (2020). Auditory cortex reflects goal-directed movement 
%    but is not necessary for behavioral adaptation in sound-cued reward tracking. Journal of Neurophysiology, 124(4), 1056–1071.

clear all,close all,clc

load unitForMLIBTesting_2_pigeon.mat

% input all waveforms (sorted spikes and noise)
c = mcheck(adc);

% input waveforms and waveform timestamps
c = mcheck(adc,'timevec',tspx);

% all of the above, plus event timestamps
c = mcheck(adc,'timevec',tspx,'tevents',tevents(cevents==13));

% now try the other file
clear all
load unitForMLIBTesting_3_pigeon.mat

% input all waveforms from channel 9; the result looks messy
c = mcheck(adc_ch9);

% now input only waveforms which have been classfied into cluster 1 during spike sorting
c = mcheck(adc_ch9,'markvec',clusters_ch9,'target',1);

% as above, but additionally provide spike timestamps
c = mcheck(adc_ch9,'markvec',clusters_ch9,'target',1,'timevec',tspx_ch9);

% as above, but also provide timestamps for a PSTH (relative to pecks on the center key) and sampling rate
c = mcheck(adc_ch9,'markvec',clusters_ch9,'target',1,'timevec',tspx_ch9,'tevents',tevents(cevents==109),'fs',fs);

% same, but for cluster 2
c = mcheck(adc_ch9,'markvec',clusters_ch9,'target',2,'timevec',tspx_ch9,'tevents',tevents(cevents==109),'fs',fs);

% Note that some subplots are only shown when the appropriate inputs to the function are provided.
% subplot(451)    shows up to 5,000 waveforms superimposed; ordinate plot ADC units (usually, mV or µV)
% subplot(452)    shows the SD and the skewness of the voltage distributions as a function of time (waveform ticks)
%                 solid lines denote mean SD (blue) and mean skewness (red); dashed red lines denote skewness of +/- 0.5
%                 ideally, skewness should be 0, but variation within these boundaries is acceptable for waveform variability
%                 standard deviation is expressed as percentage of the mean noise SD (SD in the first five ticks); the lower its variability, the better
% subplot(453)    shows the frequency distributions of the noise amplitude (first bin of all waveforms, blue), the minimum amplitudes (green), 
%                 and the positive amplitudes (red)
%                 generally, all distributions should be approximately normal and sufficiently separated from each other
% subplots(454)   and (455) zoom in on the distributions of minima and maxima shown in subplot(453)
%                 the distributions are fit with Gaussians to estimate the fraction of false negatives (missed spikes) due to a spike detection 
%                 threshold which may have been to conservative; if there is a sizable fraction of false negatives, the distribution will appear 
%                 truncated at the voltage value of the spike detection threshold; the resulting estimates of false negatives (FNs) are returned in
%                 subplot(4,5,18) as FNneg (FNs estimated from the minima distribution, negative voltage threshold) and FNpos (FNs estimated 
%                 from the maxima distribution, positive voltage threshold)
% 
% 
% subplot(456)    blue: cumulative spike count in 1-second bins as a function of experiment time
%                      black line (main diagonal) serves as visual aid to assess constancy of firing
%                 red:  average peak-to-peak amplitude of the spike waveform in 1-minute bins as a function of time, +/- 1 SD
% subplot(457)    shows the minute-wise spike counts of the unit as a function of experiment time (black bars), together with 
%                 all spike amplitudes (red dots); the black line is a regression of time vs. spike amplitude;
%                 if negative, the amplitude tends to decrease over the course of the experiment
% subplot(458)    shows a scatterplot of the firing rate (x) vs. the mean amplitude (y), calculated over 1-min bins;
%                 in theory, the correlation should be 0; a positive correlation may suggest that units are gradually lost over time
%                 (both amplitude and spike rate decrease)
% subplot(459)    shows the mean spike amplitude with 95% confidence intervals (red lines), calculated in consecutive 5-minute bins
%                 gray horizontal lines plot the overall mean firing rate and its 95% confidence interval
%                 eta^2 and the p-value come from a oneway ANOVA; a significant p-value indicates that the amplitude changes over time
%                 eta^2 is a measure of effect size indicating what fraction of amplitude variability is explained by the factor time 
%                 (ideally, this would be 0, but values <0.10 are acceptable)
%                 if the waveform distribution of one 5-min segment differs significantly from the previous one, a black dot is plotted at the top
% subplot(4,5,10) shows the color-coded standardized difference (Hedges' g) of amplitude distributions for all pairs of 5-min segments
%                 if the |g|>0.5, the cell is marked with a black dot, indicating substantially different average waveform distributions
% 
% 
% subplot(4,5,11) plots the interspike-interval distribution from 0 to 500 ms in 10-ms bins (if timevec is provided)
% subplot(4,5,12) plots the interspike-interval distribution from 0 to 50 ms in 1-ms bins  (if timevec is provided)
% subplot(4,5,13) is a density plot of the spike waveforms; the frequency of specific time-voltage combinations is color-coded;
%                 the plot allows to assess the 'purity' of the waveform distribution better than a simple waveform overlay as in subplot(451)
% subplot(4,5,14) shows the distribution of spike counts for 1-minute bins; if firing rate is constant, one would expect to obtain a 
%                 Poisson distribution (which is superimposed on the histogram); vertical dotted line denotes mean 1-min spike count
% subplot(4,5,15) shows the average waveform and highlights some waveform features extracted by mwave.m (if fs is provided)
% 
% 
% subplot(4,5,16) shows a PSTH (+/-50 ms) for an arbitrary event whose time stamps is specified with 'tevents'
%                 blue: all spikes, magenta: spikes within elimWin; the latter waveforms are superimposed in subplot(4,4,14)
%                 this plot is useful e.g. to check for potential artifacts close to event occurrence (in my case, waveforms resembling 
%                 spikes induced by key pecks of the bird I recorded from)
% subplot(4,5,17) plots the waveforms in the immediate vicinity (+/-20 ms) of the event specified in tevents to assess whether
%                 the event induced any artifacts resembling spike waveforms
% subplot(4,5,18) gives some basic information and indices on unit quality
%                 SNR(dmax) is SDs from mean of the noise band to the mean of maxima distribution; likewise for SNR(dmin)
%                           SNR is measured as the distance between mean of the noise and the mean signal amplitudes as shown in subplot(443), 
%                           divided by the standard deviation of the noise; the SD is 1%-winsorized (see mwinsor.m)
%                 SNR(d)    is simply the difference between the SNR(dmax) and SNR(dmin)
%                 SKEW      gives the maximum and minimum skewness values for the central portion of the waveforms (tick 5 to tick end-5)
%                           we'd like abs(skew)<0.3, high skewness goes along with a multi unit (i.e., spikes from several units clustered together)
%                 CV(FR)    is the coefficient of variation of the firing rates (calculated in 1-min bins); by experience, a very stable unit has CV(FR)<0.4
%                           unstable units (with decreasing amplitude) frequently have values >0.8
%                 CV(amps)  as CV(FR), but calculated for spike amplitudes; by experience, a very stable unit has CV(amps)<0.1
%                 ISI <4 ms provides the fraction of waveforms with interspike intervals less than 4 ms, we'd like this to be below 2%
%                           an ISI <4 ms indicates violation of the refractory period --> the two spikes cannot be from the same neuron
%                           note that the value of the critical ISI can be specified in the function (variable 'nospx'); a reasonable value
%                           depends on the brain region of interest
%                 FNneg/FNpos estimate the fraction of false negatives (missed spikes) if the amplitude threshold used for spike detection 
%                           was in the negative range (FNneg) and the positive range (FNpos), respectively; 
%                           the colored traffic lights (green, yellow, red) indicate whether I would consider this value to be optimal (green, candidate 
%                           for high-quality single unit), worth analyzing at all (yellow, potentially multi unit), or trash (red)
% subplot(4,5,20) plots the first, middle and last 200 spikes in the file (if timevec is provided); useful to gauge waveform stability over time

return

%% mcluster.m
% [cid,avsil] = mcluster(data,z,nclus,whichPCs,cType,whichDist,metric,T,plotting)
% performs clustering on a data matrix after running PCA and returns cluster codes cid and average silhouette values
% input 'data' is a matrix where rows are units and columns are binned firing rates or eta2/omega2 values
% This function requires the Statistics and Machine Learning Toolbox.
% 
% Recommended reading:
% Kasties, N., Starosta, S., Güntürkün, O., & Stüttgen, M. C. (2016). Neurons in the pigeon caudolateral nidopallium differentiate 
% Pavlovian conditioned stimuli but not their associated reward value in a sign-tracking paradigm. Scientific Reports, 6 (October), 35469.

clear all,close all,clc

load('data4mcluster.mat')

% The data file contains PSTHs for 98 units; for each unit, we have five PSTHs with 25 data points each. These five PSTHs were appended for clustering.
% First, let us have a look a random selection of 5 units; for each unit, we have 5 PSTHs.
exampleUnits = sort(randsample(1:98,5,'false'));      % randomly select five units
xticks = 0.1:0.2:4.9;   % sampling times (0.2 s per sample)
figure('name','Some PSTHs')
for i = 1:5
  for j = 1:5
    subplot(5,5,(i-1)*5+j),title(['Unit ',num2str(exampleUnits(i)),', stim ',num2str(j)]),hold on
    plot(xticks,zpsths(exampleUnits(i),(j-1)*25+1:25*j))
    if i==5 && j==1,xlabel('Time (s)'),ylabel('z(spikes per s)'),end
    ylim([-2,3.5])
  end
end

% It is apparent that the units respond differently to the five stims (in fact, many units prefer one of the stimuli).
% Let's see whether we can cluster the 5 units into groups, depending on their responses to ALL five stimuli.

nclus = 5;              % choose into how many clusters the data should be partitioned
whichPCs = 1:3;         % mcluster will only use these PCs for clustering
cType = 'hieragglom';   % clustering method
whichDist = 'average';  % how to measure the distance between clusters
metric = 'correlation'; % distance metric

% input matrix contains the raw PSTHs which are zscored by mcluster
[cid,avsil] = mcluster(app_psths,'zscore',nclus,whichPCs,cType,'average','correlation',[],'yes');

% now plot the z-scored PSTHs for each cluster
figure('units','normalized','Position',[.4,.55,.5,.2])
for i = 1:nclus
  subplot(1,nclus,i),title(['Cluster ',num2str(i)]),hold on
  for j = 1:5
    plot(0.1:0.2:4.9,mean(zpsths(cid==i,(j-1)*25+1:j*25)))
  end
  axis([0,5,-1.5,3])
  if i==1,xlabel('Time (s)'),ylabel('z(FR)'),end
end
legend('Stim 1','Stim 2','Stim 3','Stim 4','Stim 5','units','normalized','Position',[0.45,0,1,1])

% With 5 clusters, one can see that in four of the clusters, units prefer one of the five stimuli (color-coded) in their steady-state activity.
% In one of the clusters, units fire briefly upon stimulus onset.

return

%% mcsd.m
% [csd_r,csd_i,lfp_s] = mcsd(lfpmat,es,doplot,t,theMinus)
% performs a current source density analysis for laminar LFP recordings
% 
% Recommended reading:
% 1) Vandevelde, J. R., Yang, J.-W., Albrecht, S., Lam, H., Kaufmann, P., Luhmann, H. J., & Stüttgen, M. C. (2023). Layer- and cell-type-specific 
% differences in neural activity in mouse barrel cortex during a whisker detection task. Cerebral Cortex, 33, 1361–1382.
% 2) Yeganeh, F., Knauer, B., Backhaus, R. G., Yang, J. W., Stroh, A., Luhmann, H. J., & Stüttgen, M. C. (2022). Effects of optogenetic inhibition 
% of a small fraction of parvalbumin‑positive interneurons on the representation of sensory stimuli in mouse barrel cortex. Scientific Reports, 12, 19419.
% 3) van der Bourg, A., Yang, J.-W., Reyes-Puerta, V., Laurenczy, B., Wieckhorst, M., Stüttgen, M. C., Luhmann, H. J., & Helmchen, F. (2017). Layer-specific 
% refinement of sensory coding in developing mouse barrel cortex. Cerebral Cortex, 27, 4835–4850.

clear all,close all,clc

load('data4mcsd.mat')     % Data kindly provided by Dr Jenq-Wei Yang and Prof. Heiko Luhmann, University Medical Center Mainz, Germany.

% file contains LFPs recorded with a linear silicon probe (20 electrodes spaced 50 µm apart)
% recordings were conducted in mouse barrel cortex sampling rate 1000 Hz, unit mV
% at time point 0, a whisker stimulus was applied
% recordings range from -8 to +71 ms relative to stimulus onset
% mean_FP contains the LFP averaged over 20 trials
% also see van der Bourg et al. (2017), Figure 3

electrode_spacing = 0.05;     % inter-electrode distance in mm
doplot = 1;                   % plot the data (1) or not (0)
time_vec = -8:71;             % in ms
theMinus = 1;                 % multiply the 2nd spatial derivative with -1 or not (1/0, resp.); default 1

[csd_r,csd_i,lfp_s] = mcsd(mean_fp,electrode_spacing,doplot,time_vec,theMinus);

% subplot(221) shows the raw LFP traces, color-coded for depth from black (top) to white (bottom)
% subplot(222) shows the LFP traces both as black lines and color-coded (after smoothing)
% subplot(223) shows the raw CSD
% subplot(224) shows the CSD after 10x interpolation and spatial smoothing

% Since the CSD is the spatial derivative of the LFP, there two rows less in csd_r than in mean_FP.
set(gcf,'units','normalized','Position',[.2,.2,.6,.6])

return

%% mcvar.m
% cv = mcvar(x,mode)
% computes the coefficient of variation for interspike intervals and the Fano Factor for a spike count distribution
% 
% Recommended reading:
% Maimon, G., & Assad, J. A. (2009). Beyond Poisson: increased spike-time regularity across primate parietal cortex. Neuron, 62(3), 426–440.

clear all,close all,clc

load('unitForMLIBTesting_1_rat.mat')

cv_isi = mcvar(tspx,'diff')    % field tspx holds spike times in seconds relative to the beginning of the recording

spike_counts = histcounts(tspx,0:999); % get spike counts for consecutive 1-s bins over the first 1000 seconds of recording
cv_ff = mcvar(tspx,'abs')

return

%% mdimreduct.m
% [score,loadings] = mdimreduct(smat4dimred,standardize,plotresults);
% performs PCA-based dimensionality reduction of a sample of neurons recorded simultaneously or under identical conditions
% This function requires the Statistics and Machine Learning Toolbox.
% 
% Recommended reading:
% Ebitz, R. B., & Hayden, B. Y. (2021). The population doctrine in cognitive neuroscience. Neuron, 109(19), 3055–3068.

clear all,close all,clc

load data4mdimreduct

% nspx_s1 and nspx_s2 contain spike count data from 300 MFC units
% rows represent units, columns represent consecutive 100-ms time windows, shifted in steps of 50 ms
% elements represent average spike count over >100 trials

% first, let's look at the population trajectories separately for S1 and S2 trials
mdimreduct(nspx_s1,1,'plot');   % performs dimensionality reduction for all neurons relative to the onset of stimulus 1
mdimreduct(nspx_s2,1,'plot');   % same, but for stimulus 2

% now concatenate the two matrices so we can compare the trajectories for S1 and S2 trials
score = mdimreduct([nspx_s1,nspx_s2],1,'plot');

% in the last figure, mdimreduct plots them as one, but we'd like to plot them in a different way to facilitate comparison
% this needs a few lines of code...
nobs   = size(score,2)/2; % number of observations
colors = hot(nobs+8);

figure('name','Population trajectories for S1 and S2 trials','units','normalized','position',[.1,.3,.8,.5])
subplot(251)
for i = 1:2
  irange = (i-1)*nobs+1:i*nobs;
  plot3(score(1,irange),score(2,irange),score(3,irange),'Color',[.5,.5,.5]),hold on
end
for i = 1:nobs
  plot3(score(1,i),score(2,i),score(3,i),'.','Color',colors(i,:),'MarkerSize',12)
  plot3(score(1,i+nobs),score(2,i+nobs),score(3,i+nobs),'.','Color',colors(i,:),'MarkerSize',12)
end
title('PC 1 vs. PC 2 vs. PC 3')
xlabel('PC 1'),ylabel('PC 2'),zlabel('PC 3')
grid on

% top row - PCs plotted against each other
for i = 1:4
  subplot(2,5,i+1),hold on
  switch i
    case 1,dims = [1,2];xlabel('PC 1'),ylabel('PC 2'),title('PC 1 vs. PC 2')
    case 2,dims = [1,3];xlabel('PC 1'),ylabel('PC 3'),title('PC 1 vs. PC 3')
    case 3,dims = [2,3];xlabel('PC 2'),ylabel('PC 3'),title('PC 2 vs. PC 3')
    case 4,dims = [1,4];xlabel('PC 1'),ylabel('PC 4'),title('PC 1 vs. PC 4')
  end

  for k = 1:2
    irange = (k-1)*nobs+1:k*nobs;
    plot(score(dims(1),irange),score(dims(2),irange),'Color',[.5,.5,.5])
    for j = 1:nobs
      plot(score(dims(1),j),score(dims(2),j),'.','Color',colors(j,:),'MarkerSize',12);
      plot(score(dims(1),j+nobs),score(dims(2),j+nobs),'.','Color',colors(j,:),'MarkerSize',12);
    end
  end
  grid on
end

% bottom row - PCs plotted individually as a function of time
for i = 1:5
  subplot(2,5,i+5),title(['PC ',num2str(i)]),hold on
  for k = 1:2
    irange = (k-1)*nobs+1:k*nobs;
    plot(winpos,score(i,irange),'Color',[.5,.5,.5])
    for j = 1:nobs
      if k==1
        plot(winpos(j),score(i,j),'.','Color',colors(j,:),'MarkerSize',12)
      elseif k==2
        plot(winpos(j),score(i,j+nobs),'.','Color',colors(j,:),'MarkerSize',12)
      end
    end
  end

  plot([0,0],[-25,25],'k:')
  plot([0.12,0.12],[-25,25],'k:')
  axis([-0.45,0.45,-25,25])
  if i==1,xlabel('Time re: stimulus onset (s)'),end
end

% In the last figure, the vertical dotted lines mark stimulus onset (time = 0 s) and the time where the animal withdraws from the nose port,
% i.e., the beginning of the choice response (t=0.12 s). In PCs 1 & 2, the trajectories for the two trial types diverge only when the animal
% is already moving to the left or to the right to indicate its choice. In the other three PCs, the trajectories diverge during stimulus presentation
% already, possibly indicative of differential representation of the stimulus and/or the decision.

return

%% mep.m
% function [ep,errb,ci] = mep(fp,tevents,pre,post,varargin)
% computes and plots an evoked = event-related potential from continuous LFP data
% optionally returns tick-wise error bars (SEM) and confidence intervals
% % This function requires the Statistics and Machine Learning Toolbox.
% 
% Recommended reading:
% Yeganeh, F., Knauer, B., Backhaus, R. G., Yang, J. W., Stroh, A., Luhmann, H. J., & Stüttgen, M. C. (2022). Effects of optogenetic inhibition 
% of a small fraction of parvalbumin‑positive interneurons on the representation of sensory stimuli in mouse barrel cortex. Scientific Reports, 12, 19419.

clear all,close all,clc

load data4mep   % contains 42 minutes of continuous LFP recording, high-passed filtered at 250 Hz (4th order Butterworth)
                % recording contains light artifacts from laser stimulation 20 ms before stimulus onset

pre_ticks     = 200;
post_ticks    = 200;
tevents_ticks = round(tevents(cevents==5)*1000);   % tevents were in seconds; since sampling rate was 1000 Hz, we convert to ticks by multiplication
                                                   % cevents=5 is a useful event
[ep,errb,ci] = mep(fp,tevents_ticks,200,200);

% plot the results
figure,hold on
plot(-pre_ticks:post_ticks,ep,'k')
plot(-pre_ticks:post_ticks,ep-errb,'k:')      % EP+SEM
plot(-pre_ticks:post_ticks,ep+errb,'k:')      % EP-SEM
xlabel('Ticks re: event'),ylabel('ADC units')

% add confidence intervals to plot
plot([-pre_ticks,post_ticks],[ci(1),ci(1)],'k:')    % lower bound
plot([-pre_ticks,post_ticks],[ci(2),ci(2)],'k:')    % upper bound

return

%% mlat.m
% [mean1stspike,median1stspike,peaklat,meantspike,mediantspike] = mlat(tspx,tstims,texclude)
% provides various measures of response latency:
% - mean first-spike latency
% - median first-spike latency
% - time of peak firing rate
% - mean spike time in time window
% - median spike time in time window
% This function requires the Statistics and Machine Learning Toolbox.
% 
% Recommended reading:
% Sakata, S., & Harris, K. D. (2009). Laminar structure of spontaneous and sensory-evoked population activity in auditory cortex. Neuron, 64(3), 404–418.

clear all,close all,clc

load('unitForMLIBTesting_4_rat.mat')

tstims = tevents(cevents==1);     % get stimulus times (80-dB noise burst played via a loudspeaker)
texclude = 0.05;                          % take only spikes from the first 0.05 seconds

% get measures
[mean1stspike,median1stspike,peaklat,meantspike,mediantspike] = mlat(tspx,tstims,texclude);

% now let's compare these measures to the visual impression from a PSTH and SDF
psth = mpsth(tspx,tstims,'pre',100,'post',200,'binsz',1,'fr',1,'chart',2);
sdf  = msdf(psth,'exp',2);
subplot(212),hold on
plot(sdf(:,1),sdf(:,2),'LineWidth',2)
ylim([0,140])

plot(mean1stspike,135,'.','MarkerSize',10)
plot(median1stspike,130,'.','MarkerSize',10)
plot(peaklat,125,'.','MarkerSize',10)
plot(meantspike,120,'.','MarkerSize',10)
plot(mediantspike,115,'.','MarkerSize',10)
legend('','','Mean 1st','Median 1st','Peak lat','Mean t','Median t','Location','North')

return

%% mmi.m
% [MI,cMI1995,cMI1996] = mmi(s1,s2)
% computes mutual information (MI) between two spike count distributions
% cMI is corrected MI, uses correction formulas from two different Panzeri papers
% This function requires the Statistics and Machine Learning Toolbox.
% 
% Recommended reading:
% 1) Treves, A., & Panzeri, S. (1995). The Upward Bias in Measures of Information Derived from Limited Data Samples. Neural Computation, 7(2), 399–407.
% 2) Panzeri, S., & Treves, A. (1996). Analytical estimates of limited sampling biases in different information measures. 
%    Network: Computation in Neural Systems, 7(1), 87–107.

clear all,close all,clc

load('unitForMLIBTesting_4_rat.mat')

s1 = mnspx(tspx,tevents(cevents==201),300,0,1); % get vector of spike counts for the last 300 ms before event 201
s2 = mnspx(tspx,tevents(cevents==202),300,0,1); % get vector of spike counts for the last 300 ms before event 202

[MI,cMI1995,cMI1996] = mmi(s1,s2)   % returns MI uncorrected and corrected

% let's compare MI with the area under the ROC curve, another measure to index discriminability of two distributions (chance = 0.5)
auroc = mroc(s1,s2);

return

%% mmi_psych.m
% function [MI,d] = mmi_psych(S1R1,S1R2,S2R1,S2R2)
% computes mutual information (MI) and sensitivity (d') for behavioral data in a two-stimulus (or two-category), two-choice task
% MI allows to bring psychometric and neurometric discrimination performance to the same scale
% This function requires the Statistics and Machine Learning Toolbox.
% 
% Recommended reading:
% Vandevelde, J. R., Yang, J.-W., Albrecht, S., Lam, H., Kaufmann, P., Luhmann, H. J., & Stüttgen, M. C. (2023). Layer- and cell-type-specific 
% differences in neural activity in mouse barrel cortex during a whisker detection task. Cerebral Cortex, 33, 1361–1382.

clear all,close all,clc

% example: 100 correct and 20 incorrect S1 trials, 60 correct and 40 incorrect S2 trials
[MI,d] = mmi_psych(100,20,40,60)

% MI and d' are monotonically related:
p_hit = 0.6;p_fa = 0.3;   % true hit and false alarm rates
for i = 1:1000  % run 1,000 experiments
  ntrials1 = round(normrnd(100,40));    % number of S1 trials
  ntrials2 = round(normrnd(100,40));    % number of S2 trials
  nhits = binornd(ntrials1,p_hit);      % number of hits in this experiment
  nfa   = binornd(ntrials2,p_fa);       % number of false alarms in this experiment
  [MI(i,1),d(i,1)] = mmi_psych(nhits,ntrials1-nhits,nfa,ntrials2-nfa);  % MI and d'
end

r = corr(d, MI,'rows','complete');
figure('name','Sensitivity (d'') vs. mutual information'),hold on
scatter(d,MI,'.')
xlabel('d prime')
ylabel('MI')
text(0.2,0.5,['r=',num2str(r,'%2.2f')])

return

%% mmultregress.m
% [b,p,sr,psr] = mmultregress(tspx,tevents,predictors,winpos,winsize,varargin)
% Performs multiple regression analysis on spike counts with moving windows of specified durations and at multiple positions relative to an event.
% If interactions between predictors are to be included, they should be specified within the predictors input matrix.
% The function additionally computes semipartial correlation coefficients to indicate the unique contribution of each predictor.
% This function requires the Statistics and Machine Learning Toolbox.
% 
% Recommended reading: 
% Kim, H., Sul, J. H., Huh, N., Lee, D., & Jung, M. W. (2009). Role of striatum in updating values of chosen actions. Journal of Neuroscience, 29(47), 14701–14712.

clear all,close all,clc

load unitForMLIBTesting_6_rat.mat

% perform multiple regression analysis of spike counts on three predictors: stimulus (1 or 2), response (1 or 2), outcome (reward or not: 0 or 1)
% we will use effects coding for predictors; recode 1,2 -> -1,+1 (stims, responses) and 0,1 -> -1,+1 (outcomes)
pred_stims = stims-1;pred_stims(pred_stims==0) = -1;
pred_resps = responses-1;pred_resps(pred_resps==0) = -1;
pred_outs  = outcomes;pred_outs(pred_outs==0) = -1;

predictors = [pred_stims,pred_resps,pred_outs];

winpos  = -0.4:0.025:0.8;
winsize = [0.025,0.1,0.25];

[b,p,sr,psr] = mmultregress(tspx,tstims,predictors,winpos,winsize,'plotit','prednames',{'stimulus','response','reward'},'normalize');

% The two panels on the left show results for a winsize of 0.1 s, the two plots on the right for a winsize of 0.2 s.
% The top panels show regression (b) coefficients as functions of window position.
% The bottom panels show the fraction of regression coefficients turning out significant (p<.05).
% This particular unit from rat medial frontal cortex encodes the animal's response, starting around 0.35 s / 0.25 s after stimulus onset
% (estimated from window sizes of 0.1 and 0.2 s, respectively).
% For higher temporal precision, decrease the interval between window positions and the size of the analysis window.
% Note that results from multiple regression analysis can be tricky to interpret, see Cohen et al. (2003),"Applied multiple regression/correlation 
% analysis for the behavioral sciences (3rd ed.). Taylor & Francis Group.

return

%% mnspx.m
% nspx = mnspx(spxtimes,trigtimes,pre,post,removenans)
% counts the number of spikes in a specified time window relative to a specific event
% the output argument (nspx) can e.g. be used as one of two inputs to mroc.m or mmi.m or ttest2
% 
% Recommended reading:
% Counting spikes is done in nearly every paper which reports spike recordings.

clear all,close all,clc

load unitForMLIBTesting_1_rat.mat

% get the distribution of spike counts relative to the start of event 1 (stimulus 1) and the following 100 ms
nspx = mnspx(tspx,tevents(cevents==1),0,100);

% plot the distribution
histogram(nspx,-0.5:1:max(nspx)+1)
xlabel('N Spikes'),ylabel('Frequency')

return

%% mpsth.m
% [psth,trialspx] = mpsth(spxtimes,trigtimes,varargin)
% constructs a per-stimulus time histogram
% optional plotting of PSTH and raster display
% 
% Recommended reading:
% 1) Abeles, M. (1982). Quantification, smoothing, and confidence limits for single-units’ histograms. Journal of Neuroscience Methods, 5(4), 317–325.
% 2) Stoilova, V. v, Knauer, B., Berg, S., Rieber, E., Jäkel, F., & Stüttgen, M. C. (2020). Auditory cortex reflects goal-directed movement but is not 
%    necessary for behavioral adaptation in sound-cued reward tracking. Journal of Neurophysiology, 124(4), 1056–1071.

clear all,close all,clc

load unitForMLIBTesting_2_pigeon

% return psth and trialspx
[psth,trialspx] = mpsth(tspx,tevents(cevents==1));
plot(psth(:,1),psth(:,2))

% again, but directly plot PSTH and raster
mpsth(tspx,tevents(cevents==1),'chart',2);

% same, but scale in units of firing rate
mpsth(tspx,tevents(cevents==1),'chart',2,'fr',1);

% same, but with another bins 10 ms wide with a different time range
mpsth(tspx,tevents(cevents==1),'chart',2,'fr',1,'binsz',10,'pre',500,'post',500);

return

%% mroc.m
% [auroc,ciLoUp] = mroc(x,y,varargin)
% computes the area under the ROC curve (non-parametrically) for two samples, e.g. spike count distributions
% note that mroc.m works faithfully only with vectors containing exclusively integers
% This function requires the Statistics and Machine Learning Toolbox.
% 
% Recommended reading:
% 1) Barlow, H., Levick, W., & Yoon, M. (1971). Responses to single quanta of light in retinal ganglion cells of the cat. Vision Research, 3, 87–101.
% 2) Britten, K. H., Shadlen, M. N., Newsome, W. T., & Movshon, J. A. (1992). The analysis of visual motion: a comparison of neuronal and psychophysical 
%    performance. Journal of Neuroscience, 12(12), 4745–4765.
% 3) Parker, A. J., & Newsome, W. T. (1998). Sense and the single neuron: probing the physiology of perception. Annual Review of Neuroscience, 21, 227–277.

clear all,close all,clc

% first generate two vectors with each element being a spike count
spx1 = poissrnd(3,25,1);
spx2 = poissrnd(7,25,1);

% compute area under the ROC curve
auroc = mroc(spx2,spx1);

% same with plot and bootstrapped confidence intervals
confidence_limits = 0.95;
nboot = 1000;

[auroc,bp,ciLoUp] = mroc(spx2,spx1,'plotit',1,'ci',[confidence_limits,nboot])

return

%% msdf.m
% [sdf,kernel] = msdf(psth,ftype,w,varargin)
% convert PSTH to spike-density function; PSTH should be in units of spikes/sec
% see script msdf_test.m for more examples
% This function requires the Statistics and Machine Learning Toolbox.
% 
% Recommended reading:
% 1) Abeles, M. (1982). Quantification, smoothing, and confidence limits for single-units’ histograms. Journal of Neuroscience Methods, 5(4), 317–325.
% 2) Shinomoto, S. (2010). Estimating the Firing Rate. In S. Grün & S. Rotter (Eds.), Analysis of Parallel Spike Trains (1st ed., pp. 21–35). Springer US.

clear all,close all,clc

load unitForMLIBTesting_2_pigeon

% first generate a PSTH with mpsth; scale in units of firing rate
psth = mpsth(tspx,tevents(cevents==1),'fr',1);

% run msdf, convolve PSTH with a Gaussian kernel with a standard deviation of 50 ms
ftype = 'Gauss';
w     = 50;
sdf   = msdf(psth,ftype,w);

% run msdf, convolve PSTH with a boxcar kernel with a width of 20 ms
ftype = 'boxcar';
w     = 20;
sdf2  = msdf(psth,ftype,w);

% plot the results
figure
subplot(211)
bar(psth(:,1),psth(:,2),'k')
subplot(212),hold on
plot(sdf(:,1),sdf(:,2))
plot(sdf2(:,1),sdf2(:,2))
xlabel('Time (ms)'),ylabel('Spikes per s')

return

%% mspectan.m
% [P1,f] = mspectan(signal,Fs,n,doplot)
% computes and plots the fast fourier transform of a signal with sampling rate Fs
% This function requires the Curve Fitting Toolbox.
% 
% Recommended reading:
% 1) Destexhe, A., Bédard, C. (2015). Local Field Potentials (LFP). In: Jaeger, D., Jung, R. (eds) Encyclopedia of Computational Neuroscience. Springer, New York, NY.
% 2) Buzsáki, G., Anastassiou, C. a, & Koch, C. (2012). The origin of extracellular fields and currents--EEG, ECoG, LFP and spikes. Nature 
%    Reviews Neuroscience, 13(6), 407–420.
% 3) Einevoll, G. T., Kayser, C., Logothetis, N. K., & Panzeri, S. (2013). Modelling and analysis of local field potentials for studying the function of 
%    cortical circuits. Nature Reviews Neuroscience, 14(11), 770–785.
% 
% Many more functions for EEG analysis can be found at https://www.fieldtriptoolbox.org/.

clear all,close all,clc

load data4mspectan

% first, plot the raw data and the spectrogram
% increases in signal amplitude correspond to alpha wave when the subject closed his eyes
figure('name','Raw signal','units','normalized','position',[.2,.2,.4,.4])

subplot(211)
t = 1/fs:1/fs:(length(eeg)/fs);   % time vector (seconds)
plot(t/60,eeg)
ylabel('Voltage (A.U.)')

subplot(212)
spectrogram(eeg,500,450,2^12,fs,'yaxis'); % [s,f,t] = spectrogram(x,window,noverlap,nfft,fs) 
axis tight
ylim([0,100])
set(gca,'CLim',[-80,-40])

% second, compute and plot the amplitude spectrum
n = 2^12;      % for this signal, 2^12 is nice, but the maximum is 0.5*length(y)
[P1,f] = mspectan(eeg,fs,n,1);

return

%% mtune.m
% [avFR,semFR,fitobject] = mtune(tspx,tevents,cevents,tpost)
% computes a tuning curve of a spike response to several stimuli and fits it with a normal curve
% 
% Recommended reading:
% Henry, G. H., Dreher, B., Bishop, P. O. (1974). Orientation specificity of cells in cat striate cortex. Journal of Neurophysiology, 37(6), 1394–1409.

% First, we simulate a series of spikes with a fixed (baseline) firing rate of 5 Hertz.
% We have 7 different stimuli which increase the firing rate of the neuron to different degrees.
% Each stimulus lasts for 1 s and is presented 25 times. The inter-stimulus interval is 4 s.

clear all,close all,clc

blfr     = 5;                 % baseline firing rate (in Hertz)
nstims   = 7;                 % number of stimuli
frinc    = [4,8,10,8,4,2,1];  % firing rate changes during stimulation, for each stimulus
stimdur  = 1;                 % stimulus duration (in seconds)
isi      = 4;                 % inter-stimulus interval
trialsPerStim = 25;           % number of trials per stimulus

% generate random stimulus sequence
ntrials  = nstims*trialsPerStim;
stimseq = repmat((1:nstims)',trialsPerStim,1);
r = randperm(numel(stimseq));
cevents = stimseq(r);clear r stimseq
tevents = isi+stimdur:isi+stimdur:ntrials*(isi+stimdur);   % stimulus time stamps

% generate spike times in inter-trial intervals
tspx = [];
for i = 1:ntrials+1
  tspx_dummy = cumsum(exprnd(1/blfr,blfr*10,1)+0.003);  % spike train with 3-ms refractory period
  tspx = [tspx;tspx_dummy(tspx_dummy<isi+stimdur)+i*(isi+stimdur)];
end

% increase firing rates during stimulation
for i = 1:ntrials
  tspx_dummy = cumsum(exprnd(1/(frinc(cevents(i))),stimdur*blfr*10,1)+0.003);
  tspx = [tspx;tspx_dummy(tspx_dummy<stimdur)+i*(isi+stimdur)];
end
tspx = sort(tspx);

% take a look at the simulated neurons' peri-stimulus time histogram using the MLIB function mpsth
figure('units','normalized','position',[.2,.1,.2,.8])
for i = 1:nstims
  psth = mpsth(tspx,tevents(cevents==i),'pre',500,'post',1500,'binsz',100);
  subplot(nstims,1,i),title(['Stimulus no. ',num2str(i)]),hold on
  bar(psth(:,1)+50,psth(:,2),'k')
  plot([0,0],[0,1000],'k:')
  plot([1000*stimdur,1000*stimdur],[0,1000],'k:')
  ylim([0,blfr*max(frinc)])
end
xlabel('Time re: stimulus onset (ms)')
ylabel('Spike count')

% now construct and plot tuning curve
% note that a Gaussian fit to the tuning curve is produced and plotted only if the Curve Fitting Toolbox is available
mtune(tspx,tevents,cevents,stimdur*1000);

return

%% mvecstrength.m
% [vs,pvs] = mvecstrength(cycle_tspx,T) 
% computes the vector strength of a spike response to a periodic signal basis for computation is the cycle time histogram, i.e., 
% a PSTH relative to the onset of a sine-wave cycle which lasts until its end
% 
% Recommended reading:
% Goldberg, J. M., & Brown, P. B. (1969). Response of Binaural Neurons of Dog Superior Olivary Complex to Dichotic Tonal Stimuli: Some Physiological 
% Mechanisms of Sound Localization. Journal of Neurophysiology, 32(4), 613–636.

clear all,close all,clc

% first, we simulate a series of spikes relative to the onset of a sine wave stimulus
nspikes   = 100;                % number of spikes fired by the unit
frequency = 10;                 % frequency of sine wave in Hertz
T         = 1000/frequency;     % period length of sine wave in ms
SD        = 10;                  % standard deviation of the Gaussian kernel; the larger SD, the higher the vector strength

cycle_tspx = normrnd(T/2,SD,nspikes,1);          % generate spike times from Gaussian distribution
                                                 % this synthetic unit will fire around the peak of the sine wave
cycle_tspx(cycle_tspx<0 | cycle_tspx>T) = [];    % keep only spikes within a single cycle

cth = histcounts(cycle_tspx,0:T);                % construct cycle time histogram for plotting

[vs,pvs] = mvecstrength(cycle_tspx,T)            % calculate vector strength and its p-value

figure('units','normalized','position',[.3,.1,.4,.8])
% first plot the simulated spike train relative to the sine wave stimulus
subplot(411),title(['Simulated spike train for continuous stimulation, SD=',num2str(SD)]),hold on
ncycles = ceil(2*nspikes/frequency);      % just for illustration purposes
c = 1;
for i = 1:ncycles
  nspxthiscycle = randi(3);
  plot((i-1)*T+cycle_tspx(c:c+nspxthiscycle-1),0.1,'k.')
  c = c+nspxthiscycle;
end
plot(linspace(0,T*ncycles,500),0.2*(sin(linspace(-0.5*pi,(ncycles-1)*2*pi+1.5*pi,500))+1),'Color',ones(1,3)*0.7)
ylim([-0.1,0.6])
axis off

% plot the cycle time histogram
% add the sine wave cycle for illustration purposes
subplot(412),title(['VS=',num2str(vs,'%2.2f'),', p=',num2str(pvs,'%2.3f')]),hold on
bar(0.5:T-0.5,cth)
plot(linspace(0,T,100),0.5*max(cth).*(sin(linspace(-0.5*pi,1.5*pi,100))+1),'Color',ones(1,3)*0.7)
xlabel('Time (ms) within cycle'),ylabel('Spike count')

% now do the same with randomly distributed spike times
cycle_tspx = rand(nspikes,1)*T;
cth = histcounts(cycle_tspx,0:T);    % cycle time histogram
[vs,pvs] = mvecstrength(cycle_tspx,T)

subplot(413),title(['Simulated spike train for continuous stimulation, random spike times']),hold on
c = 1;
for i = 1:ncycles
  nspxthiscycle = randi(3);
  plot((i-1)*T+cycle_tspx(c:c+nspxthiscycle-1),0.1,'k.')
  c = c+nspxthiscycle;
end
plot(linspace(0,T*ncycles,500),0.2*(sin(linspace(-0.5*pi,(ncycles-1)*2*pi+1.5*pi,500))+1),'Color',ones(1,3)*0.7)
ylim([-0.1,0.6])
axis off

subplot(414),title(['VS=',num2str(vs,'%2.2f'),', p=',num2str(pvs,'%2.3f')]),hold on
bar(0.5:T-0.5,cth)
plot(linspace(0,T,100),0.5*max(cth).*(sin(linspace(-0.5*pi,1.5*pi,100))+1),'Color',ones(1,3)*0.7)
xlabel('Time (ms) within cycle'),ylabel('Spike count')

return

%% mwa.m
% [mat,optmax,optmin,aurocs,cis] = mwa(tspx,tevents,trefevents,winpos,winsize,varargin)
% computes the area under the ROC curve for spike count distributions for a range of analysis window durations and window positions
% 
% For an example, see Figure 3 and associated manuscript text in:
% Stüttgen MC, Schwarz C (2008) Psychophysical and neurometric detection performance under stimulus uncertainty. Nature Neuroscience 11: 1091.
% 
% This function requires the Statistics and Machine Learning Toolbox.

clear all,close all,clc

load unitForMLIBTesting_4_rat.mat

window_positions = -0.1:0.01:0.2;     % analysis window positions in seconds relative to the event - we shift analysis windows by 10 ms
window_sizes     = [0.005,0.01,0.02,0.03,0.05:0.025:0.1];     % analysis window durations in seconds, ranging from 5 to 100 ms

tstim1     = tevents(cevents==1);       % get timestamps of stimulus 1
trefevents = tstim1-0.5;                % our reference events are simply half a second before each stimulus 1 onset

[mat,optmax,optmin,aurocs,cis] = mwa(tspx,tstim1,trefevents,...
                             window_positions,window_sizes,'plotit',1,'ci',[95,500],'smoothit',1);    % smoothing affects optmax and optmin!
return

%% mwa2.m
% [mat,optmax,optmin,aurocs,cis] = mwa2(tspx,tevents,tevents2,winpos,winsize,varargin)
% similar to mwa.m, but computes the area under the ROC curve for two different conditions (e.g., stimuli) while in
% mwa.m there is only one reference distribution and the signal distribution shifts, in mwa2.m both distributions shift
% for example, with mwa2.m you can compare spike count distributions from 0 to 20 ms relative to stimulus 1 and stimulus 2, 
% then # 20 to 40 ms for stimulus 1 and stimulus 2, and so on
% 
% For an example, see Figure 3 and associated manuscript text in:
% Stüttgen MC, Schwarz C (2008) Psychophysical and neurometric detection performance under stimulus uncertainty. Nature Neuroscience 11: 1091.
% 
% This function requires the Statistics and Machine Learning Toolbox.

clear all,close all,clc

load unitForMLIBTesting_5_rat.mat

tS1 = tevents(cevents==1);
tS2 = tevents(cevents==2);

% first visualize the responses of this unit to the two different stimuli
psth1 = mpsth(tspx,tS1,'fr',1,'binsz',10,'pre',500,'post',500);
psth2 = mpsth(tspx,tS2,'fr',1,'binsz',10,'pre',500,'post',500);
figure
bar(psth1(:,1),psth1(:,2),'FaceAlpha',0.5,'EdgeAlpha',0,'BarWidth',1),hold on
bar(psth2(:,1),psth2(:,2),'FaceAlpha',0.5,'EdgeAlpha',0,'BarWidth',1)
xlabel('Time re: stimulus onset (ms)'),ylabel('Firing rate')
legend('S1','S2')

% now, use mwa2.m to compare the responses window by window
winpos  = -0.2:0.01:0.4;
winsize = [0.005,0.01,0.02,0.05,0.075,0.1,0.15,0.2];
[mat,optmax,optmin,aurocs,cis] = mwa2(tspx,tS2,tS1,winpos,winsize,'plotit',2);

% mwa2 will return an 8-by-61 matrix with 488 AUROC values, computed from -0.2 to +0.4 seconds relative to the specified events, 
% with window sizes ranging from 5 ms to 200 ms, and generates a fancy plot on top

return

%% mwave.m
% waveparms = mwave(meanwave,si,varargin)
% function computes the a range of waveform indices, such as the widths of the first and second response peaks of 
% a spike waveform (full width at half maximum), the trough-to-peak ratio and its duration etc.
% 
% Recommended reading:
% 1) McCormick, D. A., Connors, B. W., Lighthall, J. W., & Prince, D. A. (1985). Comparative Electrophysiology of Pyramidal and Sparsely Spiny Stellate 
%    Neurons of the Neocortex. Journal of Neurophysiology, 54(4), 782–806.
% 2) Trainito, C., von Nicolai, C., Miller, E. K., & Siegel, M. (2019). Extracellular Spike Waveform Dissociates Four Functionally Distinct Cell Classes in 
%    Primate Cortex. Current Biology, 29(18), 2973-2982.e5.
% 3) Quiroga, R. Q. (2009). What is the real shape of extracellular spikes? Journal of Neuroscience Methods, 177(1), 194–198.

clear all,close all,clc

load('unitForMLIBTesting_2_pigeon.mat')

fs   = 2*10^4;      % sampling frequency in Hertz, not provided in the file

waveparms = mwave(mean(adc),fs,'plot');

% waveparms are (in this order):
% - full width at half maximum of the positive (1) and the negative (2) waveform peaks, both in microseconds
% - peak-to-trough amplitude (3)
% - peak-to-rough duration (4) in microseconds
% - peak-to-trough ratio (positive peak divided over negative peak, absolute magnitudes) (5)

return

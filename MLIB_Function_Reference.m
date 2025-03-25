%% FUNCTION REFERENCE MLIB
% 
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
% Vandevelde JR, Yang JW, Albrecht S, Lam H, Kaufmann P, Luhmann HJ, Stüttgen MC (2023)
% Layer- and cell-type-specific differences in neural activity in mouse barrel cortex during a whisker detection task.
% Cerebral Cortex 33: 1361.
% 
% Yeganeh F, Knauer B, Guimaraes Backhaus R, Yang JW, Stroh A, Luhmann HJ, Stüttgen MC (2022)
% Effects of optogenetic inhibition of a small fraction of parvalbumin-positive interneurons on the representation of sensory stimuli in mouse barrel cortex.
% Scientific Reports 12: 19419.
% 
% There are several example files to illustrate the use of the MLIB functions for extracellular spike recordings.
 
% File 1: unitForMLIBTesting_1_rat.mat
% This file contains a single structure 'spx' holding waveforms of a nice single unit as well as event markers (spx.cevents) and their respective 
% time stamps (spx.tevents). Sampling rate (spx.fs) was 20 kHz. Recordings were made in the frontal cortex of a Long Evans rat performing a two-stimulus, 
% two-response conditional discrimination task.  The task is described in Stoilova et al., 2020. During the first 120 trials, all correct responses were 
% rewarded, and all incorrect responses were inconsequential. Time stamps are in seconds.
 
% File 2: unitForMLIBTesting_2_pigeon.mat
% This file contains a single structure 'spx' holding 853 waveforms of a nice single unit (spx.adc, 32 ticks) and their time stamps (spx.tspx) as well as 
% event markers (spx.cevents) and their time stamps (spx.tevents). The data is from a six-stimulus single-interval forced choice task as described in 
% Lengersdorf et al., 2014. The sampling rate is 16.67 kHz. Time stamps are in seconds.
% Codes in spx.cevents represent: 1: trial initiation, 2-7: six visual stimuli, 8: choice keys on, 9: food reward, 10: food omission, 11: punishment, 
% 12-14: key pecks to the left, center, and right pecking keys, respectively.

% File 3: unitForMLIBTesting_3_pigeon.mat
% Similar to file 2 above, but containing two spike channels (spxCh9 and spxCh10) recorded in the same session. Both spike channels hold a vector of ALL 
% detected waveforms on a given channel, i.e., spikes, noise, and artifacts.
% spxCh9 marker 0 contains trash, marker 1 is a good single unit and marker 2 is a medium-quality single unit; 
% spxCh10 marker 0 contains trash, markers 1 through 3 contain good SUs
% spx.tevents and spx.cevents contain event time stamps and their code, respectively.
% Codes represent: 97: initiation, 98-103: six visual stimuli, 104: choice keys on, 105: food, 106: food omission, 107: punishment, 
% 108-110: key pecks to the left, center, and right pecking keys, respectively.

% File 4: unitForMLIBTesting_4_rat.mat
% Similar to file 1, but the unit is from rat auditory cortex.

% File 5: unitData4mcluster.mat
% This file contains data from Kasties et al. (2016). Neurons were recorded while animals viewed one out of five visual stimuli (each stimulus presentation 
% lasted 5 seconds). For each neuron and stimulus, a PSTH was constructed, and the five PSTHs were concatenated before clustering (variable app_psths for 
% 98 neurons). Each PSTH had 25 data points (firing rate in consecutive 0.2-second bins). The upshot of the exercise using this file is that units can be 
% clustered according to their preferred stimulus.
% The variable zpsths contains the (unit-wise) z-scored PSTHs.

% File 6: data4mcsd.mat
% This file serves to exemplify the use of mcsd.m for computing a current source density analysis. The only variable mean_FP contains LFP recordings 
% from a linear silicon probe with 20 electrode sites spaced 50 microns apart. Sampling rate was 1000 Hz, only the first 80 ms following stimulus onset 
% (deflection of the principal whisker) are included. Recordings were performed in mouse barrel cortex. This data was kindly provided by Dr. Jenq-Wei Yang
% and Prof. Heiko Luhmann, University Medical Center Mainz, Germany.

% File 7: data4mdimreduct.mat
% This file serves to illustrate the use of mdimreduct.m to plot neural population trajectories (e.g., Saxena & Cunningham, Curr Opin Neurobiol 2019).
% The matrices nspx_s1 and nspx_s2 each contain spike counts from the same 300 units, recorded in rat dorsomedial frontal cortex during a two-stimulus, 
% two-choice task. Rows are units, columns are 17 100-ms time bins, shifted in steps of 50 ms, representing window positions from -400 ms to +400 ms 
% relative to the onset of stimulus 1 (nspx_s1) or stimulus 2 (nspx_s2). Only correct responses were included in averages.

% File 8: eeg_alpha.mat
% This file is used to illustrate the use of mspectan.m to calculate and plot the fast Fourier transform. An EEG recording was performed over occipital 
% cortex, the subject opened and closed his eyes every few seconds.
% The vector y contains the EEG signal which was sampled at a frequency of 1000 Hertz (Fs).
 
%% Overview over MLIB functions in alphabetical order

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
    
%% Function mcc.m
% [L_ratio,ID,A,J3,PsF,DBVI] = mcc(X,cid,whichclusters,ts,whichPCs,nocorrs)
% function to assess the quality of spike sorting (post hoc)
% computes standard quality measures (L_ratio, J3,...); see doc mcc for more information
% This function requires the Statistics and Machine Learning Toolbox.
% 
% Recommended reading:
% 1) Schmitzer-Torbert, N., Jackson, J., Henze, D., Harris, K., & Redish, A. D. (2005). Quantitative measures of cluster quality for 
%    use in extracellular recordings. Neuroscience, 131(1), 1–11.
% 2) Nicolelis, M. A. L., Dimitrov, D., Carmena, J. M., Crist, R., Lehew, G., Kralik, J. D., & Wise, S. P. (2003). Chronic, multisite, 
%    multielectrode recordings in macaque monkeys. PNAS 100(19), 11041–11046.
 
clear all,close all,clc
 
load unitForMLIBTesting_3_pigeon
 
% an electrode with two clusters, take the first three principal components
% field 'adc' contains the waveforms of both units, field 'markers' contains their cluster ID,
% field 'tspx' contains spike time stamps
[L_ratio,ID,A,J3,PsF,DBVI] = mcc(spxCh9.adc,spxCh9.markers,[1,2],spxCh9.tspx,1:3);
 
% take another electrode from which three clusters were isolated
% use the first 8 PCs for computation of L_ratio etc., omit crosscorrelation figure
[L_ratio,ID,A,J3,PsF,DBVI] = mcc(spxCh10.adc,spxCh10.markers(:,1),[1,2,3],spxCh10.tspx,1:8,'nocorrs');
 
% Figure 1 visualizes cluster separability
% subplot(331) shows the different waveforms with clusters color-coded
% subplot(332) shows the percentage of variance explained by principal components (after performing PCA on all waveforms)
% subplot(333) shows PC1 for all spikes (clusters color-coded) as a function of recording time
% subplots 4, 5, 6 show spikes in principal component space; bold lines show each cluster's mean and 2 standard deviations in 2D PC space
% subplot(337) shows the differences of all spikes' Mahalanobis distances from the two cluster means; the less overlap, the better
% subplot(338) shows a receiver operating characteristic (ROC) curve. The measure A' gives the probability to correctly categorize two random samples taken from each distribution
% subplot(339) presents some commonly used cluster separation quality indices
% 
% Figure 2 shows autocorrelograms (AC) and crosscorrelograms (CC) of the different units ACs of clean single units have very 
% few or no spikes during the first few milliseconds after each spike (refractory period).
 
return
 
%% Function mcheck.m
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
c = mcheck(spx.adc);

% input waveforms and waveform timestamps
c = mcheck(spx.adc,'timevec',spx.tspx);

% all of the above, plus event timestamps
c = mcheck(spx.adc,'timevec',spx.tspx,'tevents',spx.tevents(spx.cevents==13));

% now try the other file
clear all
load unitForMLIBTesting_3_pigeon.mat

% input all waveforms from channel 9; the result looks messy
c = mcheck(spxCh9.adc);

% now input only waveforms which have been classfied into cluster 1 during spike sorting
c = mcheck(spxCh9.adc,'markvec',spxCh9.markers,'target',1);

% as above, but additionally provide spike timestamps
c = mcheck(spxCh9.adc,'markvec',spxCh9.markers,'target',1,'timevec',spxCh9.tspx);

% as above, but also provide timestamps for a PSTH (relative to pecks on the center key) and sampling rate
c = mcheck(spxCh9.adc,'markvec',spxCh9.markers,'target',1,'timevec',spxCh9.tspx,'tevents',tevents(cevents==109),'fs',10^6/spxCh9.sampleinterval);

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

%% Function mcluster.m
% [cid,avsil] = mcluster(data,z,nclus,whichPCs,cType,whichDist,metric,T,plotting)
% performs clustering on a data matrix after running PCA and returns cluster codes cid and average silhouette values
% input 'data' is a matrix where rows are units and columns are binned firing rates or eta2/omega2 values
% This function requires the Statistics and Machine Learning Toolbox.
% 
% Recommended reading:
% Kasties, N., Starosta, S., Güntürkün, O., & Stüttgen, M. C. (2016). Neurons in the pigeon caudolateral nidopallium differentiate 
% Pavlovian conditioned stimuli but not their associated reward value in a sign-tracking paradigm. Scientific Reports, 6(October), 35469.

clear all,close all,clc

load('unitData4mcluster.mat')

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
  subplot(1,nclus,i),hold on
  for j = 1:5
    plot(0.1:0.2:4.9,mean(zpsths(cid==i,(j-1)*25+1:j*25)))
  end
  axis([0,5,-1.5,3])
  if i==1,xlabel('Time (s)'),ylabel('z(FR)'),end
end
legend('Stim 1','Stim 2','Stim 3','Stim 4','Stim 5','units','normalized','Position',[0.45,0,1,1])

return

%% Function mcsd.m
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

[csd_r,csd_i,lfp_s] = mcsd(mean_FP,electrode_spacing,doplot,time_vec,theMinus);

% subplot(221) shows the raw LFP traces, color-coded for depth from black (top) to white (bottom)
% subplot(222) shows the LFP traces both as black lines and color-coded (after smoothing)
% subplot(223) shows the raw CSD
% subplot(224) shows the CSD after 10x interpolation and spatial smoothing

% Since the CSD is the spatial derivative of the LFP, there two rows less in csd_r than in mean_FP.
set(gcf,'units','normalized','Position',[.2,.2,.6,.6])

return

%% Function mcvar.m
% cv = mcvar(x,mode)
% computes the coefficient of variation for interspike intervals and the Fano Factor for a spike count distribution
% 
% Recommended reading:
% Maimon, G., & Assad, J. A. (2009). Beyond Poisson: increased spike-time regularity across primate parietal cortex. Neuron, 62(3), 426–440.

clear all,close all,clc

load('unitForMLIBTesting_1_rat.mat')

cv_isi = mcvar(spx.tspx,'diff')    % field tspx holds spike times in seconds relative to the beginning of the recording

spike_counts = histcounts(spx.tspx,0:999); % get spike counts for consecutive 1-s bins over the first 1000 seconds of recording
cv_ff = mcvar(spx.tspx,'abs')

return

%% Function mdimreduct.m
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

%% Function mep.m
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

%% Function mlat.m
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

tspx = spx.tspx;                          % get spike times
tstims = spx.tevents(spx.cevents==1);     % get stimulus times (80-dB noise burst played via a loudspeaker)
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

%% Function mmi.m
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

s1 = mnspx(spx.tspx,spx.tevents(spx.cevents==201),300,0,1); % get vector of spike counts for the last 300 ms before event 201
s2 = mnspx(spx.tspx,spx.tevents(spx.cevents==202),300,0,1); % get vector of spike counts for the last 300 ms before event 202

[MI,cMI1995,cMI1996] = mmi(s1,s2)   % returns MI uncorrected and corrected

% let's compare MI with the area under the ROC curve, another measure to index discriminability of two distributions (chance = 0.5)
auroc = mroc(s1,s2);

return

%% Function mmi_psych.m
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

%% Function mmultregress.m
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

[b,p,sr,psr] = mmultregress(spx.tspx,tstims,predictors,winpos,winsize,'plotit','prednames',{'stimulus','response','reward'},'normalize');

% The two panels on the left show results for a winsize of 0.1 s, the two plots on the right for a winsize of 0.2 s.
% The top panels show regression (b) coefficients as functions of window position.
% The bottom panels show the fraction of regression coefficients turning out significant (p<.05).
% This particular unit from rat medial frontal cortex encodes the animal's response, starting around 0.35 s / 0.25 s after stimulus onset
% (estimated from window sizes of 0.1 and 0.2 s, respectively).
% For higher temporal precision, decrease the interval between window positions and the size of the analysis window.
% Note that results from multiple regression analysis can be tricky to interpret, see Cohen et al. (2003),"Applied multiple regression/correlation 
% analysis for the behavioral sciences (3rd ed.). Taylor & Francis Group.

return

%% Function mnspx.m
% nspx = mnspx(spxtimes,trigtimes,pre,post,removenans)
% counts the number of spikes in a specified time window relative to a specific event
% the output argument (nspx) can e.g. be used as one of two inputs to mroc.m or mmi.m or ttest2
% 
% Recommended reading:
% Counting spikes is done in nearly every paper which reports spike recordings

clear all,close all,clc

load unitForMLIBTesting_1_rat.mat

% get the distribution of spike counts relative to the start of event 1 (stimulus 1) and the following 100 ms
nspx = mnspx(spx.tspx,spx.tevents(spx.cevents==1),0,100);

% plot the distribution
histogram(nspx,-0.5:1:max(nspx)+1)
xlabel('N Spikes'),ylabel('Frequency')

return

%% Function mpsth.m
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
[psth,trialspx] = mpsth(spx.tspx,spx.tevents(spx.cevents==1));
plot(psth(:,1),psth(:,2))

% again, but directly plot PSTH and raster
mpsth(spx.tspx,spx.tevents(spx.cevents==1),'chart',2);

% same, but scale in units of firing rate
mpsth(spx.tspx,spx.tevents(spx.cevents==1),'chart',2,'fr',1);

% same, but with another bins 10 ms wide with a different time range
mpsth(spx.tspx,spx.tevents(spx.cevents==1),'chart',2,'fr',1,'binsz',10,'pre',500,'post',500);

return

%% Function mroc.m
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
spx1 = poissrnd(4,20,1);
spx2 = poissrnd(6,20,1);

% compute area under the ROC curve
auroc = mroc(spx2,spx1);

% same with plot and bootstrapped confidence intervals
confidence_limits = 0.95;
nboot = 1000;

[auroc,bp,ciLoUp] = mroc(spx2,spx1,'plotit',1,'ci',[confidence_limits,nboot])

return

%% Function msdf.m
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
psth = mpsth(spx.tspx,spx.tevents(spx.cevents==1),'fr',1);

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

%% Function mspectan.m
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

load eeg_alpha

% first, plot the raw data and the spectrogram
% increases in signal amplitude correspond to alpha wave when the subject closed his eyes
figure('name','Raw signal','units','normalized','position',[.2,.2,.4,.4])

subplot(211)
t = 1/Fs:1/Fs:(length(y)/Fs);   % time vector (seconds)
plot(t/60,y)
ylabel('Voltage (A.U.)')

subplot(212)
spectrogram(y,500,450,2^12,Fs,'yaxis'); % [s,f,t] = spectrogram(x,window,noverlap,nfft,fs) 
axis tight
ylim([0,100])
set(gca,'CLim',[-80,-40])

% second, compute and plot the amplitude spectrum
n = 2^12;      % for this signal, 2^12 is nice, but the maximum is 0.5*length(y)
[P1,f] = mspectan(y,Fs,n,1);

return

%% Function mtune.m
% [avFR,semFR,fitobject] = mtune(tspx,tevents,cevents,tpost)
% computes a tuning curve of a spike response to several stimuli and fits it with a normal curve
% 
% Recommended reading:
% Henry, G. H., Dreher, B., Bishop, P. O. (1974). Orientation specificity of cells in cat striate cortex. Journal of Neurophysiology, 37(6), 1394–1409.

% First, we simulate a series of spikes with a fixed rate of 5 Hertz.
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

% construct and plot tuning curve
% note that a Gaussian fit to the tuning curve is produced and plotted only if the Curve Fitting Toolbox is available
mtune(tspx,tevents,cevents,stimdur*1000);

return

%% Function mvecstrength.m
% [vs,pvs] = mvecstrength(tspx,T) 
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
SD        = 20;                 % the smaller the SD of the Gaussian, the higher the vector strength

tspikes = normrnd(50,SD,nspikes,1);     % generate spike times from Gaussian distribution
tspikes(tspikes<0 | tspikes>T) = [];    % keep only spikes within a single cycle

cth = histcounts(tspikes,0:100);        % construct cycle time histogram for plotting

[vs,pvs] = mvecstrength(tspikes,T)      % calculate vector strength and its p-value

figure
subplot(211),title(['VS=',num2str(vs,'%2.2f'),', p=',num2str(pvs,'%2.3f')]),hold on,bar(0.5:99.5,cth)

% now do the same with randomly distributed spike times
tspikes = rand(nspikes,1)*T;
cth = histcounts(tspikes,0:100);    % cycle time histogram
[vs,pvs] = mvecstrength(tspikes,T)

subplot(212),title(['VS=',num2str(vs,'%2.2f'),', p=',num2str(pvs,'%2.3f')]),hold on,bar(0.5:99.5,cth)
xlabel('Time (ms)'),ylabel('Spike count')

return

%% Function mwa.m
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

tevents    = spx.tevents(spx.cevents==1);
trefevents = spx.tevents(spx.cevents==1)-0.5;   % our reference events are simply half a second before each event

[mat,optmax,optmin,aurocs,cis] = mwa(spx.tspx,tevents,trefevents,...
                             window_positions,window_sizes,'plotit',1,'ci',[95,500],'smoothit',1);    % smoothing affects optmax and optmin!
return

%% Function mwa2.m
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

tS1 = spx.tevents(spx.cevents==1);
tS2 = spx.tevents(spx.cevents==2);

% first visualize the responses of this unit to the two different stimuli
psth1 = mpsth(spx.tspx,tS1,'fr',1,'binsz',10,'pre',500,'post',500);
psth2 = mpsth(spx.tspx,tS2,'fr',1,'binsz',10,'pre',500,'post',500);
figure
bar(psth1(:,1),psth1(:,2),'FaceAlpha',0.5,'EdgeAlpha',0,'BarWidth',1),hold on
bar(psth2(:,1),psth2(:,2),'FaceAlpha',0.5,'EdgeAlpha',0,'BarWidth',1)
xlabel('Time re: stimulus onset (ms)'),ylabel('Firing rate')
legend('S1','S2')

% now, use mwa2.m to compare the responses window by window
winpos  = -0.2:0.01:0.4;
winsize = [0.005,0.01,0.02,0.05,0.075,0.1,0.15,0.2];
[mat,optmax,optmin,aurocs,cis] = mwa2(spx.tspx,tS2,tS1,winpos,winsize,'plotit',2);

% mwa2 will return an 8-by-61 matrix with 488 AUROC values, computed from -0.2 to +0.4 seconds relative to the specified events, 
% with window sizes ranging from 5 ms to 200 ms, and generates a fancy plot on top

return

%% Function mwave.m
% waveParms = mwave(meanwave,si,varargin)
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

Fs   = 2*10^4;      % sampling frequency in Hertz
tick = 10^6/Fs;     % use Fs to compute the interval between consecutive samples (in microseconds)

waveParms = mwave(mean(spx.adc'),tick,'plot');

return

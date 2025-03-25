function [csd_r,csd_i,lfp_s] = mcsd(lfpmat,es,doplot,t,theMinus)
% function [csd_r,csd_i,lfp_s] = mcsd(lfpmat,es,doplot,t,theMinus)
%
% Performs a current-source density analysis from an input matrix lfpmat (containing LFP voltages) and generates an optional plot.
% The procedure is described in Reyes-Puerta et al. (2015, Cerebral Cortex), as well as in Sakata and Harris (2009, Neuron), going back to
% Freeman and Nicholson (1975, J Neurophysiol). The procedure described in Yang et al. (2017, Cerebral Cortex) is somewhat different.
% 
% INPUT
% lfpmat      2D matrix where rows are electrodes, columns are ticks, and elements are voltages (mV, usually averaged over multiple trials)
%             electrodes must be regularly spaced
% es          electrode spacing in mm (yes, mm)
%             Reyes-Puerta et al. used neighboring electrodes spaced 75 or 100 µm apart.
%             Sakata and Harris used neighboring electrodes with 50 µm inter-electrode distance.
%             Freeman and Nicholson used a stepsize of h=2; I omitted this from my code.
% doplot      to plot (1) or not to plot (0)
% t           time vector in ms (e.g. -100,-99,... for 1000 Hz)
% the Minus   multiply the 2nd spatial derivative with -1 or not (1/0, resp.)
% 
% OUTPUT
% csd_r     current source density profile of the recorded data in mV / mm^2
% csd_i     interpolated version of the CSD
% lfp_s     LFP after spatial and temporal downsampling and spatial smoothing
%
% EXAMPLE
% 
% see test_mcsd.m for an example
% 
% Big question: why are sinks sometimes positive, sometimes negative? Yang et al. include -1 as a factor, VRP possibly also.
% Differences to Yang et al.'s CSD figures can also be explained in that they calculated the CSD on the expanded LFP matrix 
% (duplicated top and bottom rows) WITHOUT spatial smoothing, but then linearly interpolated the CSD (10-fold)!
% This is correctly explained in his 2017 Cerebral Cortex paper (which also includes the -).
% Vicente clearly performed spatial smoothing on the LFP, not the CSD. Same for Sakata and Harris.
% Not having a -1 makes the code consistent with Sakata and Harris - positive values are sinks (depicted red).
% 
% HISTORY
% April 2021      adjusted subplot(224) - csd_in now ranges from electrode 2 to electrode n-1 (squares plotted on centers)
% September 2020  added argument 'theMinus' to allow for having the -1 as a factor or not
% March 2020      removed the expansion of the matrix which was previously included twice, once for spatial smoothing, once for the stepsize thing
% 
% by Maik C. Stüttgen, Sep 2017, UMC Mainz
%% duplicate uppermost and lowermost channels of lfpmat to allow spatial smoothing and CSD calculation
lfpmatexp = nan(size(lfpmat,1)+2,size(lfpmat,2));
lfpmatexp(2:end-1,:) = lfpmat;
lfpmatexp(1,:)   = lfpmatexp(2,:);
lfpmatexp(end,:) = lfpmatexp(end-1,:);
%% smooth LFPs across spatially adjacent channels to reduce high spatial-frequency noise components
% The formula is taken from Reyes-Puerta et al., who took it from Sakata and Harris. The smoothed LFP is given by:
% P(r) = 0.25 * (p(r+h) + 2*p(r) + p(r-h)),
% where p(r) is the LFP at depth r and h is the sampling interval (in mm).
% Freeman and Nicholson also give this simple spatial smoothing formula and call it a symmetrical, weighted average of the potential about
% a given point (p. 371).
lfp_s = nan(size(lfpmat));
for i = 2:size(lfpmatexp,1)-1
  lfp_s(i-1,:) = 0.25*(lfpmatexp(i+1,:) + 2*lfpmatexp(i,:) + lfpmatexp(i-1)); % simply a weighted mean
end
%% calculate the second spatial derivative
% The formula is taken from Reyes-Puerta et al., who took it from Sakata and Harris. The CSD is given by:
% D = (1/h^2) * (P(r+h) - 2*P(r) + P(r-h)).
% Freeman and Nicholson give this equation on p. 371 but say that it often gives noisy results.
% Fukui et al. (Sci Rep 2020) give the same formula but have the minus in front! They refer to Einevoll et al. (Nat Neurosci 2013) who also have it.
csd_r = nan(size(lfp_s,1)-2,size(lfp_s,2));
for i = 2:size(lfp_s,1)-1
  if theMinus==1
    csd_r(i-1,:) = -(1/(es^2)) * (lfp_s(i+1,:) - 2*lfp_s(i,:) + lfp_s(i-1,:)); % Yang et al. (2017) have the -, so do Zempeltzi and Happel (Comm Biol 2020)
  elseif theMinus==0
    csd_r(i-1,:) = (1/(es^2)) * (lfp_s(i+1,:) - 2*lfp_s(i,:) + lfp_s(i-1,:));
  end
end
% Zempeltzi and Happel write:
% CSD distributions reflect the local spatiotemporal current flow of positive ions from extracellular to intracellular space evoked by synaptic 
% populations in laminar neuronal structures. CSD activity thereby reveals the spatiotemporal sequence of neural activation across cortical layers 
% as ensembles of synaptic population activity21,53. One advantage of the CSD transformation is that it is reference-free and hence less affected 
% by far-field potentials and referencing artifacts. It allows to observe the local synaptic current flow with high spatial and temporal precision54.
% Current sinks thereby correspond to the activity of excitatory synaptic populations, while current sources mainly reflect balancing return currents.
% The CSD thus provides a functional readout of the cortical microcircuitry function, encompassing a wider, mesoscopic field of view than for instance
% single- or multiunit approaches55. Early current sinks in the auditory cortex are, therefore, indicative of thalamic input in granular layers III/IV
% and infragranular layers Vb/VI21,56–59.
%% interpolate spatially
% Reyes-Puerta et al.: Data were finally interpolated and plotted as pseudocolor images, with current sources and sinks represented
% by red (positive) and blue (negative) colors, respectively.
% Sakata et al.: Finally, data were linearly interpolated and plotted as pseudocolor images, with red (negative) and blue (positive) representing
% current sources and sinks, respectively.
% IMPORTANT: These statements are contradictory! Yang et al. have a factor of -1 in their calculation which is contained in neither paper.
upsample = 10;
csd_i = nan(size(csd_r,1)*upsample,size(csd_r,2));
for i = 1:size(csd_i,2)
  csd_i(:,i) = interp1(1:size(csd_r,1),csd_r(:,i),linspace(1,size(csd_r,1),size(csd_i,1)));
end
%% plot if requested
if doplot
  figure('position',[0.2 0.1 0.6 0.8])
  subplot(221),title('Raw LFP'),hold on
  colors = gray(size(lfpmat,1)+2);
  dist = 0+es/2*1000:1000*es:(size(lfpmat,1)-1)*1000*es+es/2*1000;  % spatial scale for LFP - we will plot each colored square centered on the electrode position
  for i = 1:size(lfpmat,1)
    plot(t,lfpmat(i,:),'Color',colors(i+1,:))
  end
  plot([0,0],get(gca,'YLim'),'k-')
  xlim([t(1),t(end)])
  xlabel('Time (ms)'),ylabel('Voltage (mV)')
  
  % plot smoothed LFP traces without expansion
  ax1 = subplot(222);title('Smoothed LFP w/o expansion'),hold on
	colormap(ax1,flipud(jet))
  imagesc(t,dist,lfp_s)  % imagesc automatically scales the sizes of C to the vectors x and y
  set(gca,'YDir','reverse')
  axis tight,c = colorbar;
  c.Label.String = 'LFP (mV)';
  for i = 1:size(lfpmat,1)
    plot(t,-100*(lfpmat(i,:)-lfpmat(i,1))+dist(i),'k')
  end
  plot([0,0],get(gca,'YLim'),'k-')
  xlim([t(1),t(end)])
  set(gca,'CLim',[-max(abs(get(gca,'CLim'))),max(abs(get(gca,'CLim')))])  % make color code symmetrical around 0
  xlabel('Time (ms)'),ylabel('Distance (µm)')
  
  % plot raw CSD image
  ax2 = subplot(223);title('CSD raw'),hold on
  if theMinus,colormap(ax2,jet)
  else colormap(ax2,flipud(jet)),end
  imagesc(t,dist(2:end-1),csd_r)
  plot([0,0],get(gca,'YLim'),'k-')
  set(gca,'YDir','reverse')
  axis tight,c = colorbar;
  c.Label.String = 'CSD (mV/mm^2)';
  xlim([t(1),t(end)])
  set(gca,'CLim',[-max(abs(get(gca,'CLim'))),max(abs(get(gca,'CLim')))])
  xlabel('Time (ms)'),ylabel('Distance (µm)')
  
  % plot interpolated CSD image
  ax3 = subplot(224);title('CSD smoothed 10x'),hold on
  if theMinus,colormap(ax3,jet)
  else colormap(ax3,flipud(jet)),end
  imagesc(t,1000*(es:es/upsample:es*(size(lfpmat,1)-1)),csd_i)
  plot([0,0],get(gca,'YLim'),'k-')
  set(gca,'YDir','reverse')
  axis tight,c = colorbar;
  c.Label.String = 'CSD (mV/mm^2)';
  xlim([t(1),t(end)])
  set(gca,'CLim',[-max(abs(get(gca,'CLim'))),max(abs(get(gca,'CLim')))])
  xlabel('Time (ms)'),ylabel('Distance (µm)')
  
  disp(' ')
  disp('Remarks on LFP and CSD plots')
  disp(['No. of electrode channels: ',num2str(size(lfpmat,1)),', spaced at ',num2str(es*1000),' µm'])
  disp('LFPs were spatially smoothed, with the top and bottom channels duplicated and then removed.')
  disp(['The raw CSD in csd_r has necessarily two channels (top and bottom) less (',num2str(size(csd_r,1)),') than the LFP input matrix.'])
  disp('The 10x interpolated CSD in csd_i has the same top and bottom channels!')
  disp('All pseudocolor plots have squares plotted centered on the electrode position.')
  disp(' ')
end
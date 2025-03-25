function waveParms = mwave(meanwave,si,varargin)
% waveParms = mwave(meanwave,si,varargin)
%
% mwave computes full width at half maximum (fwhm) for the largest positive and negative peaks of the waveform
% as well as peak-to-trough amplitude, duration, and ratio; these five values are put together in row vector 'waveParms' such that
%
% waveParms = [fwhmMax,fwhmMin,p2tAmp,p2tDur,p2tRat]
%
% ARGUMENTS
% meanwave    required - vector containing the mean spike waveform
% si          required - scalar value indicating sampling interval in microseconds
%
% plot        optional - a figure will be created with the mean wave and highlight measurements; should be third input argument
% handle      optional - to direct the plotting; should be fourth argument and is only valid if the third argument is 'plot'
%
% EXAMPLE
% waveParms = mwave(mean(s(i).spx1ValuesRaw(idx1,:)),10^6*s(i).interval,'plot','handle',h);
% 
% HISTORY
% 2022 April  added optional argument to hand over a figure handle for plotting
% 2020 Oct    changed interpolation method from 'linear' to 'spline'
% 2015 Sep    made sure that always five output values are returned
% 2015 May    added three new outputs, p2tAmp, p2tDur, and p2tRat
%             deleted argument 'pol'
% 2014 April  added new argument 'pol'
% 2014 Jan    replaced function interp with function interp1
%
% by Maik C. Stuettgen, Feb 2013
%% upsample waveform, get p2tAmp and p2tRat
if ~isa(meanwave,'double')
  meanwave = double(meanwave);
end
meanwave = meanwave-meanwave(1);       % pull start of meanwave to 0
upsamp   = 10;
dummy    = interp1(meanwave,linspace(1,length(meanwave),length(meanwave)*upsamp),'spline');
p2tAmp   = max(dummy)-min(dummy);      % peak-to-trough amplitude
p2tRat   = abs(max(dummy)/min(dummy)); % peak-to-trough ratio
%% compute full width at half maximum and peak-to-trough duration
% extract FWHM from upsampled waveform in dummy
[posPeak,posPeakIdx] = max(dummy);
[negPeak,negPeakIdx] = min(dummy);

prePosPeakIdx  = find(dummy(1:posPeakIdx)<posPeak/2,1,'last');
postPosPeakIdx = find(dummy(posPeakIdx:end)<posPeak/2,1,'first') + posPeakIdx-1;
fwhmMax        = si*abs(postPosPeakIdx - prePosPeakIdx)/upsamp;

preNegPeakIdx  = find(dummy(1:negPeakIdx)>negPeak/2,1,'last');
postNegPeakIdx = find(dummy(negPeakIdx:end)>negPeak/2,1,'first') + negPeakIdx-1;
fwhmMin        = si*abs(postNegPeakIdx - preNegPeakIdx)/upsamp;

p2tDur         = si*abs(posPeakIdx - negPeakIdx)/upsamp;

if isempty(fwhmMax),fwhmMax=nan;end
if isempty(fwhmMin),fwhmMin=nan;end
if isempty(p2tAmp),p2tAmp=nan;end
if isempty(p2tDur),p2tDur=nan;end
if isempty(p2tRat),p2tRat=nan;end

waveParms = [fwhmMax,fwhmMin,p2tAmp,p2tDur,p2tRat];
%% plot if requested
if nargin>2
  if strcmp(varargin{1},'plot')
    if nargin>3 && strcmp(varargin{2},'handle')
      subplot(varargin{3})
    else
      figure('units','normalized','position',[.375 .375 .25 .25])
    end
    t = linspace(0,numel(meanwave)*si,numel(dummy));
    interp = plot(t,dummy);hold on
    try
      plot([t(prePosPeakIdx) t(prePosPeakIdx)],[-10^3 10^3],'g:')
      plot([t(postPosPeakIdx) t(postPosPeakIdx)],[-10^3 10^3],'g:')
      plot([0 t(end)],[posPeak/2 posPeak/2],'g:')
      plot([t(prePosPeakIdx) t(postPosPeakIdx)],[posPeak/2 posPeak/2],'g')
      plot([t(preNegPeakIdx) t(preNegPeakIdx)],[-10^3 10^3],'c:')
      plot([t(postNegPeakIdx) t(postNegPeakIdx)],[-10^3 10^3],'c:')
      plot([0 t(end)],[negPeak/2 negPeak/2],'c:')
      plot([t(preNegPeakIdx) t(postNegPeakIdx)],[negPeak/2 negPeak/2],'c')
      
      % superimpose raw waveform
      raw = plot(linspace(0,t(end),numel(meanwave)),meanwave,'r');
      
      axis([0,t(end),-1.1*(max(abs([max(dummy) min(dummy)]))),1.1*(max(abs([max(dummy) min(dummy)])))])
      % scale x-axis from upsampled ticks to ms
      if max(t)<4000  % waveforms less than 4 ms
        set(gca,'XTick',0:500:t(end),'XTickLabel',num2str((0:0.5:t(end)/1000)'))
      else  % waveforms more than 4 ms
        set(gca,'XTick',0:1000:t(end),'XTickLabel',num2str((0:1.0:t(end)/1000)'))
      end
      xlabel('milliseconds'),ylabel('ADC units')
      legend([interp,raw],'interp','raw')
    catch
      disp('Plot could not be constructed.')
    end
  end
end
end
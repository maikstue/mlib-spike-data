function [P1,f] = mspectan(signal,Fs,n,doplot)
% function [P1,f] = mspectan(signal,Fs,n,doplot)
% computes and plots the discrete fourier transform of a signal with n points for a sampling rate Fs
% in essence, this is just a wrapper function for fft
% in case your Matlab license includes the signal processing toolbox, try pspectrum(signal,Fs)
% 
% n             number of points of the discrete fourier transform (DFT)
%               generally, the lower n, the faster the code, but the coarser the resolution
% doplot        flag, plot (1) or not (0)
% 
% 
% HISTORY
% Nov 2024    replaced smooth (requires Curve Fitting Toolbox) with smoothdata
% Oct 2024    added input argument n
% 
% Maik C. Stuettgen, July 19, 2017
%% input checks
if nargin<4
  error('Not enough input arguments provided.')
end
if ~isvector(signal),error('Only vector input accepted'),end

%% the works
L = length(signal);
if n>floor(length(signal)/2),error(['It is not recommended that n>0.5*length(signal) L = ',num2str(L)]),end

z  = fft(signal,n);

% compute the two-sided spectrum |P2|
P2 = abs(z/n);

% then compute the single-sided spectrum |P1| based on |P2| and the even-valued signal length |L|
P1 = P2(1:ceil(n/2)+1);

% P1(1) and P1(end) need not be multiplied by 2 because these amplitudes correspond to the zero and Nyquist frequencies, respectively, 
% and they do not have the complex conjugate pairs in the negative frequencies.
P1(2:end-1) = 2*P1(2:end-1);

% define the frequency domain |f| 
f = Fs*(0:ceil(n/2))/n;

%% plot the single-sided amplitude spectrum |P1|, if requested
if doplot

  figure('name','Single-sided Amplitude Spectrum of X(t)','units','normalized','position',[.4,.3,.2,.6])
  sgtitle(['DFT with ',num2str(n), ' points'])

  subplot(311),title('Linear abscissa'),hold on
  plot(f,P1)
  % plot(f,smooth(P1,ceil(Fs/100)))
  plot(f,smoothdata(P1,'movmean',ceil(Fs/100)))
  xlim([0,Fs/10])
  xlabel('Frequency (Hz)')
  ylabel('|P1(f)|')
  grid on

  subplot(312)
  semilogx(f,P1),hold on
  title('Logarithmic abscissa')
  % semilogx(f,smooth(P1,ceil(Fs/100)))
  semilogx(f,smoothdata(P1,'movmean',ceil(Fs/100)))
  xlabel('Frequency (Hz)')
  ylabel('|P1(f)|')
  xlim([0 Fs/10])
  grid on

  subplot(313),title('Summed over 1-Hz bins'),hold on
  binnedP1 = nan(floor(n/2),1);
  for i = 1:floor(n/2)
    binnedP1(i) = sum(P1(f>(i-1) & f<=i));
  end
  plot(0.5:floor(n/2)-0.5,binnedP1)
  xlabel('Frequency (Hz)')
  ylabel('|P1(f)|')
  xlim([0 Fs/10])
  grid on
end
end
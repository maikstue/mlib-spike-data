function [mat,optmax,optmin,aurocs,cis] = mwa2(tspx,tevents1,tevents2,winpos,winsize,varargin)
% function [mat,optmax,optmin,aurocs,cis] = mwa2(tspx,tevents,tevents2,winpos,winsize,varargin)
%
% code to compute the area under the ROC curve for two spike count distributions at different time points
% and different time intervals
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% IMPORTANT: UNLIKE MWA.M, THIS CODE SHIFTS TIME WINDOWS FOR BOTH EVENTS!
% GOOD FOR DISCRIMINATION, WHERE PSTHS FOR TWO STIMULI ARE TO BE COMPARED!
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% note that function requires mroc and mnspx in Matlab path
%
% for an example how to use mwa.m, confer Figure 3 and associated manuscript text in:
% Stuettgen MC, Schwarz C (2008) Psychophysical and neurometric detection performance under stimulus uncertainty.
% Nature Neuroscience 11(9): 1091-1098.
%
% MANDATORY INPUT ARGUMENTS
% tspx            vector of spike timestamps (in seconds)
% tevents1        vector of timestamps for event 1 (in seconds)
% tevents2        vector of timestamps for event 2 (in seconds)
% winpos          vector specifying window positions from start to end [tstart:tshift:tend],
%                 where tstart is starting time, tshift is the window displacement interval,
%                 and tend is the position of the last window (all times in seconds)
% winsize         vector specifying which window sizes to use [smalles:tshift:largest] (in seconds)
%
% OPTIONAL INPUT ARGUMENT
% smoothit        if set to 1, winsize and winpos are expanded by two rows and two columns, the matrix is then smoothed,
%                 and minima and maxima are computed on the basis of the smoothed matrix
% plotit          if set to 1 or 2 or 3, generates a figure with AUROC values color-coded (1: grayscale, 2: autumn, 3: fancy)
% after0          if set to 1, maximum and minimum aurocs are only found in time windows centered at values >=0; the default is 1
% ci              confidence interval width, should be a number >=50 and <100 (percentage); if argument is specified, CIs will be bootstrapped
%                 For each window size, take time points from 0 to the end of the session (in half-second steps). Take samples
%                 of the size of tevents (with replacement).
% nboot           number of bootstrap iterations for CIs
%
% OUTPUT
% mat             2D-matrix of AUROC values with rows denoting window sizes and columns window positions (as in the figure)
% optmax          2-element vector specifying which combination of window position (1) and size (2) yields the maximum AUROC value
% optmin          2-element vector specifying which combination of window position (1) and size (2) yields the minimum AUROC value
% aurocs          2-element vector specifying the maximum (1) and minimum (2) auroc values
% cis             matrix with two columns; each row denotes a window size, columns give lower and upper auroc bounds
%
% HISTORY
% Oct 2024        polished code
% July 2018       added optional confidence interval calculations
% Mar 2018        changed order of winpos and winsize as input arguments, bug correction for optimal combinations and auroc
%                 added smoothit as optional input argument
% mwa2.m by Maik C. Stuettgen, October 2016
% mwa.m by Maik C. Stuettgen, Summer 2013 @ Erasmus MC Rotterdam, The Netherlands
%% inputcheck
if ~isvector(tspx) || ~isvector(tevents1) || ~isvector(tevents2) || ~isnumeric(winsize) || ~isnumeric(winpos)
  error('wrong input')
end
plotit   = 0;
smoothit = 0;
after0   = 1;
ci       = 0;
nboot    = 100;
winsize  = winsize*1000;  % because mnspx works with ms - better to scale once here than multiple times below
if ~isempty(varargin)
  for i = 1:2:size(varargin,2)
    switch varargin{i}
      case 'plotit'
        plotit = varargin{i+1};
      case 'smoothit'
        smoothit = varargin{i+1};
      case 'after0'
        after0 = varargin{i+1};
      case 'ci'
        ci = varargin{i+1};
      case 'nboot'
        nboot = varargin{i+1};
    end
  end
end
%% bootstrap confidence intervals, if desired
if ci~=0
  auroc4ci_mat = nan(numel(winsize),1);  % preallocate for speed
  % for every winsize, compute the reference distribution of spike times (refDist)
  % then compute AUROC for a set of randomly selected spike counts from the session
  for i = 1:numel(winsize)
    refDist = mnspx(tspx,0.5:0.5:max([tevents1(end),tevents2(end)]),winsize(i)/2,winsize(i)/2);
    for j = 1:nboot
      auroc4ci_mat(i,j) = mroc(randsample(refDist,numel(tevents1),true),randsample(refDist,numel(tevents2),true));
    end
    clear refDist
  end
  cis = prctile(auroc4ci_mat,[(100-ci)/2,100-(100-ci)/2],2);
  clear auroc4ci_mat i j
else
  cis = [];
end
%% if smoothing is desired, adjust winpos and winsize
if smoothit==1
  diffwinpos  = [winpos(2)-winpos(1),winpos(end)-winpos(end-1)];
  diffwinsize = [winsize(2)-winsize(1),winsize(end)-winsize(end-1)];
  winpos  = [winpos(1)-diffwinpos(1),winpos,winpos(end)+diffwinpos(2)];
  winsize = [winsize(1)-diffwinsize(1),winsize,winsize(end)+diffwinsize(2)];
  if winsize(1)<=0
    winsize(1) = winsize(2)-diffwinsize(1)/2;
  end
  clear diffwinpos diffwinsize
end
%% the works
mat = nan(numel(winsize),numel(winpos));  % preallocate for speed
% for every winsize, compute auroc value for all winpos
for i = 1:numel(winsize)
  for j = 1:numel(winpos)
    mat(i,j) = mroc(mnspx(tspx,tevents1+winpos(j),winsize(i)/2,winsize(i)/2),...
                    mnspx(tspx,tevents2+winpos(j),winsize(i)/2,winsize(i)/2));
  end
end
%% smooth the matrix?
if smoothit==1
  h = 1/9*ones(3,3);  % must add up to 1
  mat = conv2(mat,h,'valid');
  winpos = winpos(2:end-1);
  winsize = winsize(2:end-1);
end
%% determine optimal winsize-winpos combinations; issue warning when n(max)>1 | n(min)>1
if ~after0
  % use the column representation of mat, mat(:), to find  value and index of the smallest element
  [aurocs(1),I] = max(mat(:));
  % The ind2sub function determines the equivalent subscript values corresponding to a single index into an array.
  [Ioptwinsize,Ioptwinpos] = ind2sub(size(mat),I);
  optmax = [winpos(Ioptwinpos),winsize(Ioptwinsize)/1000];
  clear I Ioptinwpos Ioptwinsize

  [aurocs(2),I] = min(mat(:));
  [Ioptwinsize,Ioptwinpos] = ind2sub(size(mat),I);
  optmin = [winpos(Ioptwinpos),winsize(Ioptwinsize)/1000];
  clear I Ioptinwpos Ioptwinsize
  if sum(mat(:)>=aurocs(1))>1
    disp(['MWA2 WARNING: ',num2str(sum(mat(:)>=aurocs(1))),' windows contain the maximum auroc value.'])
  end
  if sum(mat(:)<=aurocs(2))>1
    disp(['MWA2 WARNING: ',num2str(sum(mat(:)<=aurocs(2))),' windows contain the minimum auroc value.'])
  end
else
  matdummy = mat(:,find(winpos>=0,1,'first'):end);
  winposdummy = winpos(find(winpos>=0,1,'first'):end);
  [aurocs(1),I] = max(matdummy(:));
  [Ioptwinsize,Ioptwinpos] = ind2sub(size(matdummy),I);
  optmax = [winposdummy(Ioptwinpos),winsize(Ioptwinsize)/1000];
  clear I Ioptinwpos Ioptwinsize

  [aurocs(2),I] = min(matdummy(:));
  [Ioptwinsize,Ioptwinpos] = ind2sub(size(matdummy),I);
  optmin = [winposdummy(Ioptwinpos),winsize(Ioptwinsize)/1000];
  if sum(matdummy(:)>=aurocs(1))>1
    disp(['MWA2 WARNING: ',num2str(sum(matdummy(:)>=aurocs(1))),' windows contain the maximum auroc value.'])
  end
  if sum(matdummy(:)<=aurocs(2))>1
    disp(['MWA2 WARNING: ',num2str(sum(matdummy(:)<=aurocs(2))),' windows contain the minimum auroc value.'])
  end
  clear I Ioptinwpos Ioptwinsize matdummy
end
%% plot if desired
if ismember(plotit,1:3)
  figure('units','normalized','position',[.4 .5 .2 .2]),hold on
  winsize = winsize/1000;   % scale winsize back to original values for plotting
  switch plotit
    case 1
      imagesc([winpos(1),winpos(end)],[winsize(1),winsize(end)],mat)   % imagesc(x,y,C) specifies the image location. Use x and y to specify the locations of the corners 
      colormap(gray)                  % corresponding to C(1,1) and C(m,n). To specify both corners, set x and y as two-element vectors.
                                      %  To specify the first corner and let imagesc determine the other, set x and y as scalar values.
                                      % The image is stretched and oriented as applicable.
    case 2
      imagesc([winpos(1),winpos(end)],[winsize(1),winsize(end)],mat)
      colormap(autumn)
    case 3
      if numel(winsize)>1
        pcolor(winpos,winsize,mat)
      else
        imagesc(winpos,winsize,mat)
      end
        colormap(hot)
        shading interp
  end
  maxval = max(max(mat));
  minval = min(min(mat));
  r      = ceil(max([abs(0.5-maxval) abs(0.5-minval)])*10)/10;
  set(gca,'YDir','normal','Clim',[0.5-r 0.5+r])
  colorbar
  axis tight
  xlabel('window position (seconds)')
  ylabel('window size (seconds)')
  ylim([min(winsize)-min(winsize)/2,max(winsize)+min(winsize)/2])
  plot([0 0],[0 ceil(max(winsize))],'k-')
  plot(optmax(1),optmax(2),'Marker','*','MarkerEdgeColor','k')
  plot(optmin(1),optmin(2),'Marker','d','MarkerEdgeColor','k')
end
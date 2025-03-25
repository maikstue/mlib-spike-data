function [mat,optmax,optmin,aurocs,cis] = mwa(tspx,tevents,trefevents,winpos,winsize,varargin)
% function [mat,optmax,optmin,aurocs,cis] = mwa(tspx,tevents,trefevents,winpos,winsize,varargin)
%
% computes the area under the ROC curve for two spike count distributions at different time points
% and different time intervals
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% IMPORTANT: THIS CODE IS TO BE USED FOR DETECTION TASKS ONLY!
% Windows are only shifted for tevents, not for trefevents!
% For a code in which time windows shift for both events, see mwa2.m.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% note that function requires mroc and mnspx in Matlab path
%
% for an example, see Figure 3 in:
% Stüttgen MC, Schwarz C (2008) Psychophysical and neurometric detection performance under stimulus uncertainty.
% Nature Neuroscience 11(9): 1091-1098.
%
% MANDATORY INPUT ARGUMENTS
% tspx            vector of spike timestamps (in seconds)
% tevents         vector of event timestamps (in seconds)
% trefevents      vector of reference event timestamps (in seconds)
% winpos          vector specifying window positions from start to end [tstart:tshift:tend],
%                 where tstart is starting time, tshift is the window displacement interval,
%                 and tend is the position of the last window (all times in seconds)
% winsize         vector specifying which window sizes to use [smallest:tshift:largest] (in seconds)
%
% OPTIONAL INPUT ARGUMENTS
% smoothit        if set to 1, winsize and winpos are expanded by two rows and two columns per value, the filter kernel will
%                 be 3-by-3, the matrix is then smoothed, and minima and maxima are computed on the basis of the smoothed matrix
%                 the default value is 0 (no smoothing)
% plotit          if set to 1 or 2 or 3, generates a figure with AUROC values color-coded (1: grayscale, 2: autumn, 3: fancy)
% after0          if set to 1, maximum and minimum aurocs are only found in time windows centered at values >=0; the default is 1
% ci              two-element vector; first vector specifies confidence interval width in percent, should be a number >=50 and <100 
%                 if argument is specified, CIs will be bootstrapped with nboot iterations (2nd number)
%                 For each window size, take time points from 0 to the end of the session (in half-second steps). Take samples
%                 of the size of tevents (with replacement).
% nboot           number of bootstrap iterations for CIs (default is 100)
%
% OUTPUT
% mat             2D-matrix of AUROC values with rows denoting window sizes and columns window positions (as in the figure)
% optmax          2-element vector specifying which combination of window position (1) and size (2) yields the maximum AUROC value
% optmin          2-element vector specifying which combination of window position (1) and size (2) yields the minimum AUROC value
% aurocs          2-element vector specifying the maximum (1) and minimum (2) auroc values
% cis             matrix with two columns; each row denotes a window size, columns give lower and upper auroc bounds
% 
% HISTORY
% Oct 2024        polished code, set after0 default to 1
% July 2018       added optional confidence interval calculations
% April 2018      removed optional smoothing for 5-by-5-kernel; if winsize(1) for smoothing is <=0, then it will be adjusted
% Mar 2018        changed order of winpos and winsize as input arguments, bug correction for optimal combinations and aurocs
%                 added smoothit as optional argument
% Oct 2016        removed minor bugs (optimum plotted incorrectly, x and y axes were swapped; y-axis did not scale well)
%                 output vectors mat and opt are now transposed; added output auroc
%                 mwa now returns both highest and lowest AUROC values, along with their respective size-position combination
% Oct 2014        added if-structure in plotting command such that the fancy plotting option also works when only a single winsize is specified
% 
% Maik C. Stuettgen, Summer 2013 @ Erasmus MC Rotterdam, The Netherlands
%% inputcheck
if ~isvector(tspx) || ~isvector(tevents) || ~isvector(trefevents) || ~isnumeric(winsize) || ~isnumeric(winpos)
  error('wrong input')
end
plotit   = 0;
smoothit = 0;
after0   = 1;
ci       = [0,100];
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
    end
  end
end
%% bootstrap confidence intervals, if requested
if ci(1)~=0
  auroc4ci_mat = nan(numel(winsize),1);  % preallocate for speed
  % for every winsize, get a spike count distribution for many time points
  % then compute AUROCs by drawing two samples with replacement from the distributions
  for i = 1:numel(winsize)
    refDist  = mnspx(tspx,0.5:0.5:tevents(end),winsize(i)/2,winsize(i)/2);  % many events over the entire session, spaced at half-second intervals
%     refDist  = mnspx(tspx,tevents-max(winsize)/2,winsize(i)/2,winsize(i)/2);  % time points closely preceding tevents
    for j = 1:ci(2)
      auroc4ci_mat(i,j) = mroc(randsample(refDist,numel(tevents),true),randsample(refDist,numel(tevents),true));
    end
    clear refDist
  end
  cis = prctile(auroc4ci_mat,[(100-ci(1))/2,100-(100-ci(1))/2],2);
  clear auroc4ci_mat i j
else
  cis = [];
end
%% if smoothing is requested, adjust winpos and winsize
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
% for every winsize, compute the reference distribution (compDist)
% then compute AUROC for every winpos for the given winsize
for i = 1:numel(winsize)
  compDist = mnspx(tspx,trefevents-winsize(i)/2,winsize(i)/2,winsize(i)/2); % spike count distribution for comparison - need compute only once for each winsize
  for j = 1:numel(winpos)
    mat(i,j) = mroc(mnspx(tspx,tevents+winpos(j),winsize(i)/2,winsize(i)/2),compDist);
  end
end
%% smooth the matrix, if requested
if smoothit==1
  h = 1/9*ones(3,3);  % must add up to 1
  mat = conv2(mat,h,'valid');
  winpos = winpos(2:end-1);
  winsize = winsize(2:end-1);
end
%% determine optimal winpos-winsize combinations; issue warning when n(max)>1 | n(min)>1
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
    disp(['MWA WARNING: ',num2str(sum(mat(:)>=aurocs(1))),' windows contain the maximum auroc value - output set to NaN.'])
    optmax = [nan,nan];optmin = [nan,nan];aurocs = [nan,nan];
  end
  if sum(mat(:)<=aurocs(2))>1
    disp(['MWA WARNING: ',num2str(sum(mat(:)<=aurocs(2))),' windows contain the minimum auroc value - output set to NaN.'])
    optmax = [nan,nan];optmin = [nan,nan];aurocs = [nan,nan];
  end
else
  matdummy    = mat(:,find(winpos>=0,1,'first'):end);
  winposdummy = winpos(find(winpos>=0,1,'first'):end);
  [aurocs(1),I] = max(matdummy(:));
  [Ioptwinsize,Ioptwinpos] = ind2sub(size(matdummy),I);
  optmax = [winposdummy(Ioptwinpos),winsize(Ioptwinsize)/1000];
  clear I Ioptinwpos Ioptwinsize
  
  [aurocs(2),I] = min(matdummy(:));
  [Ioptwinsize,Ioptwinpos] = ind2sub(size(matdummy),I);
  optmin = [winposdummy(Ioptwinpos),winsize(Ioptwinsize)/1000];
  if sum(matdummy(:)>=aurocs(1))>1
    disp(['MWA WARNING: ',num2str(sum(matdummy(:)>=aurocs(1))),' windows contain the maximum auroc value - output set to NaN.'])
    optmax = [nan,nan];optmin = [nan,nan];aurocs = [nan,nan];
  end
  if sum(matdummy(:)<=aurocs(2))>1
    disp(['MWA WARNING: ',num2str(sum(matdummy(:)<=aurocs(2))),' windows contain the minimum auroc value - output set to NaN.'])
    optmax = [nan,nan];optmin = [nan,nan];aurocs = [nan,nan];
  end
  clear I Ioptinwpos Ioptwinsize matdummy winposdummy
end
%% plot if requested
% imagesc(x,y,C) specifies the image location. Use x and y to specify the locations of the corners corresponding to C(1,1) and C(m,n).
% To specify both corners, set x and y as two-element vectors.
%  To specify the first corner and let imagesc determine the other, set x and y as scalar values.
% The image is stretched and oriented as applicable.
if ismember(plotit,1:3)
  figure('units','normalized','position',[.4 .5 .2 .2]),hold on
  winsize = winsize/1000;   % scale winsize back to original values for plotting
  switch plotit
    case 1
      imagesc(mat)      % imagesc([winpos(1),winpos(end)],[winsize(1),winsize(end)],mat)
      colormap(gray)
    case 2
      imagesc(mat)
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

  % choose a fitting number of xticks
  if numel(winpos)<=20,xstepsz = 2;
  elseif numel(winpos)<=30,xstepsz = 3;
  elseif numel(winpos)<=40,xstepsz = 4;
  elseif numel(winpos)<=60,xstepsz = 5;
  else xstepsz = 10;
  end
  my0 = find(winpos==0);
  if ~isempty(my0)
    plot([my0 my0],[0 numel(winsize)+0.5],'k-')
    set(gca,'XTick',[fliplr(my0:-xstepsz:1),my0+xstepsz:xstepsz:numel(winpos)],'XTickLabel',num2str(winpos([fliplr(my0:-xstepsz:1),my0+xstepsz:xstepsz:numel(winpos)])'))
  else
    set(gca,'XTick',1:xstepsz:numel(winpos),'XTickLabel',num2str(winpos([1:xstepsz:numel(winpos)])'))
  end

  % choose a fitting number of yticks
  if numel(winsize)<=8
    set(gca,'YTick',1:numel(winsize),'YTickLabel',num2str(winsize'))
  elseif numel(winsize)>8 & numel(winsize)<=20
    set(gca,'YTick',1:2:numel(winsize),'YTickLabel',num2str(winsize(1:2:end)'))
  elseif numel(winsize)>20
    set(gca,'YTick',1:4:numel(winsize),'YTickLabel',num2str(winsize(1:4:end)'))
  end

  % final additions and modifications
  plot(find(winpos==optmax(1)),find(winsize==optmax(2)),'Marker','*','MarkerEdgeColor','k')
  plot(find(winpos==optmin(1)),find(winsize==optmin(2)),'Marker','o','MarkerEdgeColor','k')
  ylim([0.5,numel(winsize)+0.5])
end
end
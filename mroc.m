function [auroc,bp,ciLoUp] = mroc(x,y,varargin)
% function [auroc,bp,ciLoUp] = mroc(x,y,varargin)
%
% computes the area under the receiver operating characteristic (ROC) curve
% optionally, provides bootstrapped confidence intervals and plots the curve
%
% note that mroc assumes that x and y are arrays of integers
% if input contains non-integer rational numbers, all numbers are is converted to ranks
%
% results validated a) by simulation with normcdf/norminv b) by comparison with web-based calculators and 
% c) with code in the MES Toolbox
%
% INPUT ARGUMENTS
% x         input vector ('signal')
% y         input vector ('noise')
%           x and y are processed such that auroc>0.5 implies that x holds larger values and vice versa
%           so, if x is a spike count distribution during stimulation and y is a spike count distribution during baseline,
%           auroc>0.5 implies excitation and auroc <0.5 implies inhibition
%
% OPTIONAL INPUT ARGUMENTS
% plotit    if set to 1, generates a figure of the receiver operation characteristic curve
% ci        2-element vector; first element should be a scalar x 0<x<1 which determines the extent of the confidence interval;
%           second element should be an integer (recommended: 1000) specifying the number of bootstraps
%
% OUTPUT ARGUMENTS
% auroc     area under the ROC curve; equals the probability that a random sample taken from distribution x
%           is greater than a random sample taken from distribution y
% bp        bootstrapped p-value (2-sided)
% ciLoUp    2-element vector holding lower and upper bootstrapped confidence interval limits
%
% HISTORY
% November 2025   moved the if-case for kicking out NaNs to the top to avoid problems with the if-case about integers
% October 2024    added bootstrapped p-value to output arguments; polished code here and there; added histogram plot
% October 2021    deleted variable 'resol'
%                 added automatic conversion of non-integer values into ranks
% 
% Maik C. Stüttgen, Oct 2014
%% check input
% do we have NaNs? kick them out
if any(isnan(x)) || any(isnan(y))
  disp(['[',8,'Kicking out NaNs]',8]) % in orange
  x(isnan(x)) = [];
  y(isnan(y)) = [];
end

% do x and y really hold integers and integers only?
if (~all(x==int32(x)) || ~all(y==int32(y)) || ~isvector(x) || ~isvector(y))
  disp(['[',8,'Non-integer values provided, converting into ranks]',8]) % in orange
  z = [x;y];z = tiedrank(z);x = z(1:numel(x));y = z(numel(x)+1:end);
elseif isempty(x) || isempty(y)
  auroc  = nan;
  bp     = nan;
  ciLoUp = nan;
  return
end

%% preparation
plotit = [];
ci     = [];
ciLoUp = [];
bp     = [];
if ~isempty(varargin)
  for i = 1:2:size(varargin,2)
    switch varargin{i}
      case 'plotit'
        plotit = varargin{i+1};
      case 'ci'
        ci     = varargin{i+1};
    end
  end
end

% convert to column vectors
x = x(:);y = y(:);

%% the works
start_val = min([x;y])-1;    % overall smallest value -1
end_val   = max([x;y])+1;    % overall largest value +1
c         = start_val:end_val;         % list of criterion values
fahr      = zeros(numel(c),2);         % initialize HR and FA array
for i = 1:size(c,2)                    % loop through all criterion values
  fahr(i,1) = sum(y>c(i))/length(y);      % get p(y>c), i.e. false alarm rate
  fahr(i,2) = sum(x>c(i))/length(x);      % get p(x>c), i.e. hit rate
end
fahr     = [fahr;0,0;1,0;1,1];               % close the curve
auroc    = polyarea(fahr(:,1),fahr(:,2));    % integrate

%% compute bootstrapped p-value, if requested
if nargout>1
  auroc4bp = nan(ci(2),1);
  for i = 1:ci(2)
    a = zeros(numel(c),2);                   % initialize HR and FA array
    x1 = randsample([x;y],numel(x),true);
    y1 = randsample([x;y],numel(y),true);
    for j = 1:size(c,2)                      % loop through all criterion values
      a(j,1) = sum(y1>c(j))/length(y1);      % get p(y>c)
      a(j,2) = sum(x1>c(j))/length(x1);      % get p(x>c)
    end
    a = [a;0,0;1,0;1,1];                           % close the curve
    auroc4bp(i,1) = polyarea(a(:,1),a(:,2));       % integrate
  end
  if auroc>0.5,bp = 2*sum(auroc4bp>auroc)/ci(2);
  elseif auroc<0.5,bp = 2*sum(auroc4bp<auroc)/ci(2);
  else bp = 1;
  end
end

%% compute bootstrapped confidence interval, if requested
if nargout>2
  auroc4ci = nan(ci(2),1);
  for i = 1:ci(2)
    a = zeros(numel(c),2);                   % initialize HR and FA array
    x1 = randsample(x,numel(x),true);
    y1 = randsample(y,numel(y),true);
    for j = 1:size(c,2)                      % loop through all criterion values
      a(j,1) = sum(y1>c(j))/length(y1);      % get p(y>c)
      a(j,2) = sum(x1>c(j))/length(x1);      % get p(x>c)
    end
    a = [a;0,0;1,0;1,1];                           % close the curve
    auroc4ci(i,1) = polyarea(a(:,1),a(:,2));       % integrate
  end
  ciLoUp = prctile(auroc4ci,100*[(1-ci(1))/2 1-(1-ci(1))/2]);
end

%% plot if requested
if plotit
  figure('units','normalized','position',[.3 .5 .4 .3])
  
  subplot(121),title(['AUROC=' num2str(auroc,'%1.2f')]),hold on
  scatter(fahr(:,1),fahr(:,2),'.')
  plot(fahr(:,1),fahr(:,2))
  plot([0 1],[0 1],'k')
  xlabel('p(y>c)'),ylabel('p(x>c)')
  axis([0 1 0 1])
  
  subplot(122)
  histogram(x,-0.5:max([x;y]+0.5)),hold on
  histogram(y,-0.5:max([x;y]+0.5))
  xlabel('Counts')
  ylabel('Frequency')
  legend('x','y','Location','best')
end

end
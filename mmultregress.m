function [b,p,sr,psr] = mmultregress(tspx,tevents,predictors,winpos,winsize,varargin)
% function [b,p,sr,psr] = mmultregress(tspx,tevents,predictors,winpos,winsize,varargin)
% 
% Performs multiple regression analysis on spike counts with moving windows of specified durations and at multiple positions.
% If interactions between predictors are to be included, they should be specified within the predictors input matrix.
% Additionally computes semipartial correlation coefficients to indicate the unique contribution of each predictor.
% 
% INPUT ARGUMENTS
% tspx            vector of spike timestamps (in seconds)
% tevents         vector of event timestamps (in seconds)
% predictors      matrix of predictors (integers); a column of 1s to estimate the intercept is added automatically
%                 the code only works with 2 or 3 predictor variables
% winpos          vector specifying window positions from start to end [tstart:tshift:tend],
%                 where tstart is starting time, tshift is the window displacement interval,
%                 and tend is the position of the last window (all times in seconds)
% winsize         vector specifying which window sizes to use [smalles:tshift:largest] (in seconds)
%                 for example, if winpos=0 and winsize=0.1, the spike counting window will extent from -0.05 to +0.05 ms re: tevents
%
% OPTIONAL INPUT ARGUMENT
% udist           user-specified distribution, can be 'normal' (default) or 'poisson' (for spike counts)
% prednames       predictor names; cell array with strings, e.g. {'stimulus','response','reward'}
%                 prednames should have the same number of elements as predictors has columns; omit the intercept
% plotit          if present, generates a figure
% normalize       if present, normalizes criterion variable
% 
% OUTPUT ARGUMENTS
% b               matrix of regression coefficients, dimensions are winsize*winpos*number of regressors; the first regressor returned is always the intercept
% p               matrix of p-values corresponding to the regression coefficients in
% sr              matrix of semipartial correlation coefficients
% psr             matrix of p-values corresponding to the semipartial correlation coefficients
% 
% EXAMPLE
% load ratUnitForMLIBTesting
% 
% trialsForAnalysis = 1:spx.block1;
% dataMatrix = spx.allBData(trialsForAnalysis,:);     % analyze only the first 120 trials
% dataMatrix(isnan(dataMatrix(:,5)),:) = [];           % delete all trials without a response
% 
% % extract relevant information and recode to -1 and +1
% stimuli  = dataMatrix(:,3) - 1;
% choices  = dataMatrix(:,5) - 1;
% outcomes = dataMatrix(:,6);
% toutcomes = dataMatrix(:,7);    % time stamps (in seconds) of when the time entered the side port; reward was delivered (or not) after 0.5 s
% 
% predictors = [stimuli,choices,outcomes];
% predictors(predictors==0) = -1;
% 
% winpos  = -1:0.1:1.5;    % window positions (in s)
% winsize = [0.2,0.5];     % window durations (in s)
% 
% % [b,p,sr,psr] = mmultregress(tspx,tevents,predictors,winpos,winsize,varargin)
% [b,p,sr,psr] = mmultregress(spx.tspx,toutcomes,predictors,winpos,winsize);
% [b,p,sr,psr] = mmultregress(spx.tspx,toutcomes,predictors,winpos,winsize,'plotit');
% [b,p,sr,psr] = mmultregress(spx.tspx,toutcomes,predictors,winpos,winsize,'plotit','prednames',{'stimulus','choice','outcome'});
% [b,p,sr,psr] = mmultregress(spx.tspx,toutcomes,predictors,winpos,winsize,'plotit','prednames',{'stimulus','choice','outcome'},'normalize');
% [b,p,sr,psr] = mmultregress(spx.tspx,toutcomes,predictors,winpos,winsize,'plotit','prednames',{'stimulus','choice','outcome'},'udist','poisson');
% 
% 
% Jan 2024 by Maik C. Stüttgen, University Medical Center Mainz, Germany
%% input check
npredictors = size(predictors,2);
if npredictors>6,error('Only up to six predictors allowed.'),end
if ~isvector(tspx) || ~isvector(tevents) || ~isnumeric(winsize) || ~isnumeric(winpos),error('Wrong input'),end
if rank(predictors)~=npredictors,error(['Predictor matrix is rank-deficient (rank ',num2str(rank(predictors)),', should be ',num2str(npredictors),').']),end

plotit      = false;
normalize   = false;
udist       = 'normal';

if ~isempty(varargin)
  for i = 1:size(varargin,2)
    if ~iscell(varargin{i})
      switch varargin{i}
        case 'udist',udist  = varargin{i+1};
        case 'prednames',prednames = varargin{i+1};
          if npredictors~=numel(prednames),error('Incorrect number of predictor names specified.'),end
        case 'plotit',plotit = true;
        case 'normalize',normalize = true;
      end
    end
  end
end

if ~exist('prednames','var')  % for convenience, assign predictor names if none are specified
  for i = 1:npredictors
    prednames{i} = char(i+96);
  end
end

if (normalize && ~strcmp(udist,'normal'))
  error('Normalization does not work when Poisson or binomial response distributions are specified.')
end

%% preallocation
b = nan(numel(winsize),numel(winpos),npredictors+1);
p = nan(size(b));
sr  = nan(numel(winsize),numel(winpos),npredictors);
psr = nan(size(sr));
if npredictors>3
  disp('Semipartial correlations implemented only for up to three predictors.')
  sr = [];psr = [];
end
lm = cell(numel(winsize),numel(winpos));

%% the works
for i = 1:numel(winsize)                  % for every winsize, compute beta values for all winpos
  for j = 1:numel(winpos)

    nspx = mnspx(tspx,tevents+winpos(j),1000*winsize(i)/2,1000*winsize(i)/2);

    if normalize,nspx = zscore(nspx,0);end    % normalize if requested; flag of 0 in zscore indicates the sample variance is to be used

    % multiple regression can be performed with regress, fitlm, fitglm, and glmfit
    % fitlm does not allow to specify a distribution for the response variable, unlike fitglm and glmfit
    % [b,bint,r,rint,stats] = regress(y,X); % stats are r^2, F, p, error variance; most simple and fastest way
    % however; regress returns only a single p-value for the entire regression, but we want to have p-values for every predictor, and the other functions do that
    % for the case of 3 predictors, we compare glmfit and fitglm (the following two commands give identical results):
    % [b(i,j,:),~,stats] = glmfit(predictors,nspx,'Poisson');
    % lm{i,j}  = fitglm(tbl,['nspx ~ ',prednames{1},'+',prednames{2},'+',prednames{3}],'Distribution','poisson');
    % glmfit does not work with tables, and we don't get a model object as results, but it returns b directly and p in stats.p
    % also, glmfit has no interaction unless it is explicitly entered as predictor (which it is in this case)
    % accordingly, we will use fitglm

    switch npredictors
      case 2
        tbl = table(nspx,predictors(:,1),predictors(:,2),'VariableNames',{'nspx',prednames{1},prednames{2}});
        lm{i,j}  = fitglm(tbl,['nspx ~ ',prednames{1},'+',prednames{2}],'Distribution',udist);

        % semipartial correlations
        [~,~,residuals1_2] = regress(predictors(:,1),[ones(size(predictors,1),1),predictors(:,2)]);
        [~,~,residuals2_1] = regress(predictors(:,2),[ones(size(predictors,1),1),predictors(:,1)]);
        [sr1_2,psr1_2] = corr(nspx,residuals1_2);
        [sr2_1,psr2_1] = corr(nspx,residuals2_1);
        sr(i,j,:)  = [sr1_2,sr2_1];
        psr(i,j,:) = [psr1_2,psr2_1];

      case 3
        tbl = table(nspx,predictors(:,1),predictors(:,2),predictors(:,3),'VariableNames',{'nspx',prednames{1},prednames{2},prednames{3}});
        % lm{i,j}  = fitlm(tbl,['nspx ~ ',prednames{1},'+',prednames{2},'+',prednames{3}]);
        lm{i,j}  = fitglm(tbl,['nspx ~ ',prednames{1},'+',prednames{2},'+',prednames{3}],'Distribution',udist); % this gives the same result as fitlm for a normal distribution

        % semipartial correlations
        [~,~,residuals1_23] = regress(predictors(:,1),[ones(size(predictors,1),1),predictors(:,2:3)]);
        [~,~,residuals2_13] = regress(predictors(:,2),[ones(size(predictors,1),1),predictors(:,[1,3])]);
        [~,~,residuals3_12] = regress(predictors(:,3),[ones(size(predictors,1),1),predictors(:,1:2)]);
        [sr1_23,psr1_23] = corr(nspx,residuals1_23);
        [sr2_13,psr2_13] = corr(nspx,residuals2_13);
        [sr3_12,psr3_12] = corr(nspx,residuals3_12);
        sr(i,j,:)  = [sr1_23,sr2_13,sr3_12];
        psr(i,j,:) = [psr1_23,psr2_13,psr3_12];

    end
    b(i,j,:) = table2array(lm{i,j}.Coefficients(:,1));  % columns are regressor estimates, SE, tStat, pValue
    p(i,j,:) = table2array(lm{i,j}.Coefficients(:,4));  % columns are regressor estimates, SE, tStat, pValue
  end
end

%% plotting
if plotit
  prednames{numel(prednames)+1} = 'intercept';
  figure
  for i = 1:size(b,1)   % for all winsizes

    yminmaxtoprow = [0,0];
    yminmaxbottomrow = [0,0];
    if numel(winsize)<6
      subplot(2,numel(winsize),i),title(['Winsize ',num2str(winsize(i)),' ms']),hold on
      plot(winpos,squeeze(b(i,:,2:end)))
      plot(winpos,b(i,:,1))
      if min(b(i,:,:),[],'all')<yminmaxtoprow(1),yminmaxtoprow(1) = min(b(i,:,:),[],'all');end
      if max(b(i,:,:),[],'all')>yminmaxtoprow(2),yminmaxtoprow(2) = max(b(i,:,:),[],'all');end
      if i==1,ylabel('b coefficients'),end
      if i==size(b,1),legend(prednames,'Location','best'),end

      subplot(2,numel(winsize),i+numel(winsize)),hold on
      plot(winpos,squeeze(sr(i,:,:).^2))
      if min(sr(i,:,:).^2,[],'all')<yminmaxbottomrow(1),yminmaxbottomrow(1) = min(sr(i,:,:).^2,[],'all');end
      if max(sr(i,:,:).^2,[],'all')>yminmaxbottomrow(2),yminmaxbottomrow(2) = max(sr(i,:,:).^2,[],'all');end
      if i==1,xlabel('Window position (s)'),ylabel('Squared semipartial correlation'),end
    else
      disp('Sorry, a figure is generated only when numel(winsize)<6.')
    end
  end

  for i = 1:size(b,1),subplot(2,numel(winsize),i),axis([winpos(1),winpos(end),yminmaxtoprow]),end
  for i = 1:size(b,1),subplot(2,numel(winsize),numel(winsize)+i),axis([winpos(1),winpos(end),yminmaxbottomrow]),end
end

end

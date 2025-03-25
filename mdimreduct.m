function [score,loadings] = mdimreduct(smat4dimred,standardize,plotresults)
% function [score,loadings] = mdimreduct(smat4dimred,standardize,plotresults);
% Dimensionality reduction of a spike-count matrix via PCA.
% 
% INPUT
% smat4dimred   data matrix where rows are "variables" and columns are "observations"
%               for spike data, rows should be units/neurons, and columns are time windows in which spikes were counted
%               columns can be either vectors of spike counts across trials (same trial epoch for different consecutive trials) or
%               columns can be vectors of spike counts within trials (e.g., overlapping windows relative to stimulus onset or choice onset)
% standardize   zscore the matrix (each row / unit), yes or no (1 or 0)
% plotresults   yes if set to 'plot'
% 
% OUTPUT
% score         PC scores are the representations of X in PC space. Rows correspond to PCs, columns correspond to observations.
% loadings      loading of each unit for the first PCs
%
% NOTE OF CAUTION:
% My input matrix (rows are units, columns are observations) is transposed relative to the required input to pca.m
% Similarly, the score output is transposed such that rows are PCs (abstracted population firing patterns) and columns are observations.
%
% MORE INFORMATION ABOUT THE INPUT MATRIX
% For Sarah's extinction data, I used spike counts in 500-ms windows with the following columns:
% Columns 1:9     -0.5:0.5:2 s FS1, leave, choice, reward (correct trials only)
% Columns 10:18   as above, but for FS2
% Columns 19:27   as above, but for NS1 correct
% Columns 28:36   as above, but for NS1 error
% Columns 37:45   as above, but for NS2 correct
% Columns 46:54   as above, but for NS2 error
%% standardize, if requested
if standardize
  smat4dimred = zscore(smat4dimred,0,2);  % zscore each row (unit)
end
%% dimensionality reduction proper
% https://builtin.com/data-science/step-step-explanation-principal-component-analysis

% delete rows with NaN; rows are commonly units, columns are observations (e.g., trials or epochs within trials)
smat4dimred(sum(isnan(smat4dimred),2)>1,:) = [];

% transpose smat4dimred, because rows must be observations, columns must be variables (neurons)
[coeff,score,~,~,explained] = pca(smat4dimred','Economy',true); % setting 'Economy' to true reduces the number of returned PCs
score = score';

% get loadings of units for the first principal components
loadings = coeff;

% [coeff,score,latent,tsquared,explained,mu] = pca(X,'Economy','on');    % rows of X (n-by-p) be n observations, columns are p variables
%
%   coeff       Coeff contains PC coefficients (a.k.a. loadings; p-by-p matrix); each column contains coefficients for one PC, and PCs are sorted
%               in descending order of component variance.
%               Plotting the first column e.g. provides the first principal component; this is a weight vector for the p units.
% 	scores      Scores are representations of X in PC space. Rows of score are observations, columns correspond to PCs.
%               Plotting the first column (row in the output of this wrapper function) gives the trajectory of the first PC over observations.
%   latent      PC variances (eigenvalues of the covariance matrix of X)
%   explained   contains the percentage of the total variance explained by each principal component
%   mu          gives the estimated mean of each variable
%
% 'Economy','on': Indicator for the economy size output when the degrees of freedom, d, is smaller than the number of variables, p,
% specified as the comma-separated pair consisting of 'Economy' and a logical expression.
% pca returns only the first d elements of latent and the corresponding columns of coeff and score.
% This option can be significantly faster when the number of variables p is much larger than d.
%% plotting in 3D and 2D

if strcmp(plotresults,'plot')
  figure('position',[0.15,0.2,0.7,0.7])

  subplot(351),title('Raw'),hold on,imagesc(smat4dimred),xlabel('Observations'),ylabel('Units'),set(gca,'YDir','reverse'),axis tight

  subplot(352),title('Coefficients'),hold on,imagesc(coeff),xlabel('Observations'),ylabel('Units'),set(gca,'YDir','reverse'),axis tight

  subplot(353),title('Scores'),hold on,imagesc(score),xlabel('Observations'),ylabel('Principal components'),set(gca,'YDir','reverse'),axis tight

  subplot(354),title('Covariance matrix'),hold on,imagesc(cov(smat4dimred')),set(gca,'YDir','reverse'),axis tight

  subplot(355),title('Explained variance'),hold on,bar(explained),plot(cumsum(explained)),xlabel('Principal component'),ylabel('Percent explained variance'),axis([0,20,0,100]),grid on

  % 3D plot for PCs 1 thru 3
  subplot(356)
  nobs = size(smat4dimred,2); % number of observations
  plot3(score(1,:),score(2,:),score(3,:),'Color',[.5,.5,.5]),hold on
  colors = hot(nobs);
  for i = 1:nobs
    plot3(score(1,i),score(2,i),score(3,i),'.','Color',colors(i,:),'MarkerSize',15)
  end
  xlabel('PC 1'),ylabel('PC 2'),zlabel('PC 3')
  grid on

  % 2D plots for combinations of PCs 1 thru 4
  for i = 1:4
    subplot(3,5,i+6),hold on
    switch i
      case 1,dims = [1,2];xlabel('PC 1'),ylabel('PC 2')
      case 2,dims = [1,3];xlabel('PC 1'),ylabel('PC 3')
      case 3,dims = [2,3];xlabel('PC 2'),ylabel('PC 3')
      case 4,dims = [1,4];xlabel('PC 1'),ylabel('PC 4')
    end
    plot(score(dims(1),:),score(dims(2),:),'Color',[.5,.5,.5])
    for j = 1:nobs
      plot(score(dims(1),j),score(dims(2),j),'.','Color',colors(j,:),'MarkerSize',15);
    end
    grid on
    axis([floor(min(score(dims(1),:)))-1,ceil(max(score(dims(1),:)))+1,floor(min(score(dims(2),:)))-1,ceil(max(score(dims(2),:)))+1])
  end

  % 2D plots for PCs 1 thru 5 as a function of time windows
  for i = 1:5
    subplot(3,5,i+10),title(['PC ',num2str(i)]),hold on
    plot(score(i,:),'Color',[.5,.5,.5])
    for j = 1:nobs
      plot(j,score(i,j),'.','Color',colors(j,:),'MarkerSize',15);
    end
    ylabel(['PC ',num2str(i)])
    axis([0,nobs+1,floor(min(score(i,:)))-1,ceil(max(score(i,:)))+1])
  end
end
end

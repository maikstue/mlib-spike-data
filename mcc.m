function [L_ratio,ID,A,J3,PsF,DBVI] = mcc(X,cid,whichclusters,ts,whichPCs,nocorrs)
% function [L_ratio,ID,A,J3,PsF,DBVI] = mcc(X,cid,whichclusters,ts,whichPCs,nocorrs)
% compares two or more spike clusters regarding their separability in terms of J3, PseudoF, Davies-Bouldin validity index (DBVI)
%
% INPUT
% X               X is a matrix in which rows correspond to observations (spikes) and columns are voltage measurements
%                 mcc performs an input check and will transpose the matrix if size(X,1)<size(X,2)
% cid             column vector providing cluster ID for each observation
%                 numel(cid) should equal size(X,1)
% whichclusters   vector denoting which of the cluster IDs are to be included
% ts              vector of spike timestamps (optional)
% whichPCs        scalar/vector to specify which PCs should be included for calculation of L_ratio, ID, and A'
%                 if 'all', all PCs are included
% nocorrs         if 'nocorrs', no auto- and cross-correlograms are computed (optional)
%
% BACKGROUND INFORMATION ON J3, PSEUDO F, AND DAVIES-BOULDIN VALIDITY INDEX
% J1 reflects the average distance of points from their cluster mean and is minimized for compact clusters, whereas J2 reflects the average distance 
% between clusters and is maximized for well-separated cluster means. Hence, the ratio J3=J2/J1 is maximized for well-separated compact clusters.
% So, J3 depicts the ratio of between-cluster to within-cluster scatter; high values signify good isolation.
% DB reflects the ratio of the sum of within-cluster scatter to between-cluster separation, and low values show good isolation.
% J3(2D) should be >2, PseudoF(2D) should be >10,000, DBVI(2D) should range from 0.1-0.6
% Nicolelis et al., PNAS 2003
% Wheeler BC: Automatic discrimination of single units. In: Methods for Neural Ensemble Recordings, search in Google Books
% 
% BACKGROUND INFORMATION ON L_RATIO
% L_ratio is L/N, where N is the total number of spikes in a cluster, and L is the sum (across clusters) of cumulative chi-square distributions defined
% by the squared Mahalanobis distances of the individual spikes from their cluster. L_ratio is low (<1) for well-isolated clusters.
% "L_ratio is a measure of the amount of noise observed in the vicinity of the cluster, and isolation distance is a measure of how distant a 
% cluster is from the noise distribution." L_ratio correlates with the number of Type-1- and Type-2-errors.
% The authors do not provide guidelines to the interpretation of L_ratio, but the well-separated examples have values below 0.002, the intermediate
% examples is at 0.07, poor separation at 0.67.
% Schmitzer-Torbert and Redish, Neuroscience 2005
% 
% BACKGROUND INFORMATION ON ISOLATION DISTANCE
% Harris et al. (Neuron, 2001): "[Isolation distance is] defined to be the Mahalanobis distance from the identified cluster within which as many spikes 
% belong to the specified cluster as to others. This measure estimates the distance from this cluster to the nearest other cluster. We only considered 
% units with isolation distances >50".
% Schmitzer-Torbert and Redish, Neuroscience 2005: "sets. If a cluster contains nC cluster spikes, the Isolation Distance of the cluster is the D2 value of 
% the nCth closest noise spike. Isolation Distance is therefore the radius of the smallest ellipsoid from the cluster center containing all of the cluster 
% spikes and an equal number of noise spikes. As such, Isolation Distance estimates how distant the cluster spikes are from the other spikes recorded on 
% the same electrode. Isolation Distance is not defined for cases in which the number of cluster spikes is greater than the number of noise spikes.
% 
% HISTORY
% Nov 2024    added isolation distance and input argument 'whichPCs'
% Feb 2020    MCS reintroduced mcc code (J3, PseudoF, DB, linear discriminant analysis) to make the function more comprehensive
% Jan 2020    FJ build on mcc by MCS stripped of unnecessary code for L_ratio
% 
% Maik C. Stuettgen & Frank Jaekel, February 2016 and January 2020
%% preps
if ~isfloat(X),X = double(X);end      % make waveform matrix floating-point
if size(X,1)<size(X,2),X = X';end     % make sure that rows are spikes and columns are ticks
if numel(cid)~=size(X,1),error('input error'),end %#ok<ERTAG>   % length(cid) must match size(X,1), the number of spikes
if isempty(ts),disp('no timestamps provided, substituting...'),nots=1;ts=1:size(X,1);end % function can run with and without time stamps

if ischar(whichPCs) & strcmp(whichPCs,'all')
  whichPCs = 1:size(X,2);
elseif ~ischar(whichPCs) & numel(whichPCs)<3
  disp('Input argument whichPCs needs to have at least three entries, substituting with 1:3.')
  whichPCs = 1:3;
elseif ischar(whichPCs) & ~strcmp(whichPCs,'all')
  error('Input variable ''whichPCs'' is not specified correctly.')
end

ids = unique(cid);
if numel(ids)~=numel(whichclusters)
  ids = whichclusters(:);
  X(~ismember(cid,ids),:) = [];
  ts(~ismember(cid,ids))  = [];
  cid(~ismember(cid,ids)) = [];
end
clear whichclusters
figure('name','Cluster comparison','units','normalized','position',[.2,.2,.6,.6])

disp('--------------------------------')
disp('Note that values for all these measures depend on several factors such as how many samples per waveform are provided, ')
disp('which dimensions are used (e.g., PCs 1-3 or all PCs, waveform energy...), and whether a noise cluster is included or not.')
disp('Therefore, comparison of values across data sets / publications is not straightforward.')
disp('For the same reason, cutoff values for inclusion may be specific for data sets.')

%% subplot 331: waveforms
subplot(331),hold on
global colors
colors = parula(numel(ids)+1);
n = 100;  % how many waveforms to plot for each cluster
x = 1;
for i=1:numel(ids),plot(1,0,'Color',colors(i,:));end  % an admittedly inelegant way to make the legend appropriate
for i=ids'
  n = min([100,sum(cid==i)]);
  plot(X(randsample(find(cid==i),n),:)','Marker','none','LineStyle','-','Color',colors(x,:));
  x = x+1;
end
legend(num2str(ids),'Location','SouthEast')
axis tight off
clear i x

%% subplots 332,333: PCA variance explained and PC 1 over time
[~,score,~,~,explained] = pca(zscore(X));
% [coeff,score,latent,tsquared,explained,mu] 

subplot(332),title(['PCs 1-3 explain ',num2str(sum(explained(1:3)),'%2.1f'),'% of the variance']),hold on
bar(explained,'FaceColor','k'),hold on
plot(cumsum(explained),'r','LineWidth',2)
axis([0 min([size(X,2),10]),0,80])
xlabel('Principal components'),ylabel('Variance explained (%)')

subplot(333),hold on
x = 1;
for i=ids'
  plot(ts(cid==i)/60,score(cid==i,1),'Marker','.','MarkerSize',1,'LineStyle','none','Color',colors(x,:))
  x = x+1;
end
axis tight
xlabel('Time (min)')
ylabel('PC 1')

%% calculate L_ratio and Isolation Distance
f = score(:,whichPCs);
d = size(f,2);

% overall mean (not means of means? Nicolelis is not super clear on this)
M = mean(f,1);

n  = size(f,1);        % overall number of waveforms 
G  = size(ids,1);      % number of clusters
m  = nan(G,size(f,2)); % means for each cluster
C  = nan(d,d,G);       % covariances for each cluster
invC = nan(d,d,G);     % inverse covariance matrix
D2 = nan(G,n);         % squared Mahalanobis distance from cluster mean
D  = nan(G,n);         % Mahalanobis distance from cluster mean
N  = nan(G,1);         % number of waveforms per cluster  
for i = 1:numel(ids')
  select = (cid==ids(i));
  N(i) = sum(select);

  % compute the empirical mean
  m(i,:) = mean(f(select,:));

  % estimate covariance matrix
  C(:,:,i) = cov(f(select,:));

  % and the inverse for the Mahalanobis distance
  invC(:,:,i) = inv(C(:,:,i));

  % subtract mean from all feature vectors
  x = f-repmat(m(i,:),n,1);

  % compute distance from mean for all data points
  for j = 1:n
    D2(i,j) = x(j,:) * squeeze(invC(:,:,i)) * x(j,:)';
  end
  D(i,:) = sqrt(D2(i,:));
end

% given D2 it is easy to compute L(C)
L = nan(G,1);
L_ratio = nan(G,1);
disp('--------------------------------')
disp('L_ratio for single units should be <0.2, ideally <0.05')
for i=1:numel(ids')
  L(i)=sum(1-chi2cdf(D2(i,not(cid == ids(i))),d));
  L_ratio(i)=L(i)/N(i);
  disp(['L for cluster ',num2str(ids(i)),':       ',num2str(L(i),'%2.3f')])
  disp(['L_ratio for cluster ',num2str(ids(i)),': ',num2str(L_ratio(i),'%2.3f')])
end

% now compute isolation distance
disp('--------------------------------')
disp('Isolation Distance for single units should be >10, ideally >20')
for i=1:numel(ids')
  if N(i)>(j-N(i))
    disp(['Isolation Distance for cluster ',num2str(ids(i)),' is not defined.'])
    ID(i) = nan;
  else
    D2clust = sort(D2(i,cid==ids(i)));                            % D2clust contains only D2 values from this cluster, in ascending order
    D2noise = sort(D2(i,cid~=ids(i)));
    % D2noise = D2;D2noise(i,:) = [];D2noise(:,cid==ids(i)) = [];
    % D2noise = sort(D2noise(:));                                   % D2noise contains only D2 values from the other clusters, in ascending order
    try
      ID(i) = D2noise(find(D2noise>max(D2clust),1,'first'));
    catch
      ID(i) = ceil(max(D2clust));    % do this if no D2noise value is larger than max(D2clust)
    end
    disp(['Isolation Distance for cluster ',num2str(ids(i)),': ',num2str(ID(i),'%2.1f')])
  end
end

%% Fisher's linear discriminant analysis - commented out!
% This was used before Frank came up with the idea to use A' which uses the optimal non-linear decision bound. Linear discriminant analysis only finds
% the best linear decision bound.
pairs = nchoosek(ids,2);

% linecols = gray(size(pairs,1)+1);
% for i = 1:size(pairs,1)
%   obj = fitcdiscr(score(ismember(cid,pairs(i,:)),:),cid(ismember(cid,pairs(i,:))));
%   a = obj.Coeffs(1,2).Linear;
%   b = obj.Coeffs(1,2).Const;
%   
%   % plot the curve a*[PC1,PC2]+b=0
%   subplot(333)
%   f = @(PC1,PC2) a(1)*PC1 + a(2)*PC2 + b;
%   h = ezplot(f,[get(gca,'XLim') get(gca,'YLim')]);
%   h.Color = linecols(i,:);
%   xlabel('PC 1'),ylabel('PC 2')
%   
%   % plot the curve a*[PC3,PC4]+b=0
%   subplot(334)
%   f = @(PC3,PC4) a(3)*PC3 + a(4)*PC4 + b;
%   h = ezplot(f,[get(gca,'XLim') get(gca,'YLim')]);
%   h.Color = linecols(i,:);
%   xlabel('PC 3'),ylabel('PC 4')
%   
%   subplot(339),title('Pairs')
%   axis off
%   text(0.4,0.9-i*0.1,num2str(pairs(i,:)),'Color',linecols(i,:),'FontSize',12)
%   
%   % resubstitution error
%   resuberror = resubLoss(obj);
%   disp(' ')
%   disp('************************************************************')
%   disp(['***** CLUSTER PAIR: ' num2str(pairs(i,:)) ' ***********************************'])
%   disp('************************************************************')
%   disp(['Number of misclassified data points: ' num2str(resuberror*obj.NumObservations) ' out of ' num2str(obj.NumObservations)])
%   
%   % confusion matrix
%   cMat     = confusionmat(obj.Y,resubPredict(obj));
%   cMatFrac = cMat/obj.NumObservations;
%   disp('------------------------------------------------------------')
%   disp('Confusion matrix (absolute numbers):')
%   disp(['    |' sprintf('%6.0f',pairs(i,1)) '|' sprintf('%6.0f',pairs(i,2)) '|'])
%   disp('-------------------------')
%   disp([num2str(pairs(i,1)) '   |' sprintf('%6.0f',cMat(1,1)) '|' sprintf('%6.0f',cMat(1,2)) '|'])
%   disp('-------------------------')
%   disp([num2str(pairs(i,2)) '   |' sprintf('%6.0f',cMat(2,1)) '|' sprintf('%6.0f',cMat(2,2)) '|'])
%   disp(' ')
%   disp('------------------------------------------------------------')
%   disp('Confusion matrix (fractions): ')
%   disp(['    |' sprintf('%6.0f',pairs(i,1)) '|' sprintf('%6.0f',pairs(i,2)) '|'])
%   disp('-------------------------')
%   disp([num2str(pairs(i,1)) '   |' sprintf('%6.3f',cMatFrac(1,1)) '|' sprintf('%6.3f',cMatFrac(1,2)) '|'])
%   disp('-------------------------')
%   disp([num2str(pairs(i,2)) '   |' sprintf('%6.3f',cMatFrac(2,1)) '|' sprintf('%6.3f',cMatFrac(2,2)) '|'])
% end
%% subplots 334:336 - projection plots
subplot(334);
projection(score,cid,ids,m,C,1,2);
subplot(335);
projection(score,cid,ids,m,C,2,3);
subplot(336);
projection(score,cid,ids,m,C,3,1);

%% subplot 337: plot difference of distances
subplot(337),hold on

if G==2

  % Frank's original code to compute D1-D2 which works for two clusters only
  Diff = D(1,:)-D(2,:);
  bins = linspace(-5,5,20);
  db = mean(diff(bins))/2;
  for i=1:numel(ids')
    select = (cid == ids(i));
    % in Matlab change 'FaceAlpha' to 0.5 to make plot nicer
    [nw,bw] = hist(Diff(select),bins);
    bar(bw,nw,'FaceColor',colors(i,:),'BarWidth',0.8,'FaceAlpha',0.5)
    %   hist(Diff(select),bins,'FaceColor',colors(i,:),'BarWidth',0.5);
  end
  clear bw nw
  set(gca,'YScale','log')
  xlabel('D_1-D_2')
  ylabel('Log frequency')

elseif G>2  % this does the computation for plotting, but we do not want to plot, or do we?

  % when we have >2 clusters, we need to compute pairwise distances and pairwise A' values; another plot is warranted
  bins = linspace(-5,5,20);
  db = mean(diff(bins))/2;
  for i = 1:G
    Dnoise = D; Dnoise(i,:) = [];
    Diff = D(i,:) - Dnoise;
    select = (cid == ids(i));
    for j = 1:G-1 % for each pair
      [nw,bw] = hist(Diff(j,select),bins);
    end
  end
end

%% subplot 338: plot ROC curve-like plot and compute A'
if G==2

  % Frank's original code to compute D1-D2 which works for two clusters only
  [Diff,order] = sort(Diff);
  cs = nan(G,n);
  for i=1:numel(ids')
    cs(i,:) = cumsum(cid(order)==ids(i))/N(i);
  end
  subplot(338)
  plot(cs(1,:),1-cs(2,:));
  axis square
  xlabel('Class 1 hitrate')
  ylabel('Class 2 hitrate')
  A = sum(diff([cs(1,:),1]).*(1-cs(2,:)));
  title(['A''=' num2str(round(A*100)/100)])

elseif G>2
  for i = 1:G
    Dnoise = D; % Dnoise(i,:) = [];
    Diff = D(i,:) - Dnoise;
    for j = 1:G % for each pair
      [Diffdummy,order] = sort(Diff(j,:));
      cs = nan(G,n);
      for k=1:numel(ids')
        cs(k,:) = cumsum(cid(order)==ids(k))/N(k);
      end
      A(i,j,1) = sum(diff([cs(1,:),1]).*(1-cs(2,:)));
      clear Diffdummy order
    end
  end
  A(A<0.5) = 1-A(A<0.5);
  imagesc(A,[0.5,1])
  A(logical(eye(size(A)))) = 0.5;
  colormap(hot),colorbar
  axis tight
  set(gca,'XTick',1:G,'YTick',1:G,'Ydir','reverse','XTickLabel',{num2str(ids)},'YTickLabel',{num2str(ids)})
  xlabel('Clusters'),ylabel('Clusters')
  title('Pairwise A''')

  subplot(338),title('Pairwise A'''),hold on
  text(0,.5,num2str(A,'%10.2f'))
  axis off
end

%% subplot 339: text
[J3,PsF,DBVI] = getJ3(score,cid,ids);
subplot(339),hold on
axis off
fsize = 8;
text(0,1.0,['J3   = ',num2str(J3,'%1.3f')],'FontSize',fsize)
text(0,0.9,['PsF  = ',num2str(PsF,'%1.3f')],'FontSize',fsize)
text(0,0.8,['DBVI = ',num2str(DBVI,'%1.3f')],'FontSize',fsize)
for i = 1:numel(ids)
  text(0,0.7-(i*0.1),['L_ratio cluster ',num2str(ids(i)),' with ',num2str(sum(cid==ids(i))),' spikes: ',num2str(L_ratio(i),'%2.3f')],'FontSize',fsize,'Interpreter','none')
end
for i = 1:numel(ids)
  text(0,0.7-(i*0.1)-numel(ids)*0.1-0.1,['Isolation Distance cluster ',num2str(ids(i)),': ',num2str(ID(i),'%2.3f')],'FontSize',fsize,'Interpreter','none')
end

%% extra figure: cross-correlation histograms
% not really cross-correlograms, as the ordinate is scaled in units of events or firing rate, not correlation
if ~exist('nocorrs','var')
  if numel(ids)<6
    nPlots = numel(ids)^2;
    if nPlots==4
      figure('name','cross-correlograms','units','normalized','position',[0.6 0.6 0.2 0.3])
    elseif nPlots==9
      figure('name','cross-correlograms','units','normalized','position',[0.5 0.4 0.4 0.5])
    else
      figure('name','cross-correlograms','units','normalized','position',[0.4 0.3 0.5 0.6])
    end
    
    % plot auto-correlograms
    subplotIDs = cumsum([1,ones(1,numel(ids)-1)+numel(ids)]);
    for i = 1:numel(ids)
      subplot(numel(ids),numel(ids),subplotIDs(i))
      title(['AC Cluster ' num2str(ids(i))]),hold on
      psth = mpsth(ts(cid==ids(i)),ts(cid==ids(i)),'pre',50,'post',50,'fr',0);
      psth(psth(:,1)==0,2) = 0;
      bar(psth(:,1),psth(:,2),'k')
      axis tight
      xlabel('Time (s)'),ylabel('Counts per bin')
    end
    
    % plot cross-correlograms
    switch numel(ids)
      case 2,subplotIDs = 3;
      case 3,subplotIDs = [4,7,8];
      case 4,subplotIDs = [5,9,13,10,14:15];
      case 5,subplotIDs = [6,11,16,21,12,17,22,18,23,24];
    end
    if numel(ids)>1
      for i = 1:size(pairs,1)
        subplot(numel(ids),numel(ids),subplotIDs(i))
        title(['CC Clusters ' num2str(pairs(i,:))]),hold on
        psth = mpsth(ts(cid==pairs(i,1)),ts(cid==pairs(i,2)),'pre',50,'post',50,'fr',0);
        bar(psth(:,1),psth(:,2))
        axis tight
        xlabel('Time (s)'),ylabel('Counts per bin')
      end
    end
  else
    error('Sorry, correlograms for >5 clusters is not supported at the moment.')
  end
end

end
% NESTED FUNCTIONS
%% nested function projection(score,cid,ids,m,C,PC1,PC2)
function projection(score,cid,ids,m,C,PC1,PC2)
alpha = 0.8;
hold on
x=1;
global colors
for i=1:numel(ids')   % for all clusters
  plot(score(cid==ids(i),PC1),score(cid==ids(i),PC2),'Marker','.','MarkerSize',1,'LineStyle','none','Color',colors(x,:))
  plot(m(i,PC1),m(i,PC2),'Linewidth',3,'MarkerSize',1,'Marker','x','LineStyle','none','Color',colors(x,:)*alpha);         % plot mean of cluster
  c = [sin(linspace(0,2*pi)); cos(linspace(0,2*pi))];
  E = squeeze(C(:,:,i));
  E = E([PC1 PC2],[PC1 PC2]);
  c = chol(E)'*c;
  plot(m(i,PC1)+c(1,:),m(i,PC2)+c(2,:),'Linewidth',3,'LineStyle','-','Color',colors(x,:)*alpha);    % plot cluster sphere
  x=x+1;
end
axis tight
xlabel(['PC ' num2str(PC1)]),ylabel(['PC ' num2str(PC2)])
end

%% nested function getJ3
function [J3,PsF,DBVI] = getJ3(score,cid,ids)
% function [J3,PsF,DBVI] = getJ3(score,cid,ids)
% compute J3, PsF, and DB measures of cluster separability
% J1 reflects the average distance of points from their cluster mean and is minimized for compact clusters, whereas
% J2 reflects the average distance between clusters and is maximized for well-separated cluster means.
% Hence, the ratio J3=J2/J1 is maximized for well-separated compact clusters.
% Nicolelis et al., PNAS 2003
% Wheeler BC: Automatic discrimination of single units. In: Methods for Neural Ensemble Recordings, search in Google Books
%
% J3 depicts the ratio of between-cluster to within-cluster scatter; high values signify good isolation.
% DB reflects the ratio  of  the  sum  of  within-cluster  scatter  to  between-cluster separation, and low values show good isolation.

inputMatrix = score(:,1:2);   % restrict computation of J3 to first two principal components (the more components, the lower J3; see Nicolelis et al.)
% compute grand centroid (is essentially zero because X has been z-scored)
grandcentroid = mean(inputMatrix);

% compute cluster centroids and J2
centroids = nan(numel(ids),size(inputMatrix,2));
J2        = nan(1,numel(ids));
for i = 1:numel(ids)
  centroids(i,:) = mean(inputMatrix(cid==ids(i),:));
  J2(i) = sum(cid==ids(i))*(pdist([centroids(i,:);grandcentroid]).^2);  % distance cluster centroid to grand centroid squared, weighted by number of waveforms
end
J2 = sum(J2);

% compute J1 for each cluster
meanDist = nan(1,numel(ids));   % needed later for DBVI
J1       = nan(1,numel(ids));
for i = 1:numel(ids)
  redInputMatrix = inputMatrix(cid==ids(i),:);
  dist = nan(size(redInputMatrix,1),1);
  for j = 1:size(redInputMatrix,1)
    dist(j) = pdist([redInputMatrix(j,:);centroids(i,:)])^2;
  end
  J1(i) = sum(dist);
  meanDist(i) = mean(dist);
  clear dist j redInputMatrix
end
J1 =  sum(J1);

% compute J3 and pseudoF
J3 =  J2/J1;
PsF = ((size(inputMatrix,1)-numel(ids))/(numel(ids)-1))*J3; % adjusting J3 to the number of clusters and valid spikes in the channel

% compute the Davies-Bouldin (DB) validity  index
% DBVI = (1/numel(ids)) * sum(max(Si+Sj/dij))
% Thus, for each cluster, find the complementary cluster for which the ratio (Si+Sj)/dij is maximal.
% Si is the average Euclidean distance of all waveforms from their cluster center.
% dij is the Euclidean distance of the cluster centers.

DBratio = nan(numel(ids));
DBmax   = nan(numel(ids),1);
for i = 1:numel(ids)
  for j = 1:numel(ids)
    DBratio(i,j) = (meanDist(i)+meanDist(j)) / pdist([centroids(i,:);centroids(j,:)]'); % added transpose on June 9, 2017 - did not change much
  end
  DBratio(i,isinf(DBratio(i,:))) = -1;
  DBmax(i) = max(DBratio(i,:));
end
DBVI = mean(DBmax);

% Computation of all measures are based on the first two principal components.
% Desired ranges of measures for single-unit channels taken from Nicolelis et al. (2003), PNAS.
disp('--------------------------------')
disp(['J3 (2D)   = ',num2str(J3,'%1.3f'),' - should be >2, mean of 6'])
disp(['PsF (2D)  = ',num2str(PsF,'%1.3f'),' - should be >10,000'])
disp(['DBVI (2D) = ',num2str(DBVI,'%1.3f'),' - should range from 0.1-0.6, mean 0.3'])
end

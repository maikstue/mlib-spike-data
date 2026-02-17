function [cid,avsil] = mcluster(data,z,nclus,whichPCs,cType,whichDist,metric,T,plotting)
% function [cid,avsil] = mcluster(data,z,nclus,whichPCs,cType,whichDist,metric,T,plotting)
% 
% performs clustering on a data matrix after running PCA and returns cluster codes cid and average silhouette values
% 
% mcluster will rearrange clusters such that cluster 1 has the most elements, cluster 2 the second most etc.
%
% data      input matrix where rows are units and columns are binned firing rates or eta2/omega2 values
% z         'zscore' or anything else - zscore data row-wise or not
% nclus     desired number of clusters
% whichPCs  number of principal components to use for clustering; if 0, no pca is performed
%           by experience, the number does not really matter as long as it's at least the first three
% cType     clustering method; may be 'hieragglom' (for hierarchical agglomerative clustering) or 'kmeans'
% whichDist specifies how to measure distance between clusters - 'average', 'centroid', 'median' make sense (type doc linkage)
%           others are 'complete','single','ward','weighted'
% metric    specifies which distance measure to use to compute distances between the rows of data - 'Euclidean', 'correlation' make sense
%           others are 'squaredeuclidean', 'mahalanobis', 'cityblock'...
% T         threshold for colors in dendrogram; if [], not dendrogram will be constructed or plotted
% plotting  if 'yes', data will be plotted
%
% OUTPUT
% cid       cluster ID for each row of the data matrix
% avsil     average silhouette value for each cluster 1:nclus; see below for a brief explanation
% 
% A note on silhouette values
% A silhouette value measures how similar a point is to points in its own cluster when compared to points in other clusters.
% A high silhouette value indicates that a point is well matched to its own cluster, and poorly matched to other clusters.
% The silhouette value si for the ith point is defined as si = (bi - ai) / max(ai,bi), where ai is the average distance from the ith point 
% to the other points in the same cluster as i, and bi is the minimum average distance from the ith point to points in a different cluster, 
% minimized over the clusters. If the ith point is the only point in its cluster, then the silhouette value si is set to 1.
% towardsdatascience.com says: "The optimal number of clusters k is the one that maximizes the average silhouette over a range of possible values for k."
% 
% Maik C. Stüttgen, August 2017
%
% HISTORY
% 2024-03-05    added average silhouette values as output argument
% 2022-05-07    clusters are now sorted according to number of elements
% 2022-04-30    added extra input checks for cType, whichDist, and metric
% 2020-02-17    added input argument 'plotting'; entering [] for T skips potting the dendrogram even if plotting is 'yes'
%% input check
if size(data,1)<size(data,2),disp('Sure the data matrix is o.k.? There are more columns than rows.'),end
if ~any(strcmp(cType,{'hieragglom','kmeans'})),error('Invalid clustering method specified. Use ''kmeans'' or ''hieragglom''.'),end
if ~any(strcmp(whichDist,{'average','centroid','median','complete','single','ward','weighted'})),error('Invalid distance specified. Use e.g. ''average'' or ''centroid''.'),end
if ~any(strcmp(metric,{'Euclidean','correlation','squaredeuclidean','mahalanobis','cityblock','seuclidean','hamming','spearman'})),error('Invalid metric specified. Use e.g. ''Euclidean'' or ''correlation''.'),end
disp('------------------------------------------------------------')
if any(isnan(data(:)))
  disp(['Found ',num2str(sum(isnan(data(:,1)))),' NaNs in column 1.'])
  try
    disp(['Found ',num2str(sum(isnan(data(:,2)))),' NaNs in column 2.'])
    disp(['Found ',num2str(sum(isnan(data(:,3)))),' NaNs in column 3.'])
  end
  disp('Rows with NaNs will be ignored.')
end
%% zscore data if desired
if strcmp(z,'zscore'),disp('Data matrix is z-scored row-wise.'),data = zscore(data,0,2);end
%% perform PCA on data
if whichPCs>size(data,2),whichPCs=size(data,2);end  % we can only have so many PCs as there are columns in the data matrix
[coeff,score,~,~,explained] = pca(data);
%% linkage and clustering
if strcmp(cType,'hieragglom')
  if whichPCs>0
    disp('Hierarchical agglomerative clustering on the basis of principal components.')
    tree = linkage(score,whichDist,metric);
  else
    disp('Hierarchical agglomerative clustering on the basis of raw data.')
    tree = linkage(data,whichDist,metric);
  end
  cid = cluster(tree,'maxclust',nclus);
elseif strcmp(cType,'kmeans')
  if whichPCs>0
    disp('K-means clustering on the basis of principal components.')
    cid = kmeans(score(:,1:whichPCs),nclus);%,'Display','iter');
  else
    disp('K-means clustering on the basis of raw data.')
    cid = kmeans(data,nclus);%,'Display','iter');
  end
else
  error(['Unknown cluster method requested: ',cType])
end
%% rearrange clusters such that numel(cid==1)>numel(cid==2) etc.
nelements = histcounts(cid,0.5:1:nclus+0.5);ranks = tiedrank(nelements);
if numel(unique(ranks))<nclus % happens if there are ties
  ranks = ranks + 0.1*(rand(1,numel(ranks))-0.5); % we add small random numbers to all values and then redetermine the rank
  ranks = tiedrank(ranks);
end
ranks = tiedrank(-ranks);
c2 = nan(numel(cid),1);
for i = 1:nclus
  c2(cid==i) = ranks(i);
end
cid = c2;clear c2
for i = 1:nclus
  disp(['Cluster ',num2str(i),': ',num2str(sum(cid==i))])
end
s = silhouette(data,cid,metric);    % default is Euclidean squared
for i = 1:nclus
  avsil(i) = mean(s(cid==i));
end

%% plotting cluster results
if strcmp(plotting,'yes')

  cols  = lines(max([whichPCs,nclus])); % for clusters; use parula, lines, jet

  if whichPCs<1 & ~strcmp(z,'zscore'),figname = 'Clustering on raw data';
  elseif whichPCs<1 & strcmp(z,'zscore'),figname = 'Clustering on z-scored data';
  elseif whichPCs>0 & ~strcmp(z,'zscore'),figname = ['Clustering on raw data with ',num2str(max(whichPCs)),' PCs'];
  elseif whichPCs>0 & strcmp(z,'zscore'),figname = ['Clustering on z-scored data with ',num2str(max(whichPCs)),' PCs'];
  end

  figure('name',figname,'position',[.05 .1 .3 .8])
  clear figname

  subplot(321),title('Explained variance'),hold on
  bar(explained,'FaceColor',[0.5 0.5 0.5]),plot([0.5 50.5],[5 5],'k:'),plot([0.5 50.5],[50 50],'k:'),plot([0.5 50.5],[75 75],'k:')
  plot(cumsum(explained),'k','Marker','o')
  axis([0.5 12.5 0 100])
  xlabel('First 10 PCs'),ylabel('PEV (%)')
  
  subplot(322),title('Scores'),hold on
  for i = 1:nclus
    if numel(cid)<1000  % adjust marker size to number of elements - better to visualize densities
      scatter(score(cid==i,1),score(cid==i,2),4,'MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:)),hold on  % scatter(x,y,sz)
    else
      scatter(score(cid==i,1),score(cid==i,2),1,'MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:)),hold on  % scatter(x,y,sz)
    end
  end
  axis([floor(min(score(:,1))) ceil(max(score(:,1))) floor(min(score(:,2))) ceil(max(score(:,2)))])
  xlabel('PC 1'),ylabel('PC 2')
  l = legend(num2str((1:nclus)'));
  set(l,'FontSize',7,'Location','best')
  
  subplot(323),title('Nanmeans for clusters'),hold on
  for i = 1:nclus
    plot(nanmean(data(cid==i,:)),'Color',cols(i,:),'LineWidth',2)
  end
  xlim([0 size(data,2)+1])
  
  subplot(324),title('First PCs'),hold on
  for i = 1:min([max(whichPCs),4])    % plot no more than the first 4 PCs
    plot(coeff(:,i),'Color',cols(i,:),'LineWidth',2)
  end
  xlim([0 size(data,2)+1])
  
  uavsil = mean(avsil);
  for i = 1:nclus,n(i) = sum(cid==i);end
  wavsil = mean((avsil.*n)/sum(n));
  subplot(325),title(['Av.silhouette ',num2str(wavsil,'%2.2f'),' (w), ',num2str(uavsil,'%2.2f'),' (u)']),hold on
  silhouette(data,cid,metric);

  subplot(326),title('All data, clustered'),hold on
  for i = 1:nclus
    plot(1:size(data,2),data(cid==i,:),'Color',cols(i,:),'LineWidth',0.25)
  end
  xlim([0,size(data,2)+1])
end
%% plotting the dendrogram
if strcmp(cType,'hieragglom') && strcmp(plotting,'yes') && ~isempty(T)
  figure('position',[.35 .1 .3 .25])
  dendrogram(tree,size(data,1),'ColorThreshold',T*max(tree(:,3)));
  % If ColorThreshold is 'default', then threshold is 70% of the maximum linkage, 0.7*max(tree(:,3)).
  ylim([0 ceil(max(tree(:,3)))])
else
  disp('No dendrogram constructed.')
end
end

function [obj,shiftData] = MeanShift2Dv(obj)

localData=obj.Data;
dataCount=length(obj.Data);
shiftData=zeros(size(localData));
for localDataIte=1:dataCount
    moved=true;
    while moved
        old=repmat(localData(localDataIte,:),dataCount,1);
        dista=Batch.DistanceDisjointVec(old,obj.Data);
        weight=zeros(dataCount,length(-obj.Wind:obj.Wind));
  
        for k=-obj.Wind:obj.Wind
           weight(:,1+obj.Wind+k)=mvnpdf(dista,[0 0]+[2*pi*k 0],obj.Bandwidth);
        end
        weight=sum(weight,2);
        shiftAngular=Batch.WeightedMeanCS(obj.Data,weight);
        shiftAngular= shiftAngular(1);
        shiftLinear=sum(obj.Data(:,2).*weight);
        scaleFactor=sum(weight);
        shiftAngular=wrapTo2Pi(shiftAngular);
        shiftLinear=shiftLinear/scaleFactor;
        
        old=localData(localDataIte,:);
        new=[shiftAngular,shiftLinear];
        dista=Batch.DistanceDisjoint(old,new);
        localData(localDataIte,:)=new;
        if (dista(1) < obj.Bandwidth(1,1) * obj.StopFraction)&&(dista(2) < obj.Bandwidth(2,2) * obj.StopFraction)
            moved=false;
            shiftData(localDataIte,:)=new;
        end % dista < obj.Bandwidth * obj.StopFraction
    end % moved
end % localDataIte=1:length(localData)
clusterCenters=[];
clusterID=-ones(dataCount,1);
for localDataIte=1:dataCount
    found=false;
    [clustCentLength,~]=size(clusterCenters);
    for clustCentIter=1:clustCentLength
        dista=obj.DistanceDisjoint(clusterCenters(clustCentIter,:),shiftData(localDataIte,:));
        if (dista(1) < obj.Bandwidth(1,1) * obj.StopFraction*10)&&(dista(2) < obj.Bandwidth(2,2) * obj.StopFraction*10)
            found = true;
            clusterID(localDataIte)=clustCentIter;
            z=exp(1i*shiftData(clusterID==clustCentIter,1));
            clusterCenters(clustCentIter,:)=Batch.Centroid(shiftData(clusterID==clustCentIter,:));
            break
        end % dista < obj.Bandwidth * obj.StopFraction
        
    end % clustCentIter=1:length(clustCent)
    if ~found
        if isempty(clustCentIter)
            clustCentIter=0;
        end % isempty(clustCentIter)
        clusterCenters(clustCentIter+1,:)=shiftData(localDataIte,:);
        clusterID(localDataIte)=clustCentIter+1;
    end % ~found
end % localDataIte=1:dataCount
 
 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Prune clusters with too few data points
edges=1:max(clusterID);
counts=histc(clusterID,edges);
if sum(counts)>0
    removeClusters=edges(counts<=1); % we are removing clustsers with 1 sample
    keepClusters=edges(counts>1); % we are keeping clusters with more then one sample
    clusterCenters=clusterCenters(keepClusters,:);
    for removeClustersIter=1:length(removeClusters)
        clusterID(clusterID==removeClusters(removeClustersIter))=NaN;
    end % removeClustersIter=1:length(removeClusters)
    existingLabels=unique(clusterID);
    existingLabels=existingLabels(~isnan(existingLabels));
    % organize nicely the labels for clusters
    for clustterIDIter=1:length(clusterID)
        if ~isnan(clusterID(clustterIDIter))
            clusterID(clustterIDIter)=find(existingLabels==clusterID(clustterIDIter));
        end % ~isnan(clusterID(clustterIDIter))
    end % clustterIDIter=1:length(clusterID)
end % sum(counts)>0
%--------------------------------------------------------------------------
% compute covariances for each cluster
ClusterCenters_MLE=[];
ClusterCovariance=[];
if sum(isnan(clusterID))==length(clusterID)
    return
end
for IDIter=1:max(clusterID)
    [ClusterCenters_MLE(IDIter,:),ClusterCovariance(:,:,IDIter)]=obj.mle(obj.Data(clusterID==IDIter,:));
end % IDIter=1:max(clusterID)
%--------------------------------------------------------------------------
% remove clusters with ildefined covariance matrices (points are close to
% eachother)

keepClusters=[];
removeClusters=[];
for IDIter=1:max(clusterID)
    [~,p]=chol(ClusterCovariance(:,:,IDIter));
    if p==0 &&  ClusterCovariance(1,1,IDIter)>10^(-10) && ClusterCovariance(2,2,IDIter)>10^(-10) && ~isnan(ClusterCovariance(1,1,IDIter)) && ~isnan(ClusterCovariance(2,2,IDIter))
        keepClusters(end+1)=IDIter;
    else
        removeClusters(end+1)=IDIter;
    end
end
clusterCenters=clusterCenters(keepClusters,:);
ClusterCovariance=ClusterCovariance(:,:,keepClusters);
ClusterCenters_MLE=ClusterCenters_MLE(keepClusters,:);
for removeClustersIter=1:length(removeClusters)
    clusterID(clusterID==removeClusters(removeClustersIter))=NaN;
end % removeClustersIter=1:length(removeClusters)
%--------------------------------------------------------------------------
% organize nicely the labels for clusters
existingLabels=unique(clusterID);
existingLabels=existingLabels(~isnan(existingLabels));
for clustterIDIter=1:length(clusterID)
    if ~isnan(clusterID(clustterIDIter))
        clusterID(clustterIDIter)=find(existingLabels==clusterID(clustterIDIter));
    end % ~isnan(clusterID(clustterIDIter))
end % clustterIDIter=1:length(clusterID)

obj.ClustersIDs=clusterID;
obj.ClustersMeans=clusterCenters;
obj.ClustersMeans_MLE=ClusterCenters_MLE;
obj.ClusterCovariance=ClusterCovariance;


end


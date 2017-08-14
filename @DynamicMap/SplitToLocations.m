function obj = SplitToLocations(obj)
if isempty(obj.LocationID)
    [idx,~]=rangesearch(obj.Position,obj.Grid,obj.Radious);
    for k=1:length(idx)
        if length(idx{k})>2 
            obj.Batches(k)=obj.Batches(k).AddData(obj.ThetaRho(idx{k},:));
            obj.Batches(k)=obj.Batches(k).SetParametersMeanShift(MeanShiftKernel.Gaussian2D,DistanceType.Wrapped,0.001,-1); % just an experiemntal values
        end
    end
else
    obj.BatchesSparse=Batch.empty;
    for j=1:max(obj.LocationID)
        loc=unique(obj.Position(obj.LocationID==j,:),'rows');
        loc=loc(1,:);
        obj.BatchesSparse(j)=Batch();
        obj.BatchesSparse(j)=obj.BatchesSparse(j).SetParameters(j,loc,obj.ThetaRho(obj.LocationID==j,:),obj.Wind);
        obj.BatchesSparse(j)=obj.BatchesSparse(j).SetParametersMeanShift(0.001,-1);
    end
end
end
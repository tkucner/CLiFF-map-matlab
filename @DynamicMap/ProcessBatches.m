function obj=ProcessBatches(obj)
localBatches=obj.Batches;

obj1 = ProgressBar(total_observations, 'Title', 'Processing Batches ...', 'Total', numel(localBatches));

for i=1:numel(localBatches)
    if ~isempty(localBatches(i).Data)
        localBatches(i)=localBatches(i).MeanShift2Dv();
        if ~isempty(localBatches(i).ClustersMeans)
            localBatches(i)=localBatches(i).EMv();
        end
    end
    obj1.step([], [], []);
end
obj.Batches=localBatches;
end
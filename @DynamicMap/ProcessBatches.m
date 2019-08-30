function obj=ProcessBatches(obj)
localBatches=obj.Batches;

obj1 = ProgressBar(numel(localBatches), 'Title', 'Processing Batches ...');

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
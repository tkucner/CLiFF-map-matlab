function obj=ProcessBatches(obj)
localBatches=obj.Batches;
fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,numel(localBatches)) '\n\n']);
for i=1:numel(localBatches)
    if ~isempty(localBatches(i).Data)
        localBatches(i)=localBatches(i).MeanShift2Dv();
        if ~isempty(localBatches(i).ClustersMeans)
            localBatches(i)=localBatches(i).EMv();
        end
    end
    fprintf('\b|\n');
end
obj.Batches=localBatches;
end
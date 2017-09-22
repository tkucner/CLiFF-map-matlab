function [ ] = SaveCSV(obj,filename)

AllMeans=vertcat(obj.Batches(:).Mean);
MaxSpeed=max(AllMeans(:,2));
RES=[];
for i=1:numel(obj.Batches)
   if ~isempty(obj.Batches(i).Mean)
       [U,V]=pol2cart(obj.Batches(i).Mean(:,1),obj.Batches(i).Mean(:,2));
       
       r=length(U);
       C=repmat(obj.Batches(i).Pose,r,1);
       RES=[RES;C(:,1) C(:,2) obj.Batches(i).Mean(:,1) obj.Batches(i).Mean(:,2)];
       
   end
end

csvwrite(filename,RES)
end

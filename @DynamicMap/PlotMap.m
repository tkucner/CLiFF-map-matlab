function [ ] = PlotMap(obj,filterTreshold)
figure('Name',strcat('Dynamic map for-> ', obj.File, ' with improtance treshold = ',num2str(filterTreshold)));
grid on;
hold on;
%---------------------
% stuff for nice backgroimgund image
if ~isempty(obj.Image)
    obj.Image = flipud(obj.Image);
    imagesc([obj.GridParameters(1) obj.GridParameters(2)], [obj.GridParameters(3) obj.GridParameters(4)], obj.Image);
end
% %---------------------
% plot the nodes of the map
% plot(obj.Grid(:,1),obj.Grid(:,2),'+g','MarkerSize',3);
AllMeans=vertcat(obj.Batches(:).Mean);
MaxSpeed=max(AllMeans(:,2));
NormalisationFactor=1/MaxSpeed;
for i=1:numel(obj.Batches)
   if ~isempty(obj.Batches(i).Mean)
       [U,V]=pol2cart(obj.Batches(i).Mean(:,1),NormalisationFactor*obj.Batches(i).Mean(:,2));
       filt=obj.Batches(i).P>filterTreshold; % filtering out modes with low likelihood
       %U=U*0.2;
       %V=V*0.2;
       U=U(filt);
       V=V(filt);
       r=length(U);
       C=repmat(obj.Batches(i).Pose,r,1);
       quiver(C(:,1),C(:,2),U.*obj.GridParameters(5),V.*obj.GridParameters(5),'AutoScale','off','color','red','MaxHeadSize',5)
       %text(obj.Grid(:,1),obj.Grid(:,2),num2str(i))
   end
end

%---------------------
% stuff for nice background image
if ~isempty(obj.Image)
    set(gca,'ydir','normal');
end
%---------------------
axis([obj.GridParameters(1)-obj.GridParameters(5) obj.GridParameters(2)+obj.GridParameters(5), obj.GridParameters(3)-obj.GridParameters(5) obj.GridParameters(4)+obj.GridParameters(5)])
end


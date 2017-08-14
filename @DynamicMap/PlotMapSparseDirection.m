function [ ] = PlotMapSparseDirection(obj,filterTreshold,scale)
switch nargin
    case 2
        scale=1;
    case 1
        filterTreshold=0;
        scale=1;
end
figure('Name',strcat('Dynamic map for-> ', obj.File, ' with improtance treshold = ',num2str(filterTreshold)));
grid on;
hold on;
%---------------------

% %---------------------
% plot the nodes of the map
plot(obj.Grid(:,1),obj.Grid(:,2),'+','MarkerSize',3);

col=hsv(360);
for i=1:numel(obj.BatchesSparse)
   if ~isempty(obj.BatchesSparse(i).Mean)
       
       [U,V]=pol2cart(obj.BatchesSparse(i).Mean(:,1),obj.BatchesSparse(i).Mean(:,2));
       filt=obj.BatchesSparse(i).P>filterTreshold; % filtering out modes with low likelihood
       U=U*3;
       V=V*3;
       U=U(filt);
       V=V(filt);
       [th,~]=cart2pol(U,V);
       r=length(U);
       C=repmat(obj.BatchesSparse(i).Pose,r,1);
       for j=1:r
            c_id=mod(round(rad2deg(wrapTo2Pi(th(j)))),360);
            quiver(C(j,1),C(j,2),U(j).*obj.GridParameters(5)*scale,V(j).*obj.GridParameters(5)*scale,'AutoScale','off','color',col(c_id+1,:),'MaxHeadSize',10,'LineWidth',1)
       end
       %text(obj.Grid(:,1),obj.Grid(:,2),num2str(i))
   end
end


%---------------------
axis([obj.GridParameters(1) obj.GridParameters(2), obj.GridParameters(3) obj.GridParameters(4)])
colormap(col)
c=colorbar('Ticks',[0:0.25:1],'TickLabels',[0 90 180 270 360]);
c.Label.String = 'Orientation [deg]';
end


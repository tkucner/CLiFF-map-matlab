function [ ] = PlotMapDirection(obj,filterTreshold,scale)
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
% plot(obj.Grid(:,1),obj.Grid(:,2),'+','MarkerSize',3);
if ~isempty(obj.Image)
    obj.Image = flipud(obj.Image);
    obj.Image = (double(obj.Image)/255.0 > 0.65) - 1;
    imagesc([obj.GridParameters(1) obj.GridParameters(2)], [obj.GridParameters(3) obj.GridParameters(4)], obj.Image);
end

col=hsv(360);
for i=1:numel(obj.Batches)
   if ~isempty(obj.Batches(i).Mean)
       
       [U,V]=pol2cart(obj.Batches(i).Mean(:,1),obj.Batches(i).Mean(:,2));
       filt=obj.Batches(i).P>filterTreshold; % filtering out modes with low likelihood
       U=U*3;
       V=V*3;
       U=U(filt);
       V=V(filt);
       [th,~]=cart2pol(U,V);
       r=length(U);
       C=repmat(obj.Batches(i).Pose,r,1);
       for j=1:r
            c_id=mod(round(rad2deg(wrapTo2Pi(th(j)))),360);
            quiver(C(j,1),C(j,2),U(j).*obj.GridParameters(5)*scale,V(j).*obj.GridParameters(5)*scale,'AutoScale','off','color',col(c_id+1,:),'MaxHeadSize',10,'LineWidth',1)
       end
       %text(obj.Grid(:,1),obj.Grid(:,2),num2str(i))
   end
end

%---------------------
% stuff for nice background image
if ~isempty(obj.Image)
    set(gca,'YDir','normal');
end
%---------------------
axis([obj.GridParameters(1) obj.GridParameters(2), obj.GridParameters(3) obj.GridParameters(4)])
%colormap(col)
c=colorbar('Ticks',[0:0.25:1],'TickLabels',[0 90 180 270 360]);
c.Label.String = 'Orientation [deg]';
axis image
l = legend("arrows");
set(l, 'visible', 'off')
end


function [] = PlotUVDirection(obj,scale)
switch nargin
    case 1
        scale=1;
end
figure('Name',strcat('Colored Raw data for file->', obj.File));
grid on;
hold on;

[th,~]=cart2pol(obj.UV(:,1),obj.UV(:,2));
col=hsv(360);
for i=1:5:length(th)
    c_id=mod(round(rad2deg(wrapTo2Pi(th(i)))),360);
    %quiver(obj.Position(i,1),obj.Position(i,2),obj.UV(i,1),obj.UV(i,2),'AutoScaleFactor',0.09,'MaxHeadSize',100,'color',col(c_id+1,:),'LineWidth',1)
    quiver(obj.Position(i,1),obj.Position(i,2),obj.UV(i,1)*scale,obj.UV(i,2)*scale,'AutoScaleFactor',1,'MaxHeadSize',100,'color',col(c_id+1,:),'LineWidth',1)
     plot(obj.Position(i,1),obj.Position(i,2),'.k');
    %quiver(obj.Position(i,1),obj.Position(i,2),obj.UV(i,1),obj.UV(i,2),'AutoScale','off','MaxHeadSize',100,'color',col(c_id+1,:),'LineWidth',1)
end

axis([obj.GridParameters(1) obj.GridParameters(2), obj.GridParameters(3) obj.GridParameters(4)])
colormap(col)
c=colorbar('Ticks',[0:0.25:1],'TickLabels',[0 90 180 270 360]);
c.Label.String = 'Orientation [deg]';
end
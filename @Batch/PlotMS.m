function [] = PlotMS(obj)
figure('Name',strcat('Results after Mean Shift for batch->', num2str(obj.ID)));
plot(obj.Data(:,1),obj.Data(:,2),'.'); % raw data
grid on;
hold on;
plot(obj.ClustersMeans(:,1),obj.ClustersMeans(:,2),'rx'); % computed centers with mean shift
plot(obj.ClustersMeans_MLE(:,1),obj.ClustersMeans_MLE(:,2),'gx'); % computed centers with mean shift
colors=jet(max(obj.ClustersIDs));
for clustIDIter=1:max(obj.ClustersIDs)
    plot(obj.Data(obj.ClustersIDs==clustIDIter,1),obj.Data(obj.ClustersIDs==clustIDIter,2),'o','Color',colors(clustIDIter,:));
    legendInfo{clustIDIter}=['MS cluster ' num2str(clustIDIter)];
end % clustIDIter=1:max(clustID)
LLimit=min(obj.Data(:,2));
ULimit=max(obj.Data(:,2));
Padding = abs(ULimit-LLimit);
LLimit=LLimit-0.1*Padding;
ULimit=ULimit+0.1*Padding;
axis([0 2*pi LLimit ULimit])
legendInfo=['data','MS centers','MLE centers',legendInfo];
legend(legendInfo)
xlabel('orientation [rad]');
ylabel('amplitude [m/s]');
end
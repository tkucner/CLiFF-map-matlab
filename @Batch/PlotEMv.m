function []=PlotEMv(obj,resolution)
figure('Name',strcat('Results after Mean Shift for batch->', num2str(obj.ID)));
plot(obj.Data(:,1),obj.Data(:,2),'.'); % raw data
grid on;
hold on;
colormap jet
plot(obj.Mean(:,1),obj.Mean(:,2),'rx'); % computed centers with mean shift
LLimit=min(obj.Data(:,2));
ULimit=max(obj.Data(:,2));
Padding = abs(ULimit-LLimit);
LLimit=LLimit-0.1*Padding;
ULimit=ULimit+0.1*Padding;
[DD,VV] = meshgrid(0:resolution:2*pi,LLimit:resolution:ULimit);
[vs,ds] = size(DD);

ZZ=zeros(numel(DD),1);
[ms,~]=size(obj.Mean);
for i=1:ms
    for K=-obj.Wind:obj.Wind
        ZZ=ZZ+mvnpdf([DD(:),VV(:)],obj.Mean(i,:)+[2*pi*K 0],obj.Cov(:,:,i))*obj.P(i);
    end
end

ZZ=reshape(ZZ,length(LLimit:resolution:ULimit),length(0:resolution:2*pi));
contour(DD,VV,ZZ,40)
axis([0 2*pi LLimit ULimit])
legendInfo={'data','EM centers'};
legend(legendInfo)
xlabel('orientation [rad]');
ylabel('amplitude [m/s]');
end
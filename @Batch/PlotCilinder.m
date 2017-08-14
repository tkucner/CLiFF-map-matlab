function [] = PlotCilinder(obj,resolution)
figure('Name',strcat('Results after EM for batch->', num2str(obj.ID)));
grid on;
hold on;
colormap jet
LLimit=min(obj.Data(:,2));
ULimit=max(obj.Data(:,2));
Padding = abs(ULimit-LLimit);
LLimit=LLimit-0.1*Padding;
ULimit=ULimit+0.1*Padding;
[DD,VV] = meshgrid(0:resolution:2*pi,LLimit:resolution:ULimit);
[vs,ds] = size(DD);
ZZ=zeros(size(DD));
[ms,~]=size(obj.Mean);
for dd=1:ds
    for vv=1:vs
        for i=1:ms
            for K=-obj.Wind:obj.Wind
%                 ZZ(vv,dd)
%                 DD(vv,dd)
%                 VV(vv,dd)
%                 obj.Cov(:,:,i)
%                 obj.Mean(i,:)
%                 K
                ZZ(vv,dd)=ZZ(vv,dd)+obj.normal([DD(vv,dd),VV(vv,dd)],obj.Cov(:,:,i),obj.Mean(i,:),K);%*obj.P(i);
            end
        end
    end
end

[x, y,z] = pol2cart(DD,ZZ+1,VV);
surf(x, y,z,ZZ,'FaceColor','interp','EdgeColor','none');
%surf(x, y,z,ZZ,'EdgeColor','none');
alpha(.6)
end
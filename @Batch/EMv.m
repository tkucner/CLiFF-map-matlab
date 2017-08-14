function [obj] = EMv(obj)


        ClustersMeans=obj.ClustersMeans;
        ClusterCovariance=obj.ClusterCovariance;
  

obj.RemNaN=0;
obj.RemSmal=0;
[N,~]=size(ClustersMeans); % number of detected cluster
p=ones(N,1)*1/N; % all the distributions have equal weights
L=-obj.Wind:obj.Wind; % winding intervals
[DPCount,~]=size(obj.Data); % number of data points
r=zeros(DPCount,N,length(L)); % responsibility matrix initialized with zeros
c=zeros(size(ClusterCovariance));
for j=1:N % computing initial covairances from mean-shift
    f=diag(ClusterCovariance(:,:,j));
    f=floor(log10(f));
    c(:,:,j)=[10^(f(1)-1) 0; 0 10^(f(2)-1)];
end % j=1:N
m=ClustersMeans; % initial means
%[mSize,~]=size(m);
OldLogLikelihood=0;
delta = Inf;
epsilon =10^(-5);
delta_v=Inf(1,10);
MaxIteration=500;
Iteration=0;
while abs(delta)>epsilon && Iteration<MaxIteration%abs(sum(delta_v))>epsilon/2
    for j=1:N
        for l=1:length(L)
            r(:,j,l)=p(j)*mvnpdf(obj.Data,m(j,:)+[2*pi*L(l) 0], c(:,:,j));
        end % l=1:length(L)
    end % j=1:N
    r(r<eps)=0;
    r=r./repmat(sum(sum(r,3),2),1,N,length(L));
    r(isnan(r))=0;
    
    % M step
    %     for j=1:N
    %         m(j,:)=[0 0];
    %         for i=1:DPCount
    %             t=[0 0];
    %             for l=1:length(L)
    %                 t=t+(obj.Data(i,:)-2*pi*[L(l) 0])*r(i,j,l);
    %             end % l=1:length(L)
    %             m(j,:)=m(j,:)+t;
    %         end % i=1:DPCount
    %         m(j,:)=m(j,:)./sum(sum(r(:,j,:)));
    %     end % j=1:mSize
    
    for j=1:N
        m(j,:)=[0 0];
        t=zeros(DPCount,2);
        for l=1:length(L)

            t=t+(obj.Data-repmat(2*pi*[L(l) 0],DPCount,1)).*repmat(r(:,j,l),1,2);
        end % l=1:length(L)
        m(j,:)=sum(t,1);
        m(j,:)=m(j,:)./sum(sum(r(:,j,:)));
    end % j=1:mSize
    
    
    %     for j=1:N
    %         c(:,:,j)=zeros(2);
    %         for i=1:DPCount
    %             t=zeros(2);
    %             for l=1:length(L)
    %                 t=t+((obj.Data(i,:)-m(j,:)-2*pi*[L(l) 0]).'*(obj.Data(i,:)-m(j,:)-2*pi*[L(l) 0]))*r(i,j,l);
    %             end % l=1:length(L)
    %             c(:,:,j)=c(:,:,j)+t;
    %         end % i=1:DPCount
    %         c(:,:,j)=c(:,:,j)./sum(sum(r(:,j,:)));
    %     end % j=1:N
    
    
    for j=1:N
        c(:,:,j)=zeros(2);
        %for i=1:DPCount
        t=zeros(DPCount,2,2);
        for l=1:length(L)
            d_mod=obj.Data-repmat(m(j,:),DPCount,1)-2*pi*repmat([L(l) 0],DPCount,1);
            t(:,1,1)=t(:,1,1)+d_mod(:,1).^2.*r(:,j,l);
            t(:,2,2)=t(:,2,2)+d_mod(:,2).^2.*r(:,j,l);
            t(:,1,2)=t(:,1,2)+d_mod(:,1).*d_mod(:,2).*r(:,j,l);
            t(:,2,1)=t(:,1,2);
        end % l=1:length(L)
        c(:,:,j)=c(:,:,j)+squeeze(sum(t,1));
        %end % i=1:DPCount
        c(:,:,j)=c(:,:,j)./sum(sum(r(:,j,:)));
    end % j=1:N
    
    
    for j=1:N
        p(j)=sum(sum(r(:,j,:)))/DPCount;
    end % j=1:N
    
    %----------------------------------------------------------------------
    % Bias the covaraince matrix to avoid singularity
    for j=1:N
        [~,chol_f]=chol(c(:,:,j));
        if (chol_f~=0 && c(1,1,j)>10^(-10) && c(2,2,j)>10^(-10)) ||  (rcond(c(:,:,j))<10^(-10))
            %if c(1,1,j)>eps && c(2,2,j)>eps
            c(:,:,j)=c(:,:,j)+eye(2)*10^(-10);
        end
        %end
    end
    %---------------------------------------------------------------------
    % discarding clusters with too small covariance
    
    NewN=N;
    rem=zeros(1,N);
    rs=0;
    rn=0;
    for j=1:N
        if c(1,1,j)<10^(-9) || c(2,2,j)<10^(-9)
            rem(j)=1;
            NewN=NewN-1;
            obj.RemSmal=obj.RemSmal+1;
            rs=rs+1;
        end
        if isnan(rcond(c(:,:,j)))
            rem(j)=1;
            NewN=NewN-1;
            obj.RemNaN=obj.RemNaN+1;
            rn=rn+1;
        end
    end
    %------
    % Ridgeline analysis
    %     rem_over=obj.MergePick(m,c,p);
    %     rem=unique([rem rem_over.']);
    %     obj.RemRed=numel(rem)-rs-rn;
    %     NewN=NewN-(numel(rem)-rs-rn);
    %
    rem=logical(rem);
    if N~=NewN
        c(:,:,rem)=[];
        m(rem,:)=[];
        p(rem)=[];
        N=NewN;
        OldLogLikelihood=0;
        delta = Inf;
        %p=ones(N,1)*1/N; % all the distributions have equal weights
        p_sum=sum(p);
        for i=1:numel(p)
            p(i)=p(i)/p_sum;
        end
        r=zeros(DPCount,N,length(L)); % responsibility matrix initialized with zeros
    end
    
    %---------------------------------------------------------------------
    % discarding redundant clusters
    if isempty(m)
        return
    end
    m(:,1)=wrapTo2Pi(m(:,1));
    %----------------------------------------------------------------------
    LogLikelihood=0;
    %     for i=1:DPCount
    %         ll=0;
    %         for j=1:N
    %             for l=1:length(L)
    %                 ll=ll+p(j)*Batch.normal(obj.Data(i,:),c(:,:,j),m(j,:),L(l));
    %             end % l=1:length(L)
    %         end % j=1:N
    %         LogLikelihood=LogLikelihood+log(ll);
    %     end % i=1:DPCount
    
    ll=zeros(DPCount,1);
    for j=1:N
        for l=1:length(L)
            ll=ll+p(j)*mvnpdf(obj.Data,m(j,:)+[2*pi*L(l) 0],c(:,:,j));
        end % l=1:length(L)
    end % j=1:N
    LogLikelihood=sum(log(ll));
    
    
    delta=OldLogLikelihood-LogLikelihood;
    OldLogLikelihood=LogLikelihood;
    
    
    
    %----------------------------------------------------------------------
    %delta_v=[delta_v(2:end),delta]
    %abs(sum(delta_v))
    %     m
    %     c
    %     p
    Iteration=Iteration+1;
    %--------------------------------------------------------
    %     %after convergence prune redundant means % ridgelin anayssis
    %     if ~(abs(delta)>epsilon && Iteration<MaxIteration)%abs(sum(delta_v))>epsilon/2
    %         rem_over=obj.MergePick(m,c,p);
    %         rem=unique([rem rem_over.']);
    %         obj.RemRed=numel(rem);
    %         NewN=NewN-numel(rem);
    %         rem=logical(rem);
    %         if N~=NewN
    %             c(:,:,rem)=[];
    %             m(rem,:)=[];
    %             p(rem)=[];
    %             N=NewN;
    %             OldLogLikelihood=0;
    %             delta = Inf;
    %             %p=ones(N,1)*1/N; % all the distributions have equal weights
    %             p_sum=sum(p);
    %             for i=1:numel(p)
    %                 p(i)=p(i)/p_sum;
    %             end
    %             r=zeros(DPCount,N,length(L)); % responsibility matrix initialized with zeros
    %         end
    %         if isempty(m)
    %             return
    %         end
    %         m(:,1)=wrapTo2Pi(m(:,1));
    %
    %         LogLikelihood=0;
    %
    %         ll=zeros(DPCount,1);
    %         for j=1:N
    %             for l=1:length(L)
    %                 ll=ll+p(j)*mvnpdf(obj.Data,m(j,:)+[2*pi*L(l) 0],c(:,:,j));
    %             end % l=1:length(L)
    %         end % j=1:N
    %         LogLikelihood=sum(log(ll));
    %
    %
    %         delta=OldLogLikelihood-LogLikelihood;
    %         OldLogLikelihood=LogLikelihood;
    %
    %     end
    
    
    
    % after convergence prune redundant means % ridgelin anayssis
    if ~(abs(delta)>epsilon && Iteration<MaxIteration) && obj.ridgeline%abs(sum(delta_v))>epsilon/2 <- Check if we have converged
        fprintf('ridgeline\n');
        % find all distributions that are hiding in the slope of the others and
        % measure the distances between the picks
        dd=Inf(numel(p));
        for i=1:numel(p)
            for j=i+1:numel(p)
                %[a,x,dis,valm,PI]=RidgelineAnalysis(Mean(i,:),Cov(:,:,i),P(i),Mean(j,:),Cov(:,:,j),P(j));
                [~,~,~,~,PI]=Batch.RidgelineAnalysis(m(i,:),c(:,:,i),p(i),m(j,:),c(:,:,j),p(j));
                dPI=PI-p(i)/(p(i)+p(j));
                %dd(i,j)=numel(find(dPI(1:end-1).*dPI(2:end)<0));
                piks=numel(find(diff(dPI>=0)>0));
                if piks<=1
                    AD=Batch.DistanceWrapped(m(i,:),m(j,:));
                    dd(i,j)=AD;
                end
            end
        end
        while sum(sum(~isinf(dd)))~=0
            [~,I]=min(dd(:)); % find the picsk to merge which are very close
            [A,B]=ind2sub(size(dd),I); % find the indices
            
            N=N-1;
            if p(A)<p(B)
                dd(A,:)=[];
                dd(B,:)=Inf;
                dd(:,A)=[];
                dd(:,B)=Inf;
                [M,~] = Batch.WeightedMeanCS(m([A,B],:),p([A,B]));
                m(B,:)=M;
                %m(B,1)=wrapTo2Pi(m(B,1));
                c(:,:,B)=c(:,:,B)*(1+p(A));
                p(A)=[];
                m(A,:)=[];
                c(:,:,A)=[];
                OldLogLikelihood=0;
                %delta = Inf;
                %p=ones(N,1)*1/N; % all the distributions have equal weights
                p_sum=sum(p);
                for i=1:numel(p)
                    p(i)=p(i)/p_sum;
                end
                r=zeros(DPCount,N,length(L));
            else
                dd(B,:)=[];
                dd(A,:)=Inf;
                dd(:,B)=[];
                dd(:,A)=Inf;
                [M,~] = Batch.WeightedMeanCS(m([B,A],:),p([B,A]));
                m(A,:)=M;
                %m(A,1)=wrapTo2Pi(m(A,1));
                c(:,:,A)=c(:,:,A)*(1+p(B));
                p(B)=[];
                m(B,:)=[];
                c(:,:,B)=[];
                OldLogLikelihood=0;
                %delta = Inf;
                %p=ones(N,1)*1/N; % all the distributions have equal weights
                p_sum=sum(p);
                for i=1:numel(p)
                    p(i)=p(i)/p_sum;
                end
                r=zeros(DPCount,N,length(L));
                
            end
        end
        m(:,1)=wrapTo2Pi(m(:,1));
        %LogLikelihood=0;
        ll=zeros(DPCount,1);
        for j=1:N
            for l=1:length(L)
                ll=ll+p(j)*mvnpdf(obj.Data,m(j,:)+[2*pi*L(l) 0],c(:,:,j));
            end % l=1:length(L)
        end % j=1:N
        LogLikelihood=sum(log(ll));
        delta=OldLogLikelihood-LogLikelihood;
        OldLogLikelihood=LogLikelihood;
    end
end


obj.Mean=m;
obj.MapCluster_Direction=-1*ones(size(p));
obj.Cov=c;
obj.P=p;
end
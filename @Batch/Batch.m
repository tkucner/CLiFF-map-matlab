classdef Batch
    %BATCH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = public)
        %-input-----------------------------------------------------------%
        ID % id of the batch
        Pose % position of the batch in the cartesian cooridnate frame
        Data % Data points associated with the batch
        TrackID % id of the track to whichc data Point belongs
        Wind % winding number
        %-----------------------------------------------------------------%
        
        %-parameters for Mean Shift algorithm-----------------------------%
        Bandwidth % size of the smoothing kernel
        KernelType % type of used kernel for meanshift
        DistanceType % type of used distance
        StopFraction % stop condition, how small shift we classify as no shift
        ClustersIDs % Data points associated to clusters after running MeanSHift algorithm
        ClustersMeans % Means of the clusters after Mean Shift algorithm
        ClustersMeans_MLE % Means of the clusters after Mean Shift algorithm
        ClusterCovariance % Covariances of the clusters after Mean Shift algorithm
        %-----------------------------------------------------------------%
        
       
        %-output----------------------------------------------------------%
        Mean % Means of the data
        Cov % Covariances
        P % mixing factors
        Observed % Time for how long the point was observed
        p % weight representing our trunst in the motion pattern
        scale_p % weight representing our trunst in the motion pattern scaled so the longest observation time get score 1
        q % weight representing our trust in the location
        t % weight representing our trust in reconstructed distibution
        R % The mena resultatn length - metric shwoing how much the data is spread for each distirbution
        BIC % Bayesina information criterion informing how well the distribution fits the data
        AIC % Akaike information criterion informing how well the distribution fits the data
        %-----------------------------------------------------------------%
        
        %-cluistering-----------------------------------------------------%
        MapCluster_Direction % labels to which cluster in the map belongs the distibutions
        %-----------------------------------------------------------------%
        
        %-Diagnostics-----------------------------------------------------%
        RemNaN % number of removed nans in EM because of NaN
        RemSmal % number of remoced clusters in EM because of small cov
        RemRed % number of removed cluster in EM becasue of redundacy
        ridgeline % flag if to use ridgeline analysis
        %-----------------------------------------------------------------%
        
       
    end
    
    properties (Constant)
        ThetaLo = 0 ; % lower limit of the speed
        ThetaHi = 2*pi; % higher limit of the speed
        RhoLo = 0;
        RhoHi = 10;
    end
    
    methods (Static = false)
   
        function obj = Batch()
        end
        
        function obj = SetParameters(obj,id,pose,data,wind,ridgeline)
            obj.ID = id;
            obj.Pose = pose;
            if nargin==5
                %obj.ridgeline=true;
                obj.ridgeline=false;
            elseif nargin==6
                obj.ridgeline=ridgeline;
            end
            %obj.Data = data;
            if ~isempty(data)
                obj=obj.AddData(data);
            end
            obj.Wind = wind;
        end
        
        function obj = AddData(obj,data)
            data=data(data(:,2)~=0,:); % removing static objects from measurements
            data(:,1)=wrapTo2Pi(data(:,1));
            obj.Data = data;
        end
        
       
        
        function obj = SetParametersMeanShift(obj,stopfraction,bandwidth)
           
            obj.StopFraction = stopfraction;
            if bandwidth == -1 
                [~,C]=obj.mle(obj.Data);
                sigma_c=C(1,1);
                sigma_l=C(2,2);
                h_l=((4*sigma_l^5)/(3*length(obj.Data)))^(1/5);
                h_c=((4*sigma_c^5)/(3*length(obj.Data)))^(1/5);
                if h_l<eps || h_c<eps
                    h_l=h_l+10^-10;
                    h_c=h_c+10^-10;
                end
                obj.Bandwidth=[h_c 0;0 h_l];
            else
                obj.Bandwidth = bandwidth;
            end % nargin==3 && obj.KernelType == MeanSiftKernel.Gaussian
        end
        
        function prob = Probability(obj,Velocity)
            prob=0;
            [ms,~]=size(obj.Mean);
            for i=1:ms
                for k=-obj.Wind:obj.Wind
                    prob=prob+Batch.normal(Velocity,obj.Cov(:,:,i),obj.Mean(i,:),k)*obj.P(i);
                end
            end
        end
        
       
        
        
        [obj,shiftData] = MeanShift2Dv(obj)
        obj = EMv(obj)
       
        
        
        [] = PlotMS(obj)
        
        [] = PlotMS_jitter(obj)
        
        [] = PlotEMv(obj,resolution)
        [] = PlotPolar(obj,resolution)
        [] = PlotCilinder(obj,resolution)
        [] = PlotMSKernel(obj)
        
        
        
        [samples]=Sample(obj,N)
        [samples]=SampleRaw(obj,N)
        [samples]=SampleScaledP(obj,N)
        obj = ScoreFit(obj)
        
        
    end
    methods (Static = true)
        
        
        
        
        
        function r = DistanceWrapped(P1,P2)
            %DISTANCEWRAPPED computes distance between two points in circular
            %linear space
            switch nargin
                case 2
                    AD=abs(wrapToPi(P1(1)-P2(1)));
                    LD=abs(P1(2)-P2(2));
                    r=sqrt(AD^2+LD^2);
                otherwise
                    error('Function cl_dist requires 2 arguments')
            end
        end
        
        function r = DistanceDisjoint(P1,P2)
            %DISTANCEWRAPPED computes distance between two points in circular
            %linear space as sugested by Roy and Pauri
            switch nargin
                case 2
                    AD=abs(wrapToPi(P1(1)-P2(1)));
                    LD=abs(P1(2)-P2(2));
                    r=[AD,LD];
                otherwise
                    error('Function cl_dist requires 2 arguments')
            end
        end
        
        function r = DistanceDisjointVec(P1,P2)
            %DISTANCEWRAPPED computes distance between two points in circular
            %linear space as sugested by Roy and Pauri
            switch nargin
                case 2
                    AD=abs(wrapToPi(P1(:,1)-P2(:,1)));
                    LD=abs(P1(:,2)-P2(:,2));
                    r=[AD,LD];
                otherwise
                    error('Function cl_dist requires 2 arguments')
            end
        end
        
        
        
        function [M,C] = mle_complex(X)
            % MLE - esitmates the mean and covariance for given circualr-linear data
            % X - input, first row should be circular the second should be linear
            if ndims(X)==2
                z=exp(1i*X(:,1));
                cr_m=mod(angle(sum(z)/length(z)),2*pi);
                l_m=mean(X(:,2));
                M=[cr_m,l_m];
                c_m=sum(z)/length(z);
                r=c_m*conj(c_m);
                ee=length(z)/(length(z)-1)*(r-1/length(z));
                V_c=log(1/ee);
                V_l=cov(X(:,2));
                c=0;
                for j=1:length(X(:,1))
                    c=c+(X(j,2)-l_m)*(wrapToPi(X(j,1)-cr_m));
                end
                c=c/length(X(:,2));
                C=[V_c c; c V_l];
            else
                error('The arguments should be 2 dimensional')
            end
        end
        
        
        
        function [M,C] = mle(X)
            % MLE - esitmates the mean and covariance for given circualr-linear data
            % X - input, first row should be circular the second should be linear
            if ndims(X)==2
                C=mean(cos(X(:,1)));
                S=mean(sin(X(:,1)));
                R=sqrt(C^2+S^2);
                if C>=0
                    cr_m=atan(S/C);
                else
                    cr_m=atan(S/C)+pi;
                end
                l_m=mean(X(:,2));
                M=[wrapTo2Pi(cr_m),l_m];
                std=sqrt(-2*log(R));
                V_c=std^2;
                V_l=cov(X(:,2));
                c=0;
                for j=1:length(X(:,1))
                    c=c+(X(j,2)-l_m)*(wrapToPi(X(j,1)-cr_m));
                end
                c=c/(length(X(:,2))-1);
                C=[V_c c; c V_l];
            else
                error('The arguments should be 2 dimensional')
            end
        end
        
        function [M] = Centroid(X)
            % MLE - esitmates the mean and covariance for given circualr-linear data
            % X - input, first row should be circular the second should be linear
            if ndims(X)==2
                C=mean(cos(X(:,1)));
                S=mean(sin(X(:,1)));
                R=sqrt(C^2+S^2);
                if C>=0
                    cr_m=atan(S/C);
                else
                    cr_m=atan(S/C)+pi;
                end
                l_m=mean(X(:,2));
                M=[wrapTo2Pi(cr_m),l_m];
                
            else
                error('The arguments should be 2 dimensional')
            end
        end
        
        
        function M = WeightedMeanComplex(X,W)
            z=exp(1i*X(:,1));
            z=z.*W;
            cr_m=mod(angle(sum(z)/sum(W)),2*pi);
            l_m=sum(X(:,2).*W)/sum(W);
            M=[cr_m,l_m];
        end
        
        function [M,R] = WeightedMeanCS(X,W)
            C=sum(cos(X(:,1)).*W)/sum(W);
            S=sum(sin(X(:,1)).*W)/sum(W);
            R=sqrt(C^2+S^2);
            if C>=0
                cr_m=atan(S/C);
            elseif C<0
                cr_m=atan(S/C)+pi;
            end
            l_m=sum(X(:,2).*W)/sum(W);
            M=[wrapToPi(cr_m),l_m];
        end
        
        function M = WeightedMeanUV(X,W)
            [U,V]=pol2cart(X(:,1),X(:,2));
            Um=sum(U.*W)/sum(W);
            Vm=sum(V.*W)/sum(W);
            [cr_m,l_m]=cart2pol(Um,Vm);
            M=[cr_m,l_m];
        end
        
        function M = MeanComplex(X)
            z=exp(1i*X(:,1));
            cr_m=mod(angle(sum(z)/length(z)),2*pi);
            l_m=mean(X(:,2));
            M=[cr_m,l_m];
        end
        
        function [M,R] = MeanCS(X)
            C=mean(cos(X(:,1)));
            S=mean(sin(X(:,1)));
            R=sqrt(C^2+S^2);
            if C>=0
                cr_m=atan(S/C);
            elseif C<0
                cr_m=atan(S/C)+pi;
            end
            l_m=sum(X(:,2).*W)/sum(W);
            M=[cr_m,l_m];
        end
        
        function M = MeanUV(X)
            [U,V]=pol2cart(X(:,1),X(:,2));
            Um=mean(U);
            Vm=mean(V);
            [cr_m,l_m]=cart2pol(Um,Vm);
            M=[cr_m,l_m];
        end
        
        function [p] = normal(x,sigma,mu,k)
            %NORMAL Summary of this function goes here
            %   Detailed explanation goes here
            %x=x.';
            %mu=mu.';

            p=mvnpdf(x,mu+[2*pi*k 0],sigma);
            %p=1/(sqrt((2*pi)^2*norm(sigma)))*exp(-0.5*(x-mu-K).'*sigma^(-1)*(x-mu-K));
        end
        
        function y = mvgmmrnd(mu,sigma,p,n)
            %MVGMMRND Random vectors from a mixture of multivariate normals.
            % MU is an M-by-D matrix of means for the M component normals
            % SIGMA is a D-by-D-by-M array of covariance matrices for the
            % M component normals.
            % P is an M-by-1 vector of component mixing probabilities.
            % N is the desired number of random vectors.
            [M,d] = size(mu);
            % randomly pick from the components
            [~,compon] = histc(rand(n,1), [0; cumsum(p(:))./sum(p)]);
            % generate random vectors from the selected components with a
            % "stacked" matrix multiply
            for i = 1:M
                Rt(i,:,:) = chol(sigma(:,:,i)); % holds the transposed cholesky factors
            end
            Z = repmat(randn(n,d), [1,1,d]);
            if n==1
                y = squeeze(sum(Z.*Rt(compon,:,:),2)).' + mu(compon,:);
            else
                y = squeeze(sum(Z.*Rt(compon,:,:),2)) + mu(compon,:);
            end
        end
        
        function [a,x,dis,val,PI]= RidgelineAnalysis(mu_1,sigma_1,pi_1,mu_2,sigma_2,pi_2)
            mu_1=mu_1.';
            mu_2=mu_2.';
            % unwrap the distributions on the surface in a convinient way
            AD=wrapToPi(mu_1(1)-mu_2(1));
            mu_1(1)=0;
            mu_2(1)=AD;
            % normalize pi for 2 element GMM
            pi_norm=pi_1+pi_2;
            pi_1=pi_1/pi_norm;
            pi_2=pi_2/pi_norm;
            % compute points for ridge line
            a=0:0.01:1;
            for i=1:numel(a)
                x(:,i)=inv((1-a(i)).*inv(sigma_1)+a(i).*inv(sigma_2))* ...
                    ((1-a(i)).*inv(sigma_1)*mu_1+a(i).*inv(sigma_2)*mu_2);
            end
            %compute distance between points for normalised plot
            dis=0;
            for i=2:numel(a)
                X=[(x(:,i-1)-mu_1).';(x(:,i)-mu_1).'];
                dis=[dis dis(end)+pdist(X,'euclidean')];
            end
            %compute the densities along the ridgeline
            val=[];
            %compute profile
            for i=1:numel(a)
                val=[val pi_1*mvnpdf(x(:,i),mu_1,sigma_1)+pi_2*mvnpdf(x(:,i),mu_2,sigma_2)];
            end
            %compute PI
            for i=1:numel(a) % ugly hack fix it later
                PI(i)=1/(1+((1-a(i))/(a(i)))*(mvnpdf(x(:,i),mu_2,sigma_2)/mvnpdf(x(:,i),mu_1,sigma_1)));
            end
        end
        
        function me = MergePick(Mean,Cov,P)
            me=[];
            dd=Inf(numel(P));
            for i=1:numel(P)
                for j=i+1:numel(P)
                    %[a,x,dis,valm,PI]=RidgelineAnalysis(Mean(i,:),Cov(:,:,i),P(i),Mean(j,:),Cov(:,:,j),P(j));
                    [~,~,~,~,PI]=Batch.RidgelineAnalysis(Mean(i,:),Cov(:,:,i),P(i),Mean(j,:),Cov(:,:,j),P(j));
                    dPI=PI-P(i);
                    %dd(i,j)=numel(find(dPI(1:end-1).*dPI(2:end)<0));
                    piks=numel(find(diff(dPI>=0)>0));
                    if piks<=1
                        AD=Batch.DistanceWrapped(Mean(i,:),Mean(j,:));
                        dd(i,j)=AD;
                    end
                end
                if sum(~isinf(dd(i,:)))>0
                    [~,j]=min(dd(i,:));
                    if P(i)<P(j)
                        me=[me;i];
                    else
                        me=[me;j];
                    end
                    
                end
            end
            
        end
    end
end
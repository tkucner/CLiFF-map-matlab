classdef DynamicMap
    %DYNAMICMAP - class for handling map of dynamic within the environemtn
    %using irectional statistics.
    
    properties (SetAccess = public)
        LocationID % location ID this will be solved with batches
        TrackID % ID of the track
        TimeStamp % time stamp of observation
        Position % Position of measurment
        UV % data projected on X and Y axis
        ThetaRho % data as orientation and magnitude
        Batches % batches per location
        BatchesSparse % batches fro sparse data
        SparseP % trust factor based on the time of motion
        ScaleSparseP % trust factor based on the time of motion scaled so the longest time is 1
        SparseQ % trust factor based on the time of observation
        GridParameters % measurement discretisation [xmin xmax ymin ymax step]
        Grid % position of actual nodes
        File % file with data to load
        Radious % radious used to establish if object is within a cluster
        Image % backgraound image for the plot
        
        Wind % winding number
        
        TrustHistogramT % "histogram" showing how much we trust our reconstructed map
        TrustHistogramP % "histogram" showing how much we trust in motion model
        TrustHistogramScaleP % "histogram" showing how much we trust in motion model
        TrustHistogramQ % "histogram" showing how much we trust in the locations
        
        TrustHistogramSparse % histogram showing how many observations are supporting our model
        % Prototype % Prototype of parameters for batch
        
        %-Directional Clustering-------------------------------------------
        ClusterIDs % IDs of the clusters within the map
        ClusterMeans % mean orientation
        ClusterMembers % Batches belonging to a given cluster
        %------------------------------------------------------------------
    end
    methods(Static = false)
        function obj = DynamicMap()
        end
        
        function obj = SetParameters(obj,radious,xmin,xmax,ymin,ymax,step,wind)
            obj.Radious=radious;
            
            obj.GridParameters=[xmin,xmax,ymin,ymax,step];
            obj.Wind=wind;
            [X_eg,Y_eg] = meshgrid(xmin:step:xmax,ymin:step:ymax);
            obj.TrustHistogramT=zeros(size(X_eg));
            obj.TrustHistogramP=zeros(size(X_eg));
            obj.TrustHistogramQ=zeros(size(X_eg));
            obj.Grid=[reshape(X_eg,[],1),reshape(Y_eg,[],1)];
            obj.Batches=Batch.empty;
            % obj.Prototype=Batch(0,[0,0],[],wind);
            for j=1:length(obj.Grid)
                obj.Batches(j)=Batch();
                obj.Batches(j)=obj.Batches(j).SetParameters(j, obj.Grid(j,:),[],wind,false);        
            end
        end
        
        function obj = SetImage(obj,image)
            obj.Image = imread(image);
        end
        
        obj = SplitToLocations(obj)
        obj = ProcessBatches(obj)
        obj = ProcessBatchesSparse(obj)
        [] = SaveXML(obj,FileName)
        [] = PlotUVDirection(obj,scale)
        [] = PlotMapSparseDirection(obj,filterTreshold,scale)
        [] = PlotMapDirection(obj,filterTreshold,scale)
    end
end
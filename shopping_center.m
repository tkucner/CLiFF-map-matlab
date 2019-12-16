clear all;
close all
clc

% File conatining input data.
FILE='/home/ksatyaki/workspace/DATA/UTBM/training_atc_point1.txt';
%PATH='Data';
full_path=FILE;
% Load input data to matrix.
DATA=csvread(full_path);

%time [ms], id, x [mm], y [mm], z [mm], velocity [mm/s], angle of motion [rad], facing angle [rad]
%DATA = DATA(1:1000,:);
% Create object of dynamic map.
DM=DynamicMap();
% Load time stamps of measurements
DM.TimeStamp=DATA(:,1); 
% Load cooridantes of the measuremnts
DM.Position=DATA(:,3:4)/1000; 

DM.ThetaRho=[DATA(:,7), DATA(:,6)/1000];
[U,V] = pol2cart(DATA(:,7), DATA(:,6)/1000);
DM.UV = [U,V];

% In this example the measurments are disitributed through the
% environmentnt, in order to build a map we need to define the boundries.
% In the following 4 lines a bounding box for the data is defined.
min_x=-10.0;
max_x=3.0;
min_y=0.0;
max_y=14.0;


DM.File=FILE;
% Setting parameters for the map
DM=DM.SetParameters(0.5,min_x,max_x,min_y,max_y,0.5,1);
% Split data into batches
DM=DM.SplitToLocations();
% % Compute the parameters of the distribution
DM=DM.ProcessBatches();

% Compute PQ values
localBatches = DM.Batches;

for ib=1:numel(localBatches)
    
  [m,~]=size(localBatches(ib).Data);
  q = m/max_observations;
  p = 1.0;
  
  localBatches(ib).q=q;
  localBatches(ib).p=p;
  DM.SparseP(ib)=P;
  DM.SparseQ(ib)=Q;
  dm_2.SparseP(ib)=p;
  
  [I,J]=ind2sub(size(dm_2.TrustHistogramQ),ib);
  dm_2.TrustHistogramQ(I,J)=q;
end


%% Plot the color-coded input data
%DM.PlotUVDirection(2)

DM = DM.SetImage('/home/ksatyaki/workspace/cpp_ws/src/ompl_planners/ompl_planners_ros/maps/utbm.pgm');
%Plot resulting distribution
DM.PlotMapDirection(0.1,0.2)
%% SAVE XML
 % DONT FORGET THIS STEP
DM.SaveXML('utbm_point1.xml')

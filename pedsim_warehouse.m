clear all;
close all
clc

% File conatining input data.
% FILE='/home/ksatyaki/workspace/DATA/UTBM/training_atc_point1.txt'; UTBM
FILE='/home/ksatyaki/workspace/DATA/pedsim/warehouse/training/run1.csv';
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
% min_x=-10.0;
% max_x=3.0;
% min_y=0.0;
% max_y=14.0;

% ATC
min_x=0.0;
max_x=35.0;
min_y=0.0;
max_y=30.0;


DM.File=FILE;
% Setting parameters for the map
DM=DM.SetParameters(1.0,min_x,max_x,min_y,max_y,1.0,1);
% Split data into batches
DM=DM.SplitToLocations();
% % Compute the parameters of the distribution
DM=DM.ProcessBatches();

% Compute PQ values
localBatches = DM.Batches;

%%

max_observations = 0;
for ib=1:numel(localBatches)
  [m,~] = size(localBatches(ib).Data);
  if(max_observations < m)
    max_observations = m;
  end
end

for ib=1:numel(localBatches)
    
  [m,~]=size(localBatches(ib).Data);
  q = m/max_observations;
  p = 1.0;
  
  localBatches(ib).q=q;
  localBatches(ib).p=p;
  DM.SparseP(ib)=p;
  DM.SparseQ(ib)=p;
  
  [I,J]=ind2sub(size(DM.TrustHistogramQ),ib);
  DM.TrustHistogramQ(I,J)=q;
end

DM.Batches = localBatches

%% SAVE XML
 % DONT FORGET THIS STEP
DM.SaveXML('pedsim_warehouse.xml')

%% Plot the color-coded input data
%DM.PlotUVDirection(2)

DM = DM.SetImage('/home/ksatyaki/workspace/cpp_ws/src/ompl_planners/ompl_planners_ros/maps/atc_white.pgm');
%Plot resulting distribution
DM.PlotMapDirection(0.1,0.2)

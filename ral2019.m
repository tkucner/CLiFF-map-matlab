clear all;
close all
clc

% File conatining input data.
FILE='/home/ksatyaki/workspace/DATA/New/training_dataset_point1.txt';
%PATH='Data';
full_path=FILE;
% Load input data to matrix.
DATA=readmatrix(full_path, 'Delimiter', ' ');

fprintf("Read matrix file successfully.\n");

% unix time, x, y, vel_x, vel_y, speed, heading, id
% ALL UNITS SI
% DATA = DATA(1:1000,:);
% Create object of dynamic map.
DM=DynamicMap();
% Load time stamps of measurements
DM.TimeStamp=DATA(:,1); 
% Load cooridantes of the measuremnts
DM.Position=DATA(:,2:3); 

DM.ThetaRho=[DATA(:,7), DATA(:,6)];
DM.UV = DATA(:,4:5);

% In this example the measurments are disitributed through the
% environmentnt, in order to build a map we need to define the boundries.
% In the following 4 lines a bounding box for the data is defined.
min_x=-10.0;
max_x=2.0;
min_y=0.0;
max_y=18.0;


DM.File=FILE;
% Setting parameters for the map
DM=DM.SetParameters(0.5,min_x,max_x,min_y,max_y,0.5,1);

fprintf("Splitting data into batches.\n");
% Split data into batches
DM=DM.SplitToLocations();

fprintf("Beginning to process batches.\n");

% Compute the parameters of the distribution
DM=DM.ProcessBatches();
%%
fprintf("Computing P and Q\n");
max_q = 0.0;
for i=1:numel(DM.Batches)
  DM.Batches(i).p = 1.0;
  DM.Batches(i).q = size(DM.Batches(i).Data,1) / size(DATA,1);
  if(DM.Batches(i).q > max_q)
    max_q = DM.Batches(i).q
  end
end

for i = 1:numel(DM.Batches)
  DM.Batches(i).q = DM.Batches(i).q / max_q;
end
%% Plot the color-coded input data
%DM.PlotUVDirection(2)

%Plot resulting distribution
% DM.PlotMapDirection(0.1,0.2)

DM.SaveXML('ral2019_normalized_q.xml')

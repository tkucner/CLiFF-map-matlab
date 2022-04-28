clear all;
close all
clc

% File conatining input data.
FILE='/home/ksatyaki/workspace/matlab/CLIFF-map-matlab/experiment_csv.txt';
%PATH='Data';
full_path=FILE;
% Load input data to struct.
DATA=readmatrix(full_path, 'Delimiter', ' ');

%%

%time [ms], id, x [mm], y [mm], z [mm], velocity [mm/s], angle of motion [rad], facing angle [rad]
%DATA = DATA(1:1000,:);
% Create object of dynamic map.
DM=DynamicMap();
% Load cooridantes of the measuremnts
DM.Position=DATA(:,1:2); 

DM.ThetaRho=[DATA(:,5), DATA(:,6)];
DM.UV = [DATA(:,3), DATA(:,4)];

% In this example the measurments are disitributed through the
% environmentnt, in order to build a map we need to define the boundries.
% In the following 4 lines a bounding box for the data is defined.
min_x=-20.0;
max_x=30.0;
min_y=-20.0;
max_y=30.0;


DM.File=FILE;
% Setting parameters for the map
DM=DM.SetParameters(0.25,min_x,max_x,min_y,max_y,0.5,0);
% Split data into batches
DM=DM.SplitToLocations();
% % Compute the parameters of the distribution
DM=DM.ProcessBatches();
%% Plot the color-coded input data
%DM.PlotUVDirection(2)
DM = DM.SetImage('/home/ksatyaki/workspace/DATA/RudenkoHumanTracking/exp1_flipped.pgm');
%Plot resulting distribution
DM.PlotMapDirection(0.1,0.2)
%DM.PlotMap(0.2)
DM.SaveXML('robotlab_real_exp1-1_point5.xml')


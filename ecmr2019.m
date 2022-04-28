clear all;
close all
clc

% File conatining input data.
FILE='/home/ksatyaki/workspace/matlab/Data/two_weeks_days_nights_weekends_with_angles_plus_reversed.txt';
%PATH='Data';
full_path=FILE;

DATA = readmatrix(FILE, 'Delimiter',' ');
% Load input data to matrix.

% The data is in the following format:
% time x y theta speed flag

% Create object of dynamic map.
DM=DynamicMap();
% Load time stamps of measurements
DM.TimeStamp=DATA(:,1); 
% Load cooridantes of the measuremnts
DM.Position=DATA(:,2:3); 

DM.ThetaRho=[DATA(:,4), DATA(:,5)];
[U,V] = pol2cart(DATA(:,4), DATA(:,5));
DM.UV = [U,V];

% In this example the measurments are disitributed through the
% environmentnt, in order to build a map we need to define the boundries.
% In the following 4 lines a bounding box for the data is defined.
min_x=-8;
max_x=23;
min_y=-6;
max_y=16;


DM.File=FILE;
% Setting parameters for the map
DM=DM.SetParameters(1,min_x,max_x,min_y,max_y,1,0);
% Split data into batches
DM=DM.SplitToLocations();
% % Compute the parameters of the distribution
DM=DM.ProcessBatches();
%% Plot the color-coded input data
%DM.PlotUVDirection(2)

%Plot resulting distribution
DM.PlotMapDirection(0.1,1.5)

DM.SaveXML('ecmr_map.xml')

% ======================================================================= %
% AIR FLOW EXAMPLE                                                        %
%                                                                         %
% The following example shows how to construct CLiFF-map for cases where  %
% the measurements are clustered in a set of locations. As an example a   %
% data set contianign air flow measurments is used.                       %
%                                                                         %
%                                                                         %
% Author: Tomasz Kuncer                                                   %
% e-mail: tomasz.kucner@oru.se                                            %
% ======================================================================= %


clear all;

% File conatining input data.
FILE='air_1.csv';
PATH='Data';
full_path=fullfile(PATH,FILE);
% Load input data to matrix.
DATA=csvread(full_path);
% Create object of dynamic map.
DM=DynamicMap();
% Load time stamps of measurements
DM.TimeStamp=DATA(:,1); 
% Load cooridantes of the measuremnts
DM.Position=DATA(:,2:3); 
% Load velocity measurements in Kartesian cooridnate frame
DM.UV=DATA(:,4:5); 
% Convert measurements to the polar cooridante frame
[TH,R]=cart2pol(DM.UV(:,1),DM.UV(:,2)); 
DM.ThetaRho=[TH,R]; 
% Load locations IDs
DM.LocationID=DATA(:,6); 
DM.File=FILE;
% Setting parameters for the map
DM=DM.SetParameters(0.25,0,4,0,3.5,0.5,1);
% Split data into batches
DM=DM.SplitToLocations();
% Compute the parameters of the distribution
DM=DM.ProcessBatchesSparse();
% Plot the color-coded input data
DM.PlotUVDirection(2)
%Print resulting distribution
DM.PlotMapSparseDirection(0,2)
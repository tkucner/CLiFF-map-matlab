% ======================================================================= %
% PEDESTRIAN FLOW EXAMPLE                                                 %
%                                                                         %
% The following example shows how to construct CLiFF-map for cases where  %
% the measurements are spread over the map. Such cases are when the robot %
% is collecting velocity measurments of walkignpeople                     %
%                                                                         %
%                                                                         %
% Author: Tomasz Kuncer                                                   %
% e-mail: tomasz.kucner@oru.se                                            %
% ======================================================================= %


%clear all;

% File conatining input data.
% FILE='';
% PATH='Data';
% full_path=fullfile(PATH,FILE);
% % Load input data to matrix.
% DATA=csvread(full_path);
% Create object of dynamic map.
DM=DynamicMap();
% Load time stamps of measurements
% DM.TimeStamp=DATA(:,2); 
% Load cooridantes of the measuremnts
DM.Position=[x,y];
DM.ThetaRho=[th,s]; 
% In this example the measurments are disitributed through the
% environmentnt, in order to build a map we need to define the boundries.
% In the following 4 lines a bounding box for the data is defined.
min_x=-4.8;
max_x=4.8;
min_y=-4.8;
max_y=4.8;

%DM.File=FILE;
% Setting parameters for the map
DM=DM.SetParameters(1,min_x,max_x,min_y,max_y,1,1);
% Split data into batches
DM=DM.SplitToLocations();
% % Compute the parameters of the distribution
DM=DM.ProcessBatches();
% % Plot the color-coded input data
DM.PlotUVDirection(2)
%Plot resulting distribution
DM.PlotMapDirection(0,1.5)

DM.SaveXML('pedestrian_map.xml')
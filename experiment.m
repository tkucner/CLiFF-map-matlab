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


clear all;

% File conatining input data.
FILE='maze_s.csv';
PATH='/home/chitt/workspace/DATA';
full_path=fullfile(PATH,FILE);
% Load input data to matrix.
DATA=csvread(full_path, 1, 0);
% Create object of dynamic map.
DM=DynamicMap();
% Load time stamps of measurements
DM.TimeStamp=DATA(:,2); 
% Load cooridantes of the measuremnts
DM.Position=DATA(:,3:4); 
% Load velocity measurements in Kartesian cooridnate frame
DM.UV=DATA(:,5:6); 
% Convert measurements to the polar cooridante frame
[TH,R]=cart2pol(DM.UV(:,1),DM.UV(:,2)); 
DM.ThetaRho=[TH,R]; 
% In this example the measurments are disitributed through the
% environmentnt, in order to build a map we need to define the boundries.
% In the following 4 lines a bounding box for the data is defined.
resolution = 1.0;
min_x= 0.0;
min_y= 0.0;
max_x= 120.0;
max_y= 95.0;

DM.File=FILE;
% Setting parameters for the map
DM=DM.SetParameters(1,min_x,max_x-1,min_y,max_y-1,resolution,0);
% Split data into batches
DM=DM.SplitToLocations();
DM=DM.SetImage('/home/chitt/workspace/cpp_ws/src/smp_ros/maps/maze.pgm');

% Here we get the max observations over all the cells. 
for ib=1:numel(DM.Batches)
  m = size(unique(DM.Batches(ib).TimeStamps),1);
  % We have one tick for every track. So we take total ticks / total
  % tracks.
  DM.SparseQ(ib)=m/(size(DATA,1)/80);
  DM.SparseP(ib)=1.0;
  
  DM.Batches(ib).q=m/(size(DATA,1)/80);
  DM.Batches(ib).p(ib)=1.0;
end

% % Compute the parameters of the distribution
DM=DM.ProcessBatches();
% % Plot the color-coded input data
DM.PlotUVDirection(2)
%Plot resulting distribution
DM.PlotMapDirection(0,1.5)

DM.SaveXML('maze.xml')
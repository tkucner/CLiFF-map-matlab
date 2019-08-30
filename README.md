# CLiFF-map-matlab
The MATLAB implementation of CLiFF-map a dynamic mapping algorithm

## What is CLiFF-map?
CLiFF-map stands for Circular-Linear Flow Field map. The idea behind this representation is that you can provide a set of velocity measurements (from human tracker or anemometer) and then build a map representing the motion patterns all over the environment. If you are interested in details please read my RA-L paper[1].

## What can I find in this repository?
The repositroy can be dived into two parts.
### @Batch & @DynamicMap
In these two folders you will find classess and functions providing the functionalities of CLiFF-map.
* *@Batch* contains code necessarry to build a model of flow for one single location.
* *@DynamicMap* provides functionalities necessarry to build a complete map.

### Air flow & Pedestrian flow examples
The repository is accompanied by two examples showing how to use CLiFF-map
* *air_flow.m* - This example shows how to build a sparse map out of a set of observations collected in a set of locations. In order to show the capabilities of the method I have used a wind measurements obtained by Yuta Wada, Marco Trincavelli, Yuichiro Fukazawa and Hiroshi Ishida[2].
![alt text](https://github.com/tkucner/CLiFF-map-matlab/blob/master/img/wind_in.png "Visualisation of input data")
![alt text](https://github.com/tkucner/CLiFF-map-matlab/blob/master/img/wind_out.png "Visualisation of output map")
* *people_flow.m* - This example shows ho to build CLiFF-map if the observations are randomy distibuted all over the map and it is necessarry to discretise them so they will fit a refualr grid. This example uses data set generated with ped-sim (https://github.com/srl-freiburg/pedsim_ros).
![alt text](https://github.com/tkucner/CLiFF-map-matlab/blob/master/img/pedestrian_in.png "Visualisation of input data")
![alt text](https://github.com/tkucner/CLiFF-map-matlab/blob/master/img/pedestrian_out.png "Visualisation of output map")

## When there will be more code published?
Soon! 

Currently I am organising my research code and preparing it for release. I will try to do it as fast as possible.

## I have so many questions!
Please feel free to send me an e-mail!
tomasz.kucner@oru.se

or Chittaranjan:
ksatyaki@gmail.com

## MATLAB dependency for nice-looking ProgressBar
MatlabProgressBar: https://github.com/JAAdrian/MatlabProgressBar

or here: https://de.mathworks.com/matlabcentral/fileexchange/57895-matlabprogressbar

## References
[1] T. P. Kucner, M. Magnusson, E. Schaffernicht, V. H. Bennetts and A. J. Lilienthal, "Enabling Flow Awareness for Mobile Robots in Partially Observable Environments," in IEEE Robotics and Automation Letters, vol. 2, no. 2, pp. 1093-1100, April 2017.
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7835155&isnumber=7797562

[2] Yuta Wada, Marco Trincavelli, Yuichiro Fukazawa and Hiroshi Ishida, Collecting a Database for Studying Gas Distribution Mapping and Gas Source Localization with Mobile Robots. International Conference on Advanced Mechatronics (ICAM), 2010.
URL: http://aass.oru.se/Research/Learning/publications/2010/Wada_etal_2010_ICAM10_Collecting_a_Database_for_Studying_Gas_Distribution_Mapping_and_Gas_Source_Localization_with_Mobile_Robots.pdf

# Unscented Kalman Filter Project Starter Code
The aim of this project is to combine LIDAR and RADAR data to predict vehicles/pedestrians
using unscented Kalman Filter.
The implementation follows the Udacity unscented kalaman filter sensor fusion algorithm.

The main.cpp file is designed to run with Udacity simulator only. Hence, it does not support
a standalone execution. However, the algorithm should work as a standalone program.

Steps to execute.
1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./ExtendedKF

- Initial mean values of velocity, turn rate, and change in turn rate were set as 0.
- Initial variance of position\_x, position\_y, turn rate, and change in turn rate were set as 0.5.
Since we expect the inital value to not vary more than 0.5.
- Acceleration standard deviation was set as 0.54 since the acceleration normally varies between 0.5 and 1.0
for bicycle and we are tracking bicycle.
- Turn change rate standard deviation was set as PI/3.0.

The resulting algorithm achieves following RMSE for position\_x, position\_y, velocity\_x, and velocity\_y,
RMSE\_radar\_lidar = [0.06, 0.08, 0.28, 0.19]
RMSE\_radar = [0.17, 0.14, 0.66, 0.26]
RMSE\_lidar = [0.21, 0.2, 0.46, 0.28]

NIS\_values\_visulization.jpg provides NIS visualization when NIS is calculated for radar, lidar fusion.

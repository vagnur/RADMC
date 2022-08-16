#ifndef RADMC_COMMON_HH
#define RADMC_COMMON_HH

#include <vector>

class common {

private:

    //Light speed in cm / s^2
    double light_speed = 2.9979245800000e10;
    //Light speed in microns / s^2
    double light_speed_microns = 2.99792458e14;

public:

    //Getter for speed of light
    static double get_light_speed();

    //Getter for speed of light in microns
    static double get_light_speed_microns();

    //This method is going to remap the number of points of the old function to a new function
    //It has several options for the remap process
    //Input : int number_of_old_points -> Number of points of the old function
    //        vector<double> old_points -> Domain of the old function
    //        vector<double> old_function -> Values of the old function at each point
    //        int number_of_new_points -> Number of points of the new function
    //        vector<double> new_points -> Domain of the new function
    //        int elow,eup -> These variables are used in order to know which kind of remap we are going to use. Elow is for the lowest boundary, and eup for the highest boundary.
    //                        =0  Do not allow out of bound
    //                        =1  Take boundary value of the old function
    //                        =2  Logarithmic extrapolation
    //                        =3  Out of bound
    //                        =4  Smooth version of out of bound
    //Output : It returns a vector<double> with the values of the new function at each new point
    static std::vector<double> remap_function(int number_of_old_points, std::vector<double> old_points, std::vector<double> old_function, int number_of_new_points, std::vector<double> new_points, int elow, int eup);



};

#endif //RADMC_COMMON_HH

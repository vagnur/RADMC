#ifndef RADMC_COMMON_HH
#define RADMC_COMMON_HH

#include <vector>
#include <iostream>
#include <cmath>

/*
 * This class, like it name suggest, is a common use class. It's not intended to create an object of this class,
 * so every method of this is going to be static, in order to just use it in certain contexts.
 */

class common {

private:

    //Light speed in cm / s^2
    static constexpr double light_speed = 2.9979245800000e10;
    //Light speed in microns / s^2
    static constexpr double light_speed_microns = 2.99792458e14;
    //4 times PI
    static constexpr double four_pi = 12.5663706143591729538505735331;
    //PI/2
    static constexpr double pi_half = 1.57079632679489661923132169164;
    //PI
    static constexpr double pi = 3.14159265358979323846264338328;
    //2 times PI
    static constexpr double two_pi = 6.28318530717958647692528676656;

public:

    //Getter for speed of light
    static double get_light_speed();
    //Getter for speed of light in microns
    static double get_light_speed_microns();
    //Getter for 4 * PI
    static double get_four_times_pi();

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

    //This method is going to interpolate a remapped function. Since we first remap, and then we do the interpolation,
    //  this method is going to give us a more accurate vector.
    //Input : vector old_function -> Old function to be interpolated
    //        vector old_points -> Domain of the old function
    //        int number_of_old_points -> Number of the old points (basically, number of [x,F(x)])
    //        vector new_points -> New point to generate the interpolation
    //        int number_of_new_points -> Number of new points in the vector
    //        vector remapped -> Vector with the remapped values of the old function
    //Output : vector with the interpolated points, that is going to have a size equal to the number of new points.
    static std::vector<double> interpolation_function(std::vector<double> old_function, std::vector<double> old_points, int number_of_old_points, std::vector<double> new_points, int number_of_new_points, std::vector<double> remapped);

    //This method work by searching a value in a set of values and returning the position of the nearest lower value.
    // For example, if one is looking for the value 0.5 and the set of values is {0.1;0.4;0.7;0.95}, then the method is going
    //  to return position 2 (1 since we are in computer science domain), since 0.4 is the nearest lower value of 0.5.
    //Input : vector xx -> Vector with the set of values
    //        int n -> number of points in the vector xx
    //        double x -> the value that we are looking for
    //        int jlo -> this is like the first guess for our function. We always use the value n like the fist guess
    //Output: a int value that represent the position of the nearest lower neighbour value
    static int hunt(const std::vector<double>& xx, int n, double x, int jlo);

    //Method that is going to calculate the blackbody radiation from a body
    //Input : double temperature -> It's the temperature that we are going to use in the formula
    //        double frequency -> Frequency for the planck function
    //Output : double value that represent the radiation at frequency_i and the assigned temperature
    static double black_body_planck_function(double temperature,double frequency);

    //Method that work by dividing a string according to a character (like split in python)
    //Input : string s -> String that we want to split
    //        string del -> substring that we are going to use to do the division of the string s
    //Output : a vector with the substrings founded in the split process
    static std::vector<std::string> tokenize(std::string s, std::string del = " ");

    static void check_unity_vector(double x, double y, double z);

    static double get_pi_half();

    static double get_pi();

    static double get_two_pi();
};

#endif //RADMC_COMMON_HH

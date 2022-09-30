#ifndef RADMC_GRID_HH
#define RADMC_GRID_HH

//Include of STL libraries
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

/*
 * The object of this class its going to store the grid were the simualtion is going to run.
 * It will store information about the grid, and it has methods to found a point inside the grid.
 */

class grid{

private:
    //This string its going to store the type of grid (cartesian, polar) that we are going to use
    std::string type;
    //This vector indicates if a dimension of the grid is present or not.
    //  For example, if present_points[0] = true, means that the x dimension is present in the grid
    std::vector<bool> present_points;
    //Number of points in each dimension.
    //  *Right now we are not using it, since we can know the quantity of points in the vector with size()
    int number_of_points_x;
    int number_of_points_y;
    int number_of_points_z;
    //Grid walls in each dimension
    std::vector<double> x_points;
    std::vector<double> y_points;
    std::vector<double> z_points;
    //Delta of each dimension (used to find a ray position in the grid)
    double diference_x;
    double diference_y;
    double diference_z;

public:

    //Empty constructor
    grid(void);

    //Method to initialize a 3d cartesian regular grid
    //Input : It has no input
    //Output : It has no output, but its goint to store relevant information of the grid in the object
    //TODO : Generalizar el método para grillas regulares
    void initialize_cartesian_regular(void);

    //void initialize_spherical_regular(void);

    //Method to found a ray position in the grid
    //Input : It has no input
    //Output : An integer that represent the position of the ray in the dimension of the grid
    int found_point_cartesian_regular(double);

    //Method that calculate the delta diference in each dimension
    //Input : It has no input
    //Output : It has no output, but its going to store the calculated deltas in the object
    void calculate_points_delta(void);

    //Getter for each vector dimension of the grid
    const std::vector<double>&get_x_points() const;
    const std::vector<double>&get_y_points() const;
    const std::vector<double>&get_z_points() const;

    //Destructor of the class
    ~grid(void);

};

#endif //RADMC_GRID_HH

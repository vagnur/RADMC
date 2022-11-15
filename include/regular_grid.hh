#ifndef RADMC_REGULAR_GRID_HH
#define RADMC_REGULAR_GRID_HH

#include <general_grid.hh>

#include <fstream>
#include <iostream>

class regular_grid: public general_grid {

protected:

    //Delta of each dimension (used to find a ray position in the grid)
    double diference_x;
    double diference_y;
    double diference_z;
    //Number of points in each dimension.
    int number_of_points_x;
    int number_of_points_y;
    int number_of_points_z;
    //Grid walls in each dimension
    std::vector<double> x_points;
    std::vector<double> y_points;
    std::vector<double> z_points;

public:

    void read_grid_file();

    //Getters
    const std::vector<double>&get_x_points() const;
    const std::vector<double>&get_y_points() const;
    const std::vector<double>&get_z_points() const;
    int get_number_of_points_X() const;
    int get_number_of_points_Y() const;
    int get_number_of_points_Z() const;

};


#endif //RADMC_REGULAR_GRID_HH

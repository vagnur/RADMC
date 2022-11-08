#ifndef RADMC_CARTESIAN_REGULAR_GRID_HH
#define RADMC_CARTESIAN_REGULAR_GRID_HH

#include <regular_grid.hh>
#include <cmath>

class cartesian_regular_grid: public regular_grid {

private:



public:

    cartesian_regular_grid();

    //Method that calculate the delta diference in each dimension
    //Input : It has no input
    //Output : It has no output, but it's going to store the calculated deltas in the object
    void calculate_points_delta(void);

    //Method to found a ray position in the grid
    //Input : Vector ray_position -> Vector with the ray positions at X,Y and Z
    //Output : Vector with the grid position of the ray at X,Y and Z
    std::vector<int> found_point(const std::vector<double>& ray_position);

    void calculate_cell_volume();

    ~cartesian_regular_grid();
};


#endif //RADMC_CARTESIAN_REGULAR_GRID_HH

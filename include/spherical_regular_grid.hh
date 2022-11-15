#ifndef RADMC_SPHERICAL_REGULAR_GRID_HH
#define RADMC_SPHERICAL_REGULAR_GRID_HH

#include <regular_grid.hh>
#include <common.hh>

class spherical_regular_grid: public regular_grid {

private:

    std::vector<std::vector<std::vector<double>>> cell_volume;

public:

    spherical_regular_grid();

    void initialize_grid();

    double get_cell_volume(std::vector<int> grid_position) const;

    void adjust_theta();

    void adjust_phi();

    void calculate_cell_volume();

    std::vector<int> found_ray_position_in_grid(const std::vector<double>& ray_position);

    void calculate_points_delta();

    ~spherical_regular_grid();

};


#endif //RADMC_SPHERICAL_REGULAR_GRID_HH

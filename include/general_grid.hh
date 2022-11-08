#ifndef RADMC_GENERAL_GRID_HH
#define RADMC_GENERAL_GRID_HH

#include <vector>

class general_grid {

protected:

    std::vector<bool> present_dimensions;

public:

    virtual void initialize_grid() = 0;
    virtual void calculate_cell_volume() = 0;
    virtual std::vector<int> found_point(const std::vector<double>& ray_position) = 0;
    virtual int get_number_of_points_X() const = 0;
    virtual int get_number_of_points_Y() const = 0;
    virtual int get_number_of_points_Z() const = 0;
    virtual const std::vector<double>&get_x_points() const = 0;
    virtual const std::vector<double>&get_y_points() const = 0;
    virtual const std::vector<double>&get_z_points() const = 0;
    virtual void calculate_points_delta() = 0;
    virtual double get_cell_volume() const= 0;

};


#endif //RADMC_GENERAL_GRID_HH

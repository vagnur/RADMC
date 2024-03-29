#ifndef RADMC_GENERAL_GRID_HH
#define RADMC_GENERAL_GRID_HH

#include <vector>
#include <string>
#include <photon.hh>

class general_grid {

protected:

    std::string type;
    std::vector<bool> present_dimensions;

public:

    virtual void initialize_grid() = 0;
    virtual void read_grid_file() = 0;
    virtual void calculate_cell_volume() = 0;
    virtual std::vector<int> found_ray_position_in_grid(const std::vector<double>& ray_position) = 0;
    virtual void calculate_photon_new_position(photon &photon_i) = 0;
    virtual void calculate_photon_cell_walls(photon &photon_i) = 0;
    virtual int get_number_of_points_X() const = 0;
    virtual int get_number_of_points_Y() const = 0;
    virtual int get_number_of_points_Z() const = 0;
    virtual const std::vector<double>&get_x_points() const = 0;
    virtual const std::vector<double>&get_y_points() const = 0;
    virtual const std::vector<double>&get_z_points() const = 0;
    virtual void calculate_points_delta() = 0;
    virtual double get_cell_volume(std::vector<int> grid_position) const= 0;
    std::string get_grid_type() const;

};


#endif //RADMC_GENERAL_GRID_HH

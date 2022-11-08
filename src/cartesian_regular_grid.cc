#include <cartesian_regular_grid.hh>

cartesian_regular_grid::cartesian_regular_grid() {
    ;
}

std::vector<int> cartesian_regular_grid::found_point(const std::vector<double>& ray_position){
    //This method assumes a pair number of points
    std::vector<int> grid_points(3);
    grid_points[0] = std::floor(ray_position[0] / this -> diference_x) + (this -> number_of_points_x / 2);
    grid_points[1] = std::floor(ray_position[1] / this -> diference_y) + (this -> number_of_points_y / 2);
    grid_points[2] = std::floor(ray_position[2] / this -> diference_z) + (this -> number_of_points_z / 2);
    return grid_points;
}

void cartesian_regular_grid::calculate_points_delta(void){
    this->diference_x = this->x_points[1] - this->x_points[0];
    this->diference_y = this->y_points[1] - this->y_points[0];
    this->diference_z = this->z_points[1] - this->z_points[0];
}

void cartesian_regular_grid::calculate_cell_volume() {
    this -> cell_volume = this -> diference_x * this -> diference_y * this -> diference_z;
}

cartesian_regular_grid::~cartesian_regular_grid() {
    ;
}

#include <cartesian_regular_grid.hh>

cartesian_regular_grid::cartesian_regular_grid() {
    ;
}

void cartesian_regular_grid::initialize_grid() {
    this -> type = "cartesian";
    this -> read_grid_file();
    this -> calculate_points_delta();
    this -> calculate_cell_volume();
}

void cartesian_regular_grid::calculate_photon_cell_walls(photon &photon_i) {
    std::vector<double> cell_walls(3);
    cell_walls[0] = this -> x_points[photon_i.get_grid_position()[0] + photon_i.get_orientation()[0]];
    cell_walls[1] = this -> y_points[photon_i.get_grid_position()[1] + photon_i.get_orientation()[1]];
    cell_walls[2] = this -> z_points[photon_i.get_grid_position()[2] + photon_i.get_orientation()[2]];
    photon_i.set_walls(cell_walls);
}

void cartesian_regular_grid::calculate_photon_new_position(photon &photon_i) {
    std::vector<double> distance(3);
    std::vector<double> ray_position(3);
    std::vector<int> signs = {-1, 1};
    std::vector<int> indexes = {-1, -1, -1};
    std::vector<int> grid_position = photon_i.get_grid_position();
    int count = 0;
    double min_distance;
    this->calculate_photon_cell_walls(photon_i);
    distance[0] = (photon_i.get_cell_walls()[0] - photon_i.get_ray_position()[0]) / photon_i.get_direction()[0];
    distance[1] = (photon_i.get_cell_walls()[1] - photon_i.get_ray_position()[1]) / photon_i.get_direction()[1];
    distance[2] = (photon_i.get_cell_walls()[2] - photon_i.get_ray_position()[2]) / photon_i.get_direction()[2];
    min_distance = std::min(std::min(distance[0], distance[1]), distance[2]);
    for (int i = 0; i < 3; ++i) {
        if (distance[i] == min_distance) {
            indexes[count] = i;
            count++;
        }
    }
    ray_position[0] = photon_i.get_ray_position()[0] + min_distance * photon_i.get_direction()[0];
    ray_position[1] = photon_i.get_ray_position()[1] + min_distance * photon_i.get_direction()[1];
    ray_position[2] = photon_i.get_ray_position()[2] + min_distance * photon_i.get_direction()[2];

    //avoid bug assign cellWall to ray position
    //update grid position with signs
    for (int i = 0; i < count; i++) {
        ray_position[indexes[i]] = photon_i.get_cell_walls()[indexes[i]];
        //TODO : Revisar signo +
        grid_position[indexes[i]] = grid_position[indexes[i]] + signs[photon_i.get_orientation()[indexes[i]]];
    }

    photon_i.set_distance(distance);
    photon_i.set_min_distance(min_distance);
    photon_i.set_ray_position(ray_position);
    photon_i.set_grid_position(grid_position);
    photon_i.is_on_grid(this -> number_of_points_x, this -> number_of_points_y, this -> number_of_points_z);
}

std::vector<int> cartesian_regular_grid::found_ray_position_in_grid(const std::vector<double>& ray_position){
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

double cartesian_regular_grid::get_cell_volume(std::vector<int> grid_position) const {
    return this -> cell_volume;
}

cartesian_regular_grid::~cartesian_regular_grid() {
    ;
}

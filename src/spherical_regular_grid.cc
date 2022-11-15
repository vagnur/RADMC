#include "spherical_regular_grid.hh"

spherical_regular_grid::spherical_regular_grid() {
    ;
}

void spherical_regular_grid::initialize_grid(){
    this -> read_grid_file();
    this -> adjust_theta();
    this -> adjust_phi();
    this -> calculate_points_delta();
    this -> calculate_cell_volume();
}

void spherical_regular_grid::adjust_theta() {

    if((this -> y_points[0] != 0.0 ) && (std::abs(this -> y_points[0]) < 1.0e-4 * std::abs(this -> y_points[1] - this -> y_points[0]))){
        std::cout << "Adjusting theta(0) to exactly 0" << std::endl;
        this -> y_points[0] = 0.0;
    }

    if((this -> y_points[this -> number_of_points_y-1] != common::get_pi_half()) && (std::abs(this -> y_points[this -> number_of_points_y-1] - common::get_pi_half()) < 1.0e-4*std::abs(this -> y_points[this -> number_of_points_y-1] - this -> y_points[this -> number_of_points_y - 2]))){
        std::cout << "Adjusting theta(ny) to exactly pi/2..." << std::endl;
        this -> y_points[this -> number_of_points_y-1] = common::get_pi_half();
    }

    if((this -> y_points[this -> number_of_points_y-1] != common::get_pi()) && (std::abs(this -> y_points[this -> number_of_points_y-1] - common::get_pi()) < 1.0e-4*std::abs(this -> y_points[this -> number_of_points_y-1] - this -> y_points[this -> number_of_points_y - 2]))){
        std::cout << "Adjusting theta(ny) to exactly pi..." << std::endl;
        this -> y_points[this -> number_of_points_y-1] = common::get_pi();
    }

    for (int i = 1; i < this -> number_of_points_y - 1; ++i) {
        if((this -> y_points[i] != common::get_pi_half()) && (std::abs(this -> y_points[i]-common::get_pi_half()) < 0.5e-4*std::abs(this -> y_points[i + 1] - this -> y_points[i - 1]))){
            std::cout << "Adjusting theta(" << i <<") to exactly pi/2..." << std::endl;
        }
    }

}

void spherical_regular_grid::adjust_phi() {

    if((this -> z_points[0] != 0.0) && (std::abs(this -> z_points[0]) < 1.0e-2 )){
        std::cout << "Adjusting phi(0) to exactly 0..." << std::endl;
        this -> z_points[0] = 0.0;
    }

    if((this -> z_points[this -> number_of_points_z - 1] != common::get_two_pi()) && (std::abs(this -> z_points[this -> number_of_points_z - 1] - common::get_two_pi()) < 1.0e-2)){
        std::cout << "Adjusting phi(nz) to exactly 2*PI..." << std::endl;
        this -> z_points[this -> number_of_points_z - 1] = common::get_two_pi();
    }

}

void spherical_regular_grid::calculate_cell_volume(){
    this -> cell_volume.resize(this -> number_of_points_z - 1, std::vector<std::vector<double>>(this -> number_of_points_y - 1, std::vector<double>(this -> number_of_points_x - 1)));
    double volumen;

    for (int k = 0; k < this -> number_of_points_z - 1; ++k) {
        for (int j = 0; j < this -> number_of_points_y - 1; ++j) {
            for (int i = 0; i < this -> number_of_points_x - 1; ++i) {
                volumen = (1.0/3.0) * (std::pow(this -> x_points[i+1],3) - std::pow(this -> x_points[i],3)) * (std::cos(this -> y_points[j]) - std::cos(this -> y_points[j+1])) * (this -> z_points[k+1] - this -> z_points[k]);
                this -> cell_volume[k][j][i] = volumen;
            }
        }
    }
}

std::vector<int> spherical_regular_grid::found_ray_position_in_grid(const std::vector<double>& ray_position){
    std::vector<int> grid_points(3);
    grid_points[0] = std::floor(ray_position[0] / (this -> diference_x + this -> x_points[0]));
    grid_points[1] = std::floor(ray_position[1] / (this -> diference_y + this -> y_points[0]));
    grid_points[2] = std::floor(ray_position[2] / (this -> diference_z + this -> z_points[0]));
    return grid_points;
}

void spherical_regular_grid::calculate_points_delta(){
    this->diference_x = this->x_points[1] - this->x_points[0];
    this->diference_y = this->y_points[1] - this->y_points[0];
    this->diference_z = this->z_points[1] - this->z_points[0];
}

spherical_regular_grid::~spherical_regular_grid() {
    ;
}

double spherical_regular_grid::get_cell_volume(std::vector<int> grid_position) const{
    double volumen = (1.0/3.0) * (std::pow(this -> x_points[grid_position[0]+1],3) - std::pow(this -> x_points[grid_position[0]],3)) * (std::cos(this -> y_points[grid_position[1]]) - std::cos(this -> y_points[grid_position[1]+1])) * (this -> z_points[grid_position[2]+1] - this -> z_points[grid_position[2]]);
    return volumen;
}

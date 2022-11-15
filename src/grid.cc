#include <grid.hh>

grid::grid(void){
    ;
}
void grid::initialize_cartesian_regular(void){
    //We declare the type of grid
    this -> type = "cartesian";
    //Variables to store the metadata from the input file
    int i_format,grid_style,coord_system,grid_info,amount_x,amount_y,amount_z;
    //These variables are used to indicate which coordinates are present in the file
    bool include_x,include_y,include_z;
    //Variable to store the readed value
    double coordinate_value;
    //We read the input file for the grid
    std::ifstream input_file;
    //TODO : Pedirle al usuario una carpeta y abrir todo desde ahi. Dejar en blanco para carpeta inputs
    input_file.open("inputs/amr_grid.inp");
    //In case that the file couldn't be opened, an error message is displayed.
    if(!input_file){
        std::cerr << "Mandatory file \"amr_grid.inp\" could not be opened. Make sure that the file exists" << std::endl;
        exit(0);
    }
    //We read the metadata present in the file
    input_file >> i_format;
    input_file >> grid_style;
    input_file >> coord_system;
    input_file >> grid_info;
    //Next we read the values that indicate which coordinates are present in the file
    input_file >> include_x;
    input_file >> include_y;
    input_file >> include_z;
    this -> present_points.resize(3);
    this -> present_points[0] = include_x;
    this -> present_points[1] = include_y;
    this -> present_points[2] = include_z;
    //Now we read the number of points of each coordinate
    input_file >> amount_x;
    input_file >> amount_y;
    input_file >> amount_z;
    this -> number_of_points_x = amount_x;
    this -> number_of_points_y = amount_y;
    this -> number_of_points_z = amount_z;
    //To finish the process, we read the values of each coordinate that are present in the file
    //We make space for the values
    this -> x_points.resize(amount_x+1);
    this -> y_points.resize(amount_y+1);
    this -> z_points.resize(amount_z+1);
    //First for the X dimension
    if (include_x) {
        for (int i = 0; i <= amount_x; ++i) {
            input_file >> coordinate_value;
            //this->x_points.push_back(coordinate_value);
            this -> x_points[i] = coordinate_value;
        }
    }
    //Then for the Y dimension
    if (include_y) {
        for (int i = 0; i <= amount_y; ++i) {
            input_file >> coordinate_value;
            //this->y_points.push_back(coordinate_value);
            this -> y_points[i] = coordinate_value;
        }
    }
    //At last for the Z dimension
    if (include_z) {
        for (int i = 0; i <= amount_z; ++i) {
            input_file >> coordinate_value;
            //this->z_points.push_back(coordinate_value);
            this -> z_points[i] = coordinate_value;
        }
    }
    //We close the file and the process end
    input_file.close();
    //Then we calculate the diference of each dimension
    this->calculate_points_delta();
}

std::vector<int> grid::found_point_cartesian_regular(const std::vector<double>& ray_position){
    //This method assumes a pair number of points
    std::vector<int> grid_points(3);
    grid_points[0] = std::floor(ray_position[0] / this -> diference_x) + (this -> number_of_points_x / 2);
    grid_points[1] = std::floor(ray_position[1] / this -> diference_y) + (this -> number_of_points_y / 2);
    grid_points[2] = std::floor(ray_position[2] / this -> diference_z) + (this -> number_of_points_z / 2);
    return grid_points;
}

void grid::calculate_points_delta(void){
    this->diference_x = this->x_points[1] - this->x_points[0];
    this->diference_y = this->y_points[1] - this->y_points[0];
    this->diference_z = this->z_points[1] - this->z_points[0];
}

const std::vector<double>& grid::get_x_points() const{
    return this -> x_points;
}

const std::vector<double>& grid::get_y_points() const{
    return this -> y_points;
}

const std::vector<double>& grid::get_z_points() const{
    return this -> z_points;
}

int grid::get_number_of_points_X() const {
    return this -> number_of_points_x;
}

int grid::get_number_of_points_Y() const {
    return this -> number_of_points_y;
}

int grid::get_number_of_points_Z() const {
    return this -> number_of_points_z;
}

double grid::get_cell_volume() const{
    return this -> diference_x*diference_y*diference_z;
}

grid::~grid(void){
    ;
}

/*
void grid::initialize_spherical_regular(void){
    //We declare the type of grid
    this -> type = "spherical";
    //Variables to store the metadata from the input file
    int i_format,grid_style,coord_system,grid_info,amount_x,amount_y,amount_z;
    //This variables are used to indicate which coordinates are present in the file
    bool include_x,include_y,include_z;
    //Variable to store the readed value
    double coordinate_value;
    //We read the input file for the grid
    std::ifstream input_file;
    //TODO : Pedirle al usuario una carpeta y abrir todo desde ahi. Dejar en blanco para carpeta inputs
    input_file.open("inputs/amr_grid.inp");
    //In case that the file couldn't be opened, an error message is displayed.
    if(!input_file){
        std::cerr << "Mandatory file \"amr_grid.inp\" could not be opened. Make sure that the file exists" << std::endl;
        exit(0);
    }
    //We read the metadata present in the file
    input_file >> i_format;
    input_file >> grid_style;
    input_file >> coord_system;
    input_file >> grid_info;
    //Next we read the values that indicate which coordinates are present in the file
    input_file >> include_x;
    input_file >> include_y;
    input_file >> include_z;
    //Now we read the number of points of each coordinate
    input_file >> amount_x;
    input_file >> amount_y;
    input_file >> amount_z;
    //To finish the process, we read the values of each coordinate that are present in the file
    //First for the X dimension
    if (include_x) {
        for (int i = 0; i <= amount_x; ++i) {
            input_file >> coordinate_value;
            this->x_points.push_back(coordinate_value);
        }
    }
    //Then for the Y dimension
    if (include_y) {
        for (int i = 0; i <= amount_y; ++i) {
            input_file >> coordinate_value;
            this->y_points.push_back(coordinate_value);
        }
    }
    //At last for the Z dimension
    if (include_z) {
        for (int i = 0; i <= amount_z; ++i) {
            input_file >> coordinate_value;
            this->z_points.push_back(coordinate_value);
        }
    }
    //We close the file and the process end
    input_file.close();
    this->test_vectors();
    this->calculate_points_delta();
}
 */
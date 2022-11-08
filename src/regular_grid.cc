#include "regular_grid.hh"

void regular_grid::initialize_grid() {
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
    this -> present_dimensions.resize(3);
    this -> present_dimensions[0] = include_x;
    this -> present_dimensions[1] = include_y;
    this -> present_dimensions[2] = include_z;
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
}

const std::vector<double>& regular_grid::get_x_points() const{
    return this -> x_points;
}

const std::vector<double>& regular_grid::get_y_points() const{
    return this -> y_points;
}

const std::vector<double>& regular_grid::get_z_points() const{
    return this -> z_points;
}

int regular_grid::get_number_of_points_X() const {
    return this -> number_of_points_x;
}

int regular_grid::get_number_of_points_Y() const {
    return this -> number_of_points_y;
}

int regular_grid::get_number_of_points_Z() const {
    return this -> number_of_points_z;
}

double regular_grid::get_cell_volume() const {
    return cell_volume;
}
#include <grid.hh>

grid::grid(void){
    ;
}
void grid::initialize_cartesian_regular(void){
    //We declare the type of grid
    this -> type = "cartesian";
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
    //TODO : delete this test
    this->test_vectors();
    this->calculate_points_delta();
}

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

void grid::test_vectors(void){
    for (int i = 0; i < this->x_points.size(); ++i) {
        std::cout << this->x_points[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < this->y_points.size(); ++i) {
        std::cout << this->y_points[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < this->z_points.size(); ++i) {
        std::cout << this->z_points[i] << " ";
    }
    std::cout << std::endl;
}

double grid::found_point(double x){
    if(this->type == "cartesian"){
        int number_of_points = this->x_points.size();
        return std::floor(x/this->diference_x) + (number_of_points/2);
    }
    if(this->type == "spherical"){
        return std::floor(x/this->diference_x);
    }
}

void grid::calculate_points_delta(void){
    this->diference_x = this->x_points[1] - this->x_points[0];
    this->diference_y = this->y_points[1] - this->y_points[0];
    this->diference_z = this->z_points[1] - this->z_points[0];
    std::cout << this->diference_x << std::endl;
}

grid::~grid(void){
    ;
}
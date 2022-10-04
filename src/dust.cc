#include <dust.hh>

dust::dust(void){
    ;
}

void dust::read_dust_species_density(int number_of_point_x,int number_of_point_y,int number_of_point_z){
    //iformat es 1 y no se da mayor info
    //number of cells is the total number of cells in the grid
    //number of species indicates the number of dust species present in the simulation
    int iformat,number_of_cells,number_of_species;
    //Vector that is going to store de density of each dust specie
    std::vector<double> specie_densities;
    //Initializing a 3D vector to store the density of a specie
    //vector<vector<vector<double>>> specie_density(number_of_point_z, vector<vector<double>>(number_of_point_y, vector<double>(number_of_point_x)));
    std::vector<std::vector<std::vector<double>>> specie_density;
    double density;
    //We read the input file for the dust species density
    // TODO : Pedir carpeta, dejar en blanco para inputs
    std::ifstream input_file;
    input_file.open("inputs/dust_density.inp");
    //In case that the file couldn't be opened, an error message is displayed.
    if(!input_file){
        std::cerr << "Mandatory file \"dust_density.inp\" could not be opened. Make sure that the file exists" << std::endl;
        exit(0);
    }
    //Now we read the header information
    input_file >> iformat;
    this -> iformat_dust_density = iformat;
    input_file >> number_of_cells;
    this -> number_of_cells = number_of_cells;
    input_file >> number_of_species;
    this -> number_of_species = number_of_species;
    this -> dust_species_information.resize(number_of_species);
    //Then we read the densities of each specie
    for (int specie = 0; specie < number_of_species; ++specie) {
        //We make space in the 3d vector
        specie_density.resize(number_of_point_z, std::vector<std::vector<double>>(number_of_point_y, std::vector<double>(number_of_point_x)));
        for (int i = 0; i < number_of_point_z; ++i) {
            for (int j = 0; j < number_of_point_y; ++j) {
                for (int k = 0; k < number_of_point_x; ++k) {
                    //And we store the density present in the point [x,y,z]
                    //TODO : Esto funciona para grillas regulares cartesianas 3D. El método no es generalizable creo, pero hay que hacer el método para otras grillas
                    input_file >> density;
                    specie_density[i][j][k] = density;
                }
            }
        }
        //Now we assign the density to the dust specie
        this -> dust_species_information[specie].set_density(specie_density);
        //Then we clear the 3d vector in order to pass to the next specie
        //TODO : ¿Será necesario borrar la memoria de cada dimension?
        //TODO : Investigar y ver efecto en rendimiento
        specie_density.clear();
    }
    //We finish by clossing the input file
    input_file.close();
}

void dust::read_opacities_meta(const std::vector<double>& frequencies){
    //Metadata from the "dustopac.inp" file.
    //Iformat is always 2.
    //number of species is the number of species present in the file
    //TODO : Comprobar que sean las mismas ingresadas anteriormente
    int iformat,number_of_species;
    //Data for each dust specie
    //input style can be 1 or 10
    //iquantum is always 0 for now
    int input_style, iquantum;
    //name of the specie name
    std::string specie_name;
    //First we open the input file
    std::ifstream input_file;
    //TODO : Generalizar la carpeta de entrada
    input_file.open("inputs/dustopac.inp");
    //We read the iformat from the file
    input_file >> iformat;
    //We read the number of species from the file
    input_file >> number_of_species;
    for (int i = 0; i < number_of_species; ++i) {
        //TODO : El archivo de entrada puede tener comentarios y elementos raros, hay que ver como generalizar
        //TODO : el método para que funcione con archivos con comentarios y sin comentarios
        //Then we read the input style
        input_file >> input_style;
        //We read the commentary from the parameter
        //input_file >> specie_name;
        //We read the iquantum
        input_file >> iquantum;
        //We read the commentary from the parameter
        //input_file >> specie_name;
        //Now we read the name of the dust specie, and then we read the relevant file with the method read_opacities
        input_file >> specie_name;
        this -> read_opacities(i,input_style,specie_name, frequencies);
        //Then we read the separation line (only because the format of the input file)
        //NOTE : We are not reading the specie name, but it's not necessary to make a new string
        input_file >> specie_name;
    }
    input_file.close();
}

void dust::read_opacities(int specie_position,int input_style, std::string specie_name, const std::vector<double>& frequencies){
    //TODO : Por ahora asumo que input_style es 1, pero hay que generalizarlo para que pueda ser 10
    //Metadata from the dustkappa_*.inp file
    //iformat indicates the information present in the file
    //number of lambdas indicate the number of wavelength points in the file. Note that is not necessary the same
    //  as the file wavelength_micron.inp (in general it's not the same)
    int iformat,number_of_lambdas;
    //We are going to use this variable to read the data from the file
    double value;
    //We open the input file for the dust specie
    std::ifstream input_file;
    std::string file_name = "inputs/dustkappa_" + specie_name + ".inp";
    input_file.open(file_name);
    //We read the iformat
    input_file >> iformat;
    //We read the number of wavelength points
    input_file >> number_of_lambdas;
    //Vector to store the wavelength points
    std::vector<double> lambda(number_of_lambdas);
    //Vector to store the absorption opacity
    std::vector<double> kappa_absorption(number_of_lambdas);
    //Vector to store (if iformat = 2) the scattering opacity
    std::vector<double> kappa_scattering(number_of_lambdas);
    //Vector to store (if iformat = 3) the mean scattering angle
    std::vector<double> g(number_of_lambdas);
    //We proceed to read the information present in the file
    for (int i = 0; i < number_of_lambdas; ++i) {
        //First we read the wavelength point
        input_file >> value;
        lambda[i] = value;
        //Then we read the absorption opacity
        input_file >> value;
        kappa_absorption[i] = value;
        //If input style its 2, we add the kappa_scattering column and we read its value
        if(iformat == 2 or iformat == 3){
            input_file >> value;
            kappa_scattering[i] = value;
        }
        //If input style its 3, we add the kappa_scattering column and g column and we read its value
        if(iformat == 3){
            input_file >> value;
            g[i] = value;
        }
    }
    //We set each vector to the dust specie object
    this -> dust_species_information[specie_position].set_lambda(lambda);
    this -> dust_species_information[specie_position].set_kappa_absorption(kappa_absorption);
    if(iformat == 2 or iformat == 3){
        this -> dust_species_information[specie_position].set_kappa_scattering(kappa_scattering);
    }
    if(iformat == 3){
        this -> dust_species_information[specie_position].set_g(g);
    }
    //TODO : El original comienza a confirmar datos y propiedades de los datos leídos
    //Then we need to convert the wavelenght points to frequency points
    std::vector<double> frequency = this -> convert_lambda_to_frequency(lambda,number_of_lambdas);
    this -> dust_species_information[specie_position].set_frequency(frequency);
    //TODO : Ahora mismo, tenemos distinta cantidad de puntos en los vectores, por lo que es necesario
    //TODO : remaper los puntos al dominio de las frequencias
    this -> remap_opacities_values(specie_position,frequencies,iformat);
}



void dust::initialize_specie_temperature(int number_of_points_x,int number_of_points_y,int number_of_points_z){
    for (int i = 0; i < this -> number_of_species; ++i) {
        this -> dust_species_information[i].initialize_temperature(number_of_points_x,number_of_points_y,number_of_points_z);
    }
}

//TODO : Este metodo es general ya que se usa en varias partes del codigo, evaluar si dejar en common
std::vector<double> dust::convert_lambda_to_frequency(std::vector<double> lambda, int number_of_lambdas){
    //The vector that it's going to store the calculated frequency points
    std::vector<double> frequency(number_of_lambdas);
    //The point that we are going to calculate
    double frequency_value;
    //The value 2.99792458E+14 comes from the value of speed of light in micrometers
    //double light_speed_microns = common::get_light_speed_microns();
    for (int i = 0; i < number_of_lambdas; ++i) {
        //In order to convert from micrometers to Hz, we need to do: Hz = light speed in micrometers / wavelenght in micrometers
        frequency_value = 2.9979e14 / lambda[i];
        frequency[i] = frequency_value;
    }
    return frequency;
}

const std::vector<dust_species>& dust::get_dust_species(void) const{
    return this -> dust_species_information;
}

int dust::get_number_of_dust_species(void) const{
    return this -> number_of_species;
}

dust::~dust(void){
    ;
}

void dust::remap_opacities_values(int specie, const std::vector<double>& frequencies, int iformat) {
    std::vector<double> specie_frequencies;
    std::vector<double> kappa_absorption;
    std::vector<double> kappa_scattering;
    std::vector<double> g;
    std::vector<double> kappa_absorption_remapped;
    std::vector<double> kappa_scattering_remapped;
    std::vector<double> g_remapped;
    std::vector<double> kappa_absorption_interpoled;
    std::vector<double> kappa_scattering_interpoled;
    std::vector<double> g_interpoled;
    specie_frequencies = this -> dust_species_information[specie].get_frequency();
    //Remap for kappa absorption
    kappa_absorption = this -> dust_species_information[specie].get_absoprtion();
    kappa_absorption_remapped = common::remap_function(kappa_absorption.size(),specie_frequencies,kappa_absorption,frequencies.size(),frequencies,2,1);
    this -> dust_species_information[specie].set_kappa_absorption_remapped(kappa_absorption_remapped);
    kappa_absorption_interpoled = common::interpolation_function(kappa_absorption,specie_frequencies,kappa_absorption.size(),frequencies,frequencies.size(),kappa_absorption_remapped);
    this -> dust_species_information[specie].set_kappa_absorption_interpoled(kappa_absorption_interpoled);
    //Remap for kappa scattering
    if(iformat == 2 or iformat == 3) {
        kappa_scattering = this->dust_species_information[specie].get_scattering();
        kappa_scattering_remapped = common::remap_function(kappa_scattering.size(),specie_frequencies,kappa_scattering,frequencies.size(),frequencies,2,1);
        this -> dust_species_information[specie].set_kappa_scattering_remapped(kappa_scattering_remapped);
        kappa_scattering_interpoled = common::interpolation_function(kappa_scattering,specie_frequencies,kappa_scattering.size(),frequencies,frequencies.size(),kappa_absorption_remapped);
        this -> dust_species_information[specie].set_kappa_scattering_interpoled(kappa_scattering_interpoled);
    }
    if(iformat == 3){
        g = this->dust_species_information[specie].get_g();
        g_remapped = common::remap_function(g.size(),specie_frequencies,g,frequencies.size(),frequencies,1,1);
        this -> dust_species_information[specie].set_g_remapped(g_remapped);
        g_interpoled = common::interpolation_function(g,specie_frequencies,g.size(),frequencies,frequencies.size(),g_remapped);
        this -> dust_species_information[specie].set_g_interpoled(g_interpoled);
    }
}

std::vector<dust_species> &dust::get_dust_species_to_change(void) {
    return this -> dust_species_information;
}

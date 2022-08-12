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
        specie_density.clear();
    }
    //We finish by clossing the input file
    input_file.close();
}

void dust::read_opacities_meta(void){
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
        //First we read the separation line (only because the format of the input file)
        //NOTE : We are not reading the specie name, but it's not necessary to make a new string
        input_file >> specie_name;
        //Then we read the input style
        input_file >> input_style;
        //We read the commentary from the parameter
        input_file >> specie_name;
        //We read the iquantum
        input_file >> iquantum;
        //We read the commentary from the parameter
        input_file >> specie_name;
        //Now we read the name of the dust specie, and then we read the relevant file with the method read_opacities
        input_file >> specie_name;
        read_opacities(i,input_style,specie_name);
    }
    input_file.close();
}

void dust::read_opacities(int specie_position,int input_style, std::string specie_name){
    //TODO : Por ahora asumo que input_style es 1, pero hay que generalizarlo para que pueda ser 10
    //Metadata from the dustkappa_*.inp file
    //iformat indicates the information present in the file
    //number of lambas indicate the number of wavelength points in the file. Note that is not necessary the same
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
        //If input style its 2, we add the kappa_scattering column and we read it's value
        if(input_style == 2){
            input_file >> value;
            kappa_scattering[i] = value;
        }
        //If input style its 3, we add the kappa_scattering column and g column and we read it's value
        if(input_style == 3){
            input_file >> value;
            g[i] = value;
        }
    }
    //We set each vector to the dust specie object
    this -> dust_species_information[specie_position].set_lambda(lambda);
    this -> dust_species_information[specie_position].set_kappa_absorption(kappa_absorption);
    if(iformat == 2){
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
    //TODO : No quiero copiar el código original, porque raro, pero no sé cómo se llaman las técnicas de remapeo :c
}

std::vector<double> dust::convert_lambda_to_frequency(std::vector<double> lambda, int number_of_lambdas){
    //The vector that it going to store the calculated frequency points
    std::vector<double> frequency(number_of_lambdas);
    //The point that we are going to calculate
    double frequency_value;
    //The value 2.99792458E+14 comes from the value of speed of light in micrometers
    double convert_factor = 2.99792458e14;
    for (int i = 0; i < number_of_lambdas; ++i) {
        //In order to convert from micrometers to Hz, we need to do: Hz = 2.99792458E+14 / wavelenght in micrometers
        frequency_value = convert_factor / lambda[i];
        frequency[i] = frequency_value;
    }
    return frequency;
}

std::vector<double> dust::remap_function(int number_of_old_points, std::vector<double> old_points, std::vector<double> old_function, int number_of_new_points, std::vector<double> new_points, int elow, int eup){
    //This vector it's going to store the remapped values
    std::vector<double> new_function(number_of_new_points);
    //In order to iterate over the vectors, we need to know the step, the start and the end of each vector (old and new)
    int old_step,old_start,old_end;
    int new_step,new_start,new_end;
    //This values are used in the logarithmic interpolation
    double x,y;
    //First, we need to know if the old points and the new points are ascending or descending.
    //Since this are monotonic functions, we need to compare only 2 points
    //Ascending case for old points
    if(old_points[1] >= old_points[0]){
        old_step = 1;
        old_start = 0;
        old_end = number_of_old_points - 1;
    }
    //Descending case for old points
    else{
        old_step = - 1;
        old_start = number_of_old_points - 1;
        old_end = 0;
    }
    //Ascending case for new points
    if(new_points[1] >= new_points[0]){
        new_step = 1;
        new_start = 0;
        new_end = number_of_new_points - 1;
    }
    //Descending case for new points
    else{
        new_step = -1;
        new_start = number_of_new_points - 1;
        new_end = 0;
    }
    //Now we iterate over the new values, from the lower to the highest value
    for (int i = new_start; i <= new_end; i = i + new_step) {
        //If the new point is lower than the lowest boundary...
        if(new_points[i] < old_points[old_start]){
            //We have several options
            if(elow == 0){
                //If elow = 0, means that the method doesn't allow a point out of bound
                std::cerr << "Out of bound point in remapping function (lower boundary)" << std::endl;
                exit(0);
            }
            else if(elow == 1){
                //If elow = 1, means that we take boundary value of old function
                new_function[i] = old_function[old_start];
            }
            else if(elow == 2){
                //If elow = 2, means that we are going to do a logarithmic extrapolation
                if(number_of_old_points > 1){
                    if(old_function[old_start+old_step] > 0 && old_function[old_start] > 0){
                        x = old_function[old_start+old_step]/old_function[old_start];
                        y = (std::log(new_points[i])-std::log(old_points[old_start])) / (std::log(old_points[old_start+old_step])-std::log(old_points[old_start]));
                        new_function[i] = old_function[old_start] * std::pow(x,y);
                    }
                    else{
                        new_function[i] = 0.0;
                    }
                }
                else{
                    new_function[i] = old_function[old_start];
                }
            }
            else if(elow == 3){
                //If elow = 3, means that if the value it's out of bound it's going to be 0
                new_fuction[i] = 0.0;
            }
            else if(elow == 4){
                //If elow = 4, we do a smoother version of elow = 3
                if(number_of_new_points > 1){
                    x = (new_points[i])-old_points[old_start]) / (old_points[old_start+old_step]-old_points[old_start]);
                    if(x > -1.0){
                        new_function[i] = (1-std::abs(x))*old_function[old_start];
                    }
                    else{
                        new_function[i] = 0.0;
                    }
                }
                else{
                    new_function[i] = old_function[old_start];
                }
            }
        }
        //If the new point is higher than the highest boundary...
        else if(new_points[i] > old_points[old_end]){
            //Same as before, we have the same options
            if(eup == 0){
                std::cerr << "Out of bound point in remapping function (higher boundary)" << std::endl;
                exit(0);
            }
            else if(eup == 1){
                new_function[i] = old_function[old_end];
            }
            else if(eup == 2){
                if(number_of_old_points > 1){
                    if(old_function[old_end-old_step] > 0.0 && old_function[old_end] > 0.0){
                        x = old_function[old_end-old_step] / old_function[old_end];
                        y = (std::log(new_points[i]) - std::log(old_points[old_end])) / (std::log(old_points[old_end-old_step])-std::log(old_points[old_end]));
                        new_function[i] = old_function[old_end] * std::pow(x,y);
                    }
                    else{
                        new_function[i] = 0.0;
                    }
                }
                else{
                    new_function[i] = old_function[old_end];
                }
            }
            else if(eup == 3){
                new_function[i] = 0.0;
            }
            else if(eup == 4){
                if(number_of_old_points > 1){
                    x = (new_points[i] - old_points[old_end]) / (old_points[old_end-old_step]-old_points[old_end]);
                    if(x > -1.0){
                        new_function[i] = (1-std::abs(x))*old_function[old_end];
                    }
                    else{
                        new_function[i] = 0.0;
                    }
                }
                else{
                    new_function[i] = old_function[old_end];
                }
            }
        }
        else{
            //In any other case, we are in the domain of the old points
            //Here we have 2 options
            if(av == 0){
                //TODO : Acá se llama a la función hunt (hay que implementarla u________________u)
                // x = hunt(<parameters>)
            }
        }

    }
}

dust::~dust(void){
    ;
}

*** Please tell me who you are.
Run
        git config --global user.email "you@example.com"
git config --global user.name "Your Name"
to set your account's default identity.
Omit --global to set the identity only in this repository.
fatal: empty ident name (for <>) not allowed

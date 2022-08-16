#include <iostream>
#include <map>
#include <vector>
#include <fstream>

std::map<std::string,float> read_main_file();
std::vector<std::string> tokenize(std::string s, std::string del = " ");

int main() {
    std::cout << "Hello, World!" << std::endl;
    /*
     * PROBANDO DIMENSIONES DEL VECTOR EN 3D
    std::vector<std::vector<std::vector<int>>> prueba;// = {{{1,2},{3,4}},{{5,6},{7,8}}};
    prueba.resize(2,std::vector<std::vector<int>>(2,std::vector<int>(2,1)));
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                std::cout << prueba[i][j][k];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    prueba.clear();
    prueba.resize(2,std::vector<std::vector<int>>(2,std::vector<int>(2,3)));
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                std::cout << prueba[i][j][k];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    */
}

std::vector<std::string> tokenize(std::string s, std::string del){
    int start = 0;
    int end = s.find(del);
    std::string word;
    std::vector<std::string> words;
    while (end != -1) {
        word = s.substr(start, end - start);
        words.push_back(word);
        //std::cout << s.substr(start, end - start) << std::endl;
        start = end + del.size();
        end = s.find(del, start);
    }
    word = s.substr(start, end - start);
    words.push_back(word);
    //std::cout << s.substr(start, end - start) << std::endl;
    return words;
}

std::map<std::string,float> read_main_file(){
    //First we create a map with the default values
    std::map<std::string,float> simulation_parameters;
    //###PHOTONS###
    //Number of photons packages to be used in the Monte Carlo simulation
    simulation_parameters["nphot"] = 100000.0;
    //Number of photons packages to be used for the scattering Monte Carlo simulation
    simulation_parameters["nphot_scat"] = 10000.0;
    //Number of photon packages for the scattering Monte Carlo simulations, done during spectrum-calculation
    simulation_parameters["nphot_spec"] = 10000.0;
    //The number of photon packages for the Monte Carlo simulations for the mcmono
    simulation_parameters["nphot_mono"] = 100000.0;
    //Indicates if all photon packages will have the same energy
    simulation_parameters["mc_weighted_photons"] = 1.0;
    //Indicates if MC calculate the photon motion inside cells more efficiently
    simulation_parameters["optimized_motion"] = 0.0;
    //A starting value of the random seed for the Monte Carlo simulation
    //###MONTE CARLO SIMULATION###
    //TODO : Creo que esto debería iniciar desde un modulo random o algo asi
    simulation_parameters["iseed"] = -17933201.0;
    //Can make a faster MC simulation at cost of being less accurate
    simulation_parameters["ifast"] = 0.0;
    //This is the fraction by which the energy in each cell may increase before the temperature is recalculated in the Monte Carlo simulation
    simulation_parameters["enthres"] = 0.01;
    //If 0, then all stars are treated as point-sources. If 1, then all stars are treated as finite-size spheres
    simulation_parameters["istar_sphere"] = 0.0;
    //The temperatures are determined in the Monte Carlo method using tabulated pre-computed integrals
    simulation_parameters["ntemp"] = 1000.0;
    //The lowest pre-calculated temperature
    simulation_parameters["temp0"] = 0.01;
    //The highest pre-calculated temperature.
    simulation_parameters["temp1"] = 1.0e5;
    //###DUST###
    //Indicates if dust files must be included
    simulation_parameters["incl_dust"] = 1.0;
    //If set to 0, then the temperatures of all coexisting dust species are always forced to be the same. If 1, then each dust species is thermally independent of the other.
    simulation_parameters["itempdecoup"] = 1.0;
    //###LINES###
    //Indicates if line files must be included
    simulation_parameters["incl_lines"] = 0.0;
    //TODO : Esta parte tiene que ver con los archivos de lineas que aun no usamos
    simulation_parameters["lines_mode"] = 1.0;
    simulation_parameters["lines_maxdoppler"] = 0.3;
    simulation_parameters["lines_partition_ntempint"] = 1000.0;
    simulation_parameters["lines_partition_temp0"] = 0.1;
    simulation_parameters["lines_partition_temp1"] = 1.0e5;
    simulation_parameters["lines_show_pictograms"] = 0.0;
    //###CAMERA####
    //TODO : Esta parte tiene que ver con las imagenes lo que queda para mas adelante
    simulation_parameters["camera_tracemode"] = 1.0;
    simulation_parameters["camera_nrrefine"] = 100.0;
    simulation_parameters["camera_refine_criterion"] = 1.0;
    simulation_parameters["camera_incl_stars"] = 1.0;
    simulation_parameters["camera_starsphere_nrpix"] = 20.0;
    simulation_parameters["camera_spher_cavity_relres"] = 0.05;
    simulation_parameters["camera_localobs_projection"] = 1.0;
    simulation_parameters["camera_min_dangle"] = 0.05;
    simulation_parameters["camera_max_dangle"] = 0.3;
    simulation_parameters["camera_min_dr"] = 0.003;
    simulation_parameters["camera_diagnostics_subpix"] = 0.0;
    simulation_parameters["camera_secondorder"] = 0.0;
    simulation_parameters["camera_interpol_jnu"] = 0.0;
    //TODO : CATEGORIZAR SEGUN CORRESPONDA
    //Determines whether the output of space-dependent data will be in ASCII form (rto_style=1) or binary form (rto_style=3)
    simulation_parameters["rto_style"] = 1.0;
    //TODO : Verificar y entender uso de este parametro
    simulation_parameters["tgas_eq_tdust"] = 0.0;
    //TODO : Entender este parametro, tiene que ver con el scattering pero hay que revisar
    simulation_parameters["scattering_mode_max"] = 1.0;
    //Variable to store each line from the file
    std::string line;
    //We read the main input file
    std::ifstream input_file;
    //TODO : Pedirle al usuario una carpeta y abrir todo desde ahi. Dejar en blanco para carpeta inputs
    input_file.open("inputs/radmc3d.inp");
    std::vector<std::string> values;
    //In case that the file couldn't be opened, an error message is displayed.
    if(!input_file){
        std::cerr << "Mandatory file \"radmc3d.inp\" could not be opened. Make sure that the file exists" << std::endl;
        exit(0);
    }
    //For each line we obtain the information and we store it in the map
    while (std::getline(input_file, line)) {
        values = tokenize(line);
        simulation_parameters[values[0]] = std::stof(values[2]);
    }
    return simulation_parameters;
}
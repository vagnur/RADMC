#include <monte_carlo.hh>

monte_carlo::monte_carlo(void) {
    ;
}

void monte_carlo::do_monte_carlo_therm_regular_cartesian_grid() {
    //TODO : Cómo determinamos el scattering mode???
    //TODO : RADMC original toma las decisiones a partir de los datos de entrada, creo que prefiero consultarlo
    //TODO : como una entrada o en el radmc3d.inp
    int scattering_mode = 2;
    //At first, we read the main radmc3d.inp file. This file contains several information about the simulation
    //  and its parameters
    std::map<std::string,double> simulation_parameters = read_main_file();

    // ************ SETUP OF THE SIMULATION ENVIRONMENT **********************
    //We read the regular cartesian grid from the "amr_grid.inp" file
    this -> grid_object.initialize_cartesian_regular();
    //We read the frequencies from the "wavelength_micron.inp" file, and we calculate
    //  the mean intesity.
    //Note that freq_nu = frequencies_object.get_frequencies()
    //  and freq_dnu = frequencies_object.get_mean_intensity()
    this -> frequencies_object.read_frequencies();
    this -> frequencies_object.calculate_mean_intensity();
    //We read the star sources from the "stars.inp" and we calculate the spectrum and its cumulative function,
    //  the luminosities, the energy, and we invoke the function "jitter_stars" in order to avoid the star being
    //  exactly in the position of a wall of the grid
    this -> stars_object.read_stars();
    this -> stars_object.calculate_spectrum(this -> frequencies_object.get_frequencies());
    this -> stars_object.calculate_cumulative_spectrum(this -> frequencies_object.get_mean_intensity());
    this -> stars_object.calculate_total_luminosities(this -> frequencies_object.get_mean_intensity());
    this -> stars_object.calculate_energy(simulation_parameters["nphot"]);
    this -> stars_object.fix_luminosities();
    this -> stars_object.jitter_stars(this -> grid_object.get_x_points(),this -> grid_object.get_y_points(),this -> grid_object.get_z_points());
    //We read the dust information.
    //First, we read the "dust_density.inp" file, to obtain the density of the species in the grid
    this -> dust_object.read_dust_species_density(this -> grid_object.get_number_of_points_X(),this -> grid_object.get_number_of_points_Y(),this -> grid_object.get_number_of_points_Z());
    //Then, we read the opacities meta file.
    //  This function also read each dustkappa_* file for each specie name.
    //  It's going to remap and interpolate the readed values according to the frequencies domain
    this -> dust_object.read_opacities_meta(this -> frequencies_object.get_frequencies());
    //The last process for the dust is to initialize the temperatures of each specie in the grid
    this -> dust_object.initialize_specie_temperature(this -> grid_object.get_number_of_points_X(),this -> grid_object.get_number_of_points_Y(),this -> grid_object.get_number_of_points_Z());
    //We calculate the temperatures DB. These values are precalculated temperatures that we are going to use in the simulation.
    this -> emissivity_object.generate_emissivity_table(simulation_parameters,this -> dust_object.get_number_of_dust_species(),this -> dust_object.get_dust_species(),this -> frequencies_object.get_frequencies(),this->frequencies_object.get_mean_intensity());
    this -> emissivity_object.compute_derivate(this -> dust_object.get_number_of_dust_species(),simulation_parameters["ntemp"],this -> frequencies_object.get_mean_intensity());
    //With this, we end the setup process, and we proceed to run the photons
    // *********************** SIMULATION LOGIC *************************
    //First, we initialize the random number generators
    // Will be used to obtain a seed for the random number engine
    std::random_device rd;
    // Standard mersenne_twister_engine seeded with rd()
    std::mt19937 generator(rd());
    //We create a uniform distribution between 0 and 1. With this, we are going to
    //generate random number in that distribution in order to obtain
    //the star source, the photon frequency, the tau path and the scattering status of each photon
    std::uniform_real_distribution<> uniform_zero_one_distribution(0.0, 1.0);

    //Now we set up the photons
    //TODO : Cambiar la logica del programa
    //TODO : El foton sólo debería tener atributos, no métodos (o al menos no todos los métodos considerados hasta ahora)
    this -> photons.resize(simulation_parameters["nphot"]);
    for (int i = 0; i < this -> photons.size(); ++i) {
    //for (int i = 0; i < 10; ++i) {
        std::cout << "Launching photon: " << i+1 << std::endl;
        //TODO : New lleva el objeto al heap. Acorde a lo leído, el heap tiene más espacio, pero es necesario borrar
        this -> photons[i] = photon(generator,
                                    uniform_zero_one_distribution,
                                    this -> dust_object.get_number_of_dust_species(),
                                    this -> stars_object.get_number_of_stars(),
                                    this -> frequencies_object.get_number_frequency_points(),
                                    this -> stars_object.get_cumulative_luminosity(),
                                    this -> stars_object.get_stars_information());
        //TODO : calcular la posición del foton respecto a la pos del rasho
        this -> photons[i].setGridPosition(this -> grid_object.found_point_cartesian_regular(this -> photons[i].getRayPosition()));
        //TODO : Acá inicia el ciclo del camino completo, es saltar al sgte evento y verificar si estamos dentro de la grilla. Si está, hacer los eventos relacionados
        std::cout << "walking to the next event" << std::endl;
        this -> photons[i].walk_next_event(generator,
                                           uniform_zero_one_distribution,
                                           this -> dust_object.get_number_of_dust_species(),
                                           this -> dust_object.get_dust_species_to_change(),
                                           this -> stars_object.get_stars_information(),
                                           this -> grid_object.get_number_of_points_X(),
                                           this -> grid_object.get_number_of_points_Y(),
                                           this -> grid_object.get_number_of_points_Z(),
                                           this -> grid_object.get_x_points(),
                                           this -> grid_object.get_y_points(),
                                           this -> grid_object.get_z_points());
        //While the photon is on the grid
        while(this -> photons[i].get_on_grid_condition()){
            std::cout << "The photon " << i + 1 << " is on the grid" << std::endl;
            //We need to do a scattering event...
            if(this -> photons[i].get_is_scattering_condition()){
                std::cout << "The photon " << i + 1 << " is going to scatter" << std::endl;
                //If we do a scattering event, we need to know the type of scattering
                if(scattering_mode == 1){
                    this -> photons[i].get_random_direction(generator,uniform_zero_one_distribution);
                }
                if(scattering_mode == 2){
                    int scattering_specie = this -> photons[i].find_specie_to_scattering(generator,uniform_zero_one_distribution,this -> dust_object.get_number_of_dust_species());
                    this -> photons[i].get_henvey_greenstein_direction(generator,uniform_zero_one_distribution,this -> dust_object.get_dust_species()[scattering_specie].get_g_interpoled());
                }
            }
            //Or an absorption event
            else{
                std::cout << "The photon " << i + 1 << " is going to do absorption" << std::endl;
                this -> photons[i].do_absorption_event(generator,
                                                       uniform_zero_one_distribution,
                                                       this -> dust_object.get_number_of_dust_species(),
                                                       this -> dust_object.get_dust_species_to_change(),
                                                       this -> stars_object.get_stars_information(),
                                                       this -> emissivity_object.get_db_temp(),
                                                       this -> emissivity_object.get_db_logenertemp(),
                                                       this -> emissivity_object.get_db_enertemp(),
                                                       simulation_parameters["ntemp"],
                                                       this -> frequencies_object.get_number_frequency_points(),
                                                       this -> grid_object.get_cell_volume(),
                                                       this -> emissivity_object.get_db_cumulnorm());
                this -> photons[i].get_random_direction(generator,uniform_zero_one_distribution);
            }
            std::cout << "The photon " << i + 1 << " is generating new tau path" << std::endl;
            this -> photons[i].get_tau_path(generator,uniform_zero_one_distribution);
            std::cout << "The photon " << i + 1 << " is walking to the next event" << std::endl;
            this -> photons[i].walk_next_event(generator,
                                               uniform_zero_one_distribution,
                                               this -> dust_object.get_number_of_dust_species(),
                                               this -> dust_object.get_dust_species_to_change(),
                                               this -> stars_object.get_stars_information(),
                                               this -> grid_object.get_number_of_points_X(),
                                               this -> grid_object.get_number_of_points_Y(),
                                               this -> grid_object.get_number_of_points_Z(),
                                               this -> grid_object.get_x_points(),
                                               this -> grid_object.get_y_points(),
                                               this -> grid_object.get_z_points());
        }
    }

    

    /*
    void photon::walk_full_path(){
        this -> walk_next_event();
        while (this -> on_grid){
            if (this -> is_scattering){
                if(scattering_mode == 1){
                    this -> get_random_direction();
                }
                else if (scattering_mode == 2){
                    int specie = this -> find_specie_to_scattering();
                    this -> getHenveyGreensteinDirection();
                }
            }
            else{
                this -> do_absoprtion_event();
                this -> get_random_direction();
            }
            this -> get_tau_path();
            this -> walk_next_event();
        }
    }
*/

}

std::map<std::string,double> monte_carlo::read_main_file(void){
    //First we create a map with the default values
    std::map<std::string,double> simulation_parameters;
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
        values = common::tokenize(line);
        simulation_parameters[values[0]] = std::stof(values[2]);
    }
    return simulation_parameters;
}

monte_carlo::~monte_carlo(void) {
    ;
}

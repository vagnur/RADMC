#include <monte_carlo.hh>

monte_carlo::monte_carlo(void) {
    ;
}

void monte_carlo::do_monte_carlo_therm_regular_cartesian_grid() {
    //TODO : Cómo determinamos el scattering mode???
    //TODO : RADMC original toma las decisiones a partir de los datos de entrada, creo que prefiero consultarlo
    //TODO : como una entrada o en el radmc3d.inp
    int scattering_mode = 1;
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
    for (unsigned int i = 0; i < this -> photons.size(); ++i) {
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
}

void monte_carlo::do_monte_carlo_therm_regular_cartesian_grid_2() {
    //TODO : Cómo determinamos el scattering mode???
    //TODO : RADMC original toma las decisiones a partir de los datos de entrada, creo que prefiero consultarlo
    //TODO : como una entrada o en el radmc3d.inp
    int scattering_mode = 1;
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

    //We initialize the photons of the simulation
    this -> intialize_cartesian_regular_photons(generator,uniform_zero_one_distribution,simulation_parameters["nphot"]);
    //We launch each photon
    this -> launch_photons(generator,uniform_zero_one_distribution,scattering_mode,simulation_parameters["nphot"],simulation_parameters["ntemp"]);
    //When the photons end, we calculate the temperature of each dust specie from the energy in each cell of the grid
    this -> calculate_dust_temperature();
    //To end the process, we write the temperatures output file
    this -> write_temperatures_file();

    /*
    for (int i = 0; i < simulation_parameters["nphot"]; ++i) {
    //for (int i = 0; i < 10; ++i) {
        if(i%1000 == 0) {
            std::cout << "Launching photon: " << i + 1 << std::endl;
        }
        this -> walk_next_event(this -> photons[i], generator, uniform_zero_one_distribution);
        while(this -> photons[i].get_on_grid_condition()){
            if(this -> photons[i].get_is_scattering_condition()) {
                if(scattering_mode == 1){
                    this->photons[i].get_random_direction(generator, uniform_zero_one_distribution);
                }
                else if(scattering_mode == 2){
                    int scattering_specie = this -> photons[i].find_specie_to_scattering(generator,uniform_zero_one_distribution,this -> dust_object.get_number_of_dust_species());
                    this -> photons[i].get_henvey_greenstein_direction(generator,uniform_zero_one_distribution,this -> dust_object.get_dust_species()[scattering_specie].get_g_interpoled());
                }
            }
            else{
                this -> do_absorption_event(generator,uniform_zero_one_distribution,this -> photons[i],simulation_parameters["ntemp"]);
                this -> photons[i].get_random_direction(generator,uniform_zero_one_distribution);
            }
            this -> photons[i].get_tau_path(generator,uniform_zero_one_distribution);
            this -> walk_next_event(this -> photons[i], generator, uniform_zero_one_distribution);
        }
    }
     */

    /*
    double cumuled_energy, temp;
    //When the photons end, we calculate the temperature of each dust specie from the energy in each cell of the grid
    for (int iSpecie = 0; iSpecie < dust_object.get_number_of_dust_species(); ++iSpecie) {
        for (int i = 0; i < grid_object.get_number_of_points_X(); ++i) {
            for (int j = 0; j < grid_object.get_number_of_points_Y(); ++j) {
                for (int k = 0; k < grid_object.get_number_of_points_Z(); ++k) {
                    cumuled_energy = dust_object.get_dust_species()[iSpecie].get_cumulative_energy()[k][j][i] / (dust_object.get_dust_species()[iSpecie].get_densities()[k][j][i] * grid_object.get_cell_volume());
                    if (cumuled_energy <= 0){
                        dust_object.set_null_temperature(iSpecie,i,j,k);
                    }
                    else{
                        temp = emissivity_object.compute_dust_temp_energy(cumuled_energy,iSpecie);
                        dust_object.set_specie_temperature_at_position(iSpecie,i,j,k,temp);
                    }
                }
            }
        }
    }
    */

}

void monte_carlo::write_temperatures_file(){
    std::ofstream temperatures_file("dust_temperature.dat");
    for (int iSpecie = 0; iSpecie < dust_object.get_number_of_dust_species(); ++iSpecie) {
        for (int i = 0; i < grid_object.get_number_of_points_X(); ++i) {
            for (int j = 0; j < grid_object.get_number_of_points_Y(); ++j) {
                for (int k = 0; k < grid_object.get_number_of_points_Z(); ++k) {
                    temperatures_file << dust_object.get_dust_species()[iSpecie].get_temperature()[k][j][i] << "\n";
                }
            }
        }
    }
    temperatures_file.close();
}

void monte_carlo::calculate_dust_temperature(){
    double cumuled_energy, temp;
    //When the photons end, we calculate the temperature of each dust specie from the energy in each cell of the grid
    for (int iSpecie = 0; iSpecie < this->dust_object.get_number_of_dust_species(); ++iSpecie) {
        for (int i = 0; i < this->grid_object.get_number_of_points_X(); ++i) {
            for (int j = 0; j < this->grid_object.get_number_of_points_Y(); ++j) {
                for (int k = 0; k < this->grid_object.get_number_of_points_Z(); ++k) {
                    cumuled_energy = this->dust_object.get_dust_species()[iSpecie].get_cumulative_energy()[k][j][i] / (this->dust_object.get_dust_species()[iSpecie].get_densities()[k][j][i] * this->grid_object.get_cell_volume());
                    if (cumuled_energy <= 0){
                        dust_object.set_null_temperature(iSpecie,i,j,k);
                    }
                    else{
                        temp = this->emissivity_object.compute_dust_temp_energy(cumuled_energy,iSpecie);
                        this->dust_object.set_specie_temperature_at_position(iSpecie,i,j,k,temp);
                    }
                }
            }
        }
    }
}

void monte_carlo::launch_photons(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,int scattering_mode, int number_of_photons, int number_of_temperatures){
    //We launch each photon
    for (int i = 0; i < number_of_photons; ++i) {
        //for (int i = 0; i < 10; ++i) {
        if(i%1000 == 0) {
            std::cout << "Launching photon: " << i + 1 << std::endl;
        }
        this -> walk_next_event(this -> photons[i], generator, uniform_zero_one_distribution);
        while(this -> photons[i].get_on_grid_condition()){
            if(this -> photons[i].get_is_scattering_condition()) {
                if(scattering_mode == 1){
                    this->photons[i].get_random_direction(generator, uniform_zero_one_distribution);
                }
                else if(scattering_mode == 2){
                    int scattering_specie = this -> photons[i].find_specie_to_scattering(generator,uniform_zero_one_distribution,this -> dust_object.get_number_of_dust_species());
                    this -> photons[i].get_henvey_greenstein_direction(generator,uniform_zero_one_distribution,this -> dust_object.get_dust_species()[scattering_specie].get_g_interpoled());
                }
            }
            else{
                this -> do_absorption_event(generator,uniform_zero_one_distribution,this -> photons[i],number_of_temperatures);
                this -> photons[i].get_random_direction(generator,uniform_zero_one_distribution);
            }
            this -> photons[i].get_tau_path(generator,uniform_zero_one_distribution);
            this -> walk_next_event(this -> photons[i], generator, uniform_zero_one_distribution);
        }
    }
}

void monte_carlo::do_absorption_event(std::mt19937& generator, std::uniform_real_distribution<>& uniform_real_distribution, photon& photon_i, int number_of_temperatures){
    int ix,iy,iz,ray_inu;
    std::vector<double> temp_local(dust_object.get_number_of_dust_species());
    ix = photon_i.get_grid_position()[0];
    iy = photon_i.get_grid_position()[1];
    iz = photon_i.get_grid_position()[2];
    this -> divide_absorved_energy(photon_i);
    this -> add_temperature_decoupled(photon_i);
    for (int i = 0; i < dust_object.get_number_of_dust_species(); ++i) {
        //this -> tempLocal[i] = dust_species_information[i].get_temperature()[iz][iy][ix];
        temp_local[i] = dust_object.get_dust_species()[i].get_temperature()[iz][iy][ix];
        //this->tempLocal[i] = temperatures[i][iz][iy][ix];
    }
    photon_i.set_temp_local(temp_local);
    ray_inu = this -> pickRandomFreqDb(generator,uniform_real_distribution,photon_i,number_of_temperatures);
    photon_i.set_ray_inu(ray_inu);
}

int monte_carlo::pickRandomFreqDb(std::mt19937& generator, std::uniform_real_distribution<>& uniform_real_distribution, photon& photon_i, int number_of_temperatures){
    //float sumTemp1, sumTemp2;
    std::vector<double> enerCum(dust_object.get_number_of_dust_species() + 1);
    int numCumul = frequencies_object.get_number_frequency_points() + 1;
    int intplt = 1;
    double rn;
    int iSpec = 0;
    int iTemp = 0;
    int inuPick = 0;
    //this->enerCum[0] = 0;
    enerCum[0] = 0;
    for (int i = 1; i < (dust_object.get_number_of_dust_species() + 1); ++i) {
        enerCum[i] = enerCum[i - 1] + photon_i.get_ener_part()[i - 1];
    }
    if (dust_object.get_number_of_dust_species() > 1) {
        rn = uniform_real_distribution(generator);
        iSpec = common::hunt(photon_i.get_ener_cum(), dust_object.get_number_of_dust_species() + 1, rn, dust_object.get_number_of_dust_species());
        if ((iSpec < 0 || (iSpec > dust_object.get_number_of_dust_species() - 1))) {
            std::cerr << "ERROR : Specie found out of range ..." << std::endl;
            //exit(0);
            //printf("exit\n");
            //exit(1);
        }
    }
    iTemp = common::hunt(emissivity_object.get_db_temp(), number_of_temperatures, photon_i.get_temp_local()[iSpec], number_of_temperatures);
    double eps = 0.0;

    if (iTemp >= number_of_temperatures - 1) {
        //printf("ERROR: Too high temperature discovered\n");
        std::cerr << "ERROR : Too high temperature discovered" << std::endl;
        //exit(0);
    }
    if (iTemp <= -1) {
        iTemp = 0;
    } else {
        eps = (photon_i.get_temp_local()[iSpec] - emissivity_object.get_db_temp()[iTemp]) / (emissivity_object.get_db_temp()[iTemp + 1] - emissivity_object.get_db_temp()[iTemp]);
        if ((eps > 1) || (eps < 0)) {
            std::cerr << "ERROR : In picking new random frequency, eps out of range" << std::endl;
            //exit(0);
            //printf("ERROR: in pick randomFreq, eps<0 or eps>1\n");
            //printf("exit\n");
            //exit(1);
        }
    }
    std::vector<double> dbCumul(frequencies_object.get_number_frequency_points()+1);
    if (intplt == 1) {
        for (int inu = 0; inu < numCumul; ++inu) {
            dbCumul[inu] = (1.0 - eps) * emissivity_object.get_db_cumulnorm()[iSpec][iTemp][inu] + eps * emissivity_object.get_db_cumulnorm()[iSpec][iTemp + 1][inu];
            //this->dbCumul[inu] =
            //        (1.0 - eps) * dbCumulNorm[iSpec][iTemp][inu] + eps * dbCumulNorm[iSpec][iTemp + 1][inu];
        }
        rn = uniform_real_distribution(generator);
        //rn = 0.20160651049169737;
        inuPick = common::hunt(dbCumul, frequencies_object.get_number_frequency_points(), (double) rn, frequencies_object.get_number_frequency_points());
        //free(dbCumul);
    }
    else {
        //Now, within this species/size, find the frequency
        if (eps > 0.5) {
            if (iTemp < number_of_temperatures - 2) {
                iTemp++;
            }
        }
        rn = uniform_real_distribution(generator);
        //verify dbCumulNorm
        inuPick = common::hunt(emissivity_object.get_db_cumulnorm()[iSpec][iTemp - 1], frequencies_object.get_number_frequency_points(),
                               rn,frequencies_object.get_number_frequency_points());
    }
    return inuPick;
}

void monte_carlo::add_temperature_decoupled(photon& photon_i){
    int ix,iy,iz;
    ix = photon_i.get_grid_position()[0];
    iy = photon_i.get_grid_position()[1];
    iz = photon_i.get_grid_position()[2];
    double cumen, temperature;
    for (int iSpec = 0; iSpec < dust_object.get_number_of_dust_species(); ++iSpec) {
        //TODO : Obtener cumulEner de la especie ispec y densidad en ispec y cell
        //cumen = cumulEner[iSpec][iz][iy][ix] / (densities[iSpec][iz][iy][ix] * cellVolumes);
        cumen = dust_object.get_dust_species()[iSpec].get_cumulative_energy()[iz][iy][ix] / (dust_object.get_dust_species()[iSpec].get_densities()[iz][iy][ix] * grid_object.get_cell_volume());

        //cumen = dust_species_information[iSpec].get_cumulative_energy()[iz][iy][ix] / (dust_species_information[iSpec].get_densities()[iz][iy][ix]*cellVolumes);
        //TODO : Este vector inicia en 0 para cada especie en el objeto de polvo...
        //TODO : Hay que cachar si puedo cargar en cada protón y luego sumar todo o depende de las temps anteriores
        //const std::vector<double>& dbTemp, std::vector<std::vector<double>>& dbLogEnerTemp, const std::vector<std::vector<double>>& dbEnerTemp, int number_of_temperatures, double energy, int iSpec
        //temperature = this ->computeDusttempEnergyBd(dbTemp, dbLogEnerTemp, dbEnerTemp, number_of_temperatures, cumen,iSpec);
        temperature = emissivity_object.compute_dust_temp_energy(cumen, iSpec);
        //dust_species_information[iSpec].set_temperature_at_position(ix,iy,iz,temperature);
        dust_object.set_specie_temperature_at_position(iSpec,ix,iy,iz,temperature);
        //temperatures[iSpec][iz][iy][ix] = this->computeDusttempEnergyBd(dbTemp, dbLogEnerTemp, dbEnerTemp,
        //                                                               number_of_temperatures, cumen, iSpec);
        //dustTemperature->temperatures[iSpec][iz][iy][ix] = computeDusttempEnergyBd(emissivityDb, cumen, iSpec);
    }
}

void monte_carlo::divide_absorved_energy(photon& photon_i){
    int ix,iy,iz;
    double alpha_A = 0;
    std::vector<double> ener_part(dust_object.get_number_of_dust_species());
    ix = photon_i.get_grid_position()[0];
    iy = photon_i.get_grid_position()[1];
    iz = photon_i.get_grid_position()[2];
    if (dust_object.get_number_of_dust_species() == 1) {
        //TODO : Obtener energia de la estrella
        ener_part[0] = stars_object.get_stars_information()[photon_i.get_star_source()].get_energy();
        //this -> enerPart[0] = star_information[this -> star_source].get_energy();

        //this->enerPart[0] = star_energies[this->star_source];
    }
    else {
        for (int i = 0; i < dust_object.get_number_of_dust_species(); ++i) {
            //TODO : Obtener densidad de la especie i
            ener_part[i] = dust_object.get_dust_species()[i].get_densities()[iz][iy][ix]  * dust_object.get_dust_species()[i].get_kappa_absorption_interpoled()[photon_i.get_ray_inu()];
            //this -> enerPart[i] = dust_specie_information[i].get_densities()[iz][iy][ix] * dust_specie_information[i].get_kappa_absorption_interpoled()[this -> ray_inu];
            //this->enerPart[i] = density[i][iz][iy][ix] * kappa_A[i][this->ray_inu];
            //photon->enerPart[i] = dustDensity->densities[i][iz][iy][ix] * dustOpacity->kappaA[i][photon->iFrequency];
            alpha_A += ener_part[i];

        }
        for (int i = 0; i < dust_object.get_number_of_dust_species(); ++i) {
            ener_part[i] = stars_object.get_stars_information()[photon_i.get_star_source()].get_energy() * ener_part[i] / alpha_A;
            //this -> enerPart[i] = star_information[this -> star_source].get_energy() * this -> enerPart[i] / alphaA;
            //this->enerPart[i] = star_energies[this->star_source] * this->enerPart[i] / alphaA;
        }
    }
    photon_i.set_ener_part(ener_part);
}

void monte_carlo::walk_next_event(photon& photon_i,std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution){
    double minor_distance, fraction, add_tmp, dum;
    bool carry_on = true;
    while (carry_on) {
        //First we obtain the actual ray and grid position of the photon
        photon_i.set_prev_ray_position();
        photon_i.set_prev_grid_position();

        //Then we calculate the minor distance to move to the next cell of the grid
        minor_distance = photon_i.advance_next_position(this -> grid_object.get_number_of_points_X(),
                                                        this -> grid_object.get_number_of_points_Y(),
                                                        this -> grid_object.get_number_of_points_Z(),
                                                        this -> grid_object.get_x_points(),
                                                        this -> grid_object.get_y_points(),
                                                        this -> grid_object.get_z_points());
        photon_i.calculate_opacity_coefficients(minor_distance,
                                                      this -> dust_object.get_number_of_dust_species(),
                                                      this -> dust_object.get_dust_species());
        if (photon_i.get_tau_path_gone() + photon_i.get_dtau() > photon_i.get_tau_path_total()) {
            fraction = (photon_i.get_tau_path_total() - photon_i.get_tau_path_gone()) / photon_i.get_dtau();
            photon_i.update_ray_position(fraction);

            //TODO : Acá está la energia
            //double dum = (1.0 - this->albedo) * (this->tau_path_total - this->tau_path_gone) *
            //             stars_information[this->star_source].get_energy() / this->alpha_A_total;
            dum = (1.0 - photon_i.get_albedo()) * (photon_i.get_tau_path_total() - photon_i.get_tau_path_gone()) *
                    stars_object.get_stars_information()[photon_i.get_star_source()].get_energy() / photon_i.get_alpha_A_total();
            //for (int i = 0; i < number_of_species; ++i) {
            for (int i = 0; i < dust_object.get_number_of_dust_species(); ++i) {
                //add_tmp = dum * this->alpha_A_specie[i];
                add_tmp = dum * photon_i.get_alpha_A_specie()[i];
                //TODO : como recibimos y cambiamos el vector?
                //dust_specie_information[i].add_energy(this -> grid_position[0],this -> grid_position[1],this -> grid_position[2],add_tmp);
                dust_object.add_energy_specie(i,photon_i.get_grid_position(),add_tmp);
                //cumulative_energy_specie[i][this->grid_position[2]][this->grid_position[1]][this->grid_position[0]] += add_tmp;
            }
            carry_on = false;
        } else {
            //dum = (1.0 - this->albedo) * this->dtau * stars_information[this->star_source].get_energy() / this->alpha_A_total;
            dum = (1.0 - photon_i.get_albedo()) * photon_i.get_dtau() * stars_object.get_stars_information()[photon_i.get_star_source()].get_energy() / photon_i.get_alpha_A_total();
            //for (int i = 0; i < number_of_species; ++i) {
            for (int i = 0; i < dust_object.get_number_of_dust_species(); ++i) {
                //add_tmp = dum * this->alpha_A_specie[i];
                add_tmp = dum * photon_i.get_alpha_A_specie()[i];
                //dust_specie_information[i].add_energy(this -> prev_grid_position[0],this -> prev_grid_position[1],this -> prev_grid_position[2],add_tmp);
                dust_object.add_energy_specie(i,photon_i.get_prev_grid_position(),add_tmp);
                //cumulative_energy_specie[i][this->grid_position[2]][this->grid_position[1]][this->grid_position[0]] += add_tmp;
            }
            //this->tau_path_gone = this->tau_path_gone + this->dtau;
            photon_i.update_tau_path_gone();
            //carry_on = this->on_grid;
            carry_on = photon_i.get_on_grid_condition();
        }
    }
    double rn = uniform_zero_one_distribution(generator);
    //this->is_scattering = rn < this->albedo;
    photon_i.set_scattering_state(rn);
}

void monte_carlo::intialize_cartesian_regular_photons(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution, int number_of_photons){
    this -> photons.resize(number_of_photons);
    int star, ray_frequency;
    for (unsigned int i = 0; i < this -> photons.size(); ++i) {
        star = stars_object.identify_star(generator,uniform_zero_one_distribution);
        ray_frequency = frequencies_object.get_random_frequency(generator,uniform_zero_one_distribution,stars_object.get_stars_information()[star].get_cumulative_spectrum());
        this -> photons[i] = photon(generator,
                                    uniform_zero_one_distribution,
                                    dust_object.get_number_of_dust_species(),
                                    frequencies_object.get_number_frequency_points(),
                                    star,
                                    ray_frequency,
                                    stars_object.get_stars_information()[star].get_star_position());
        this -> photons[i].setGridPosition(this -> grid_object.found_point_cartesian_regular(this -> photons[i].getRayPosition()));
    }
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

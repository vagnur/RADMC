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
    std::map<std::string, double> simulation_parameters = read_main_file();

    // ************ SETUP OF THE SIMULATION ENVIRONMENT **********************

    //We read the regular cartesian grid from the "amr_grid.inp" file
    this->m_grid.initialize_cartesian_regular();

    //We read the frequencies from the "wavelength_micron.inp" file, and we calculate
    //  the mean intesity.
    //Note that freq_nu = m_frequencies.get_frequencies()
    //  and freq_dnu = m_frequencies.get_mean_intensity()
    this->m_frequencies.read_frequencies();
    this->m_frequencies.calculate_mean_intensity();

    //We read the star sources from the "stars.inp" and we calculate the spectrum and its cumulative function,
    //  the luminosities, the energy, and we invoke the function "jitter_stars" in order to avoid the star being
    //  exactly in the position of a wall of the grid
    this->m_stars.read_stars();
    this->m_stars.calculate_spectrum(this->m_frequencies.get_frequencies());
    this->m_stars.calculate_cumulative_spectrum(this->m_frequencies.get_mean_intensity());
    this->m_stars.calculate_total_luminosities(this->m_frequencies.get_mean_intensity());
    this->m_stars.calculate_energy(simulation_parameters["nphot"]);
    this->m_stars.fix_luminosities();
    this->m_stars.jitter_stars(this->m_grid.get_x_points(), this->m_grid.get_y_points(),
                               this->m_grid.get_z_points());

    //We read the dust information.
    //First, we read the "dust_density.inp" file, to obtain the density of the species in the grid
    this->m_dust.read_dust_species_density(this->m_grid.get_number_of_points_X(),
                                           this->m_grid.get_number_of_points_Y(),
                                           this->m_grid.get_number_of_points_Z());
    //Then, we read the opacities meta file.
    //  This function also read each dustkappa_* file for each specie name.
    //  It's going to remap and interpolate the readed values according to the frequencies domain
    this->m_dust.read_opacities_meta(this->m_frequencies.get_frequencies());
    //The last process for the dust is to initialize the temperatures of each specie in the grid
    this->m_dust.initialize_specie_temperature(this->m_grid.get_number_of_points_X(),
                                               this->m_grid.get_number_of_points_Y(),
                                               this->m_grid.get_number_of_points_Z());

    //We calculate the temperatures DB. These values are precalculated temperatures that we are going to use in the simulation.
    this->m_emissivity.generate_emissivity_table(simulation_parameters,
                                                 this->m_dust.get_number_of_dust_species(),
                                                 this->m_dust.get_dust_species(),
                                                 this->m_frequencies.get_frequencies(),
                                                 this->m_frequencies.get_mean_intensity());
    this->m_emissivity.compute_derivate(this->m_dust.get_number_of_dust_species(),
                                        simulation_parameters["ntemp"],
                                        this->m_frequencies.get_mean_intensity());
    //With this, we end the setup process, and we proceed to run the m_photons

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

    //We initialize the m_photons of the simulation
    this->initialize_cartesian_regular_photons(generator, uniform_zero_one_distribution, simulation_parameters["nphot"]);
    //We launch each photon
    this->launch_photons(generator, uniform_zero_one_distribution, scattering_mode, simulation_parameters["nphot"],
                         simulation_parameters["ntemp"]);
    //When the m_photons end, we calculate the temperature of each dust specie from the energy in each cell of the grid
    this->calculate_dust_temperature();
    //To end the process, we write the temperatures output file
    //TODO : Cambiar nombre métdo a dust_temperatures
    this->write_dust_temperatures_file();
}

void monte_carlo::write_dust_temperatures_file(){
    std::ofstream temperatures_file("dust_temperature.dat");
    for (int iSpecie = 0; iSpecie < m_dust.get_number_of_dust_species(); ++iSpecie) {
        for (int i = 0; i < m_grid.get_number_of_points_X(); ++i) {
            for (int j = 0; j < m_grid.get_number_of_points_Y(); ++j) {
                for (int k = 0; k < m_grid.get_number_of_points_Z(); ++k) {
                    temperatures_file << m_dust.get_dust_species()[iSpecie].get_temperature()[k][j][i] << "\n";
                }
            }
        }
    }
    temperatures_file.close();
}

void monte_carlo::calculate_dust_temperature(){
    double cumulated_energy, temp;
    //When the m_photons end, we calculate the temperature of each dust specie from the energy in each cell of the grid
    for (int iSpecie = 0; iSpecie < this->m_dust.get_number_of_dust_species(); ++iSpecie) {
        for (int i = 0; i < this->m_grid.get_number_of_points_X(); ++i) {
            for (int j = 0; j < this->m_grid.get_number_of_points_Y(); ++j) {
                for (int k = 0; k < this->m_grid.get_number_of_points_Z(); ++k) {
                    cumulated_energy = this->m_dust.get_dust_species()[iSpecie].get_cumulative_energy()[k][j][i] / (this->m_dust.get_dust_species()[iSpecie].get_densities()[k][j][i] * this->m_grid.get_cell_volume());
                    if (cumulated_energy <= 0){
                        m_dust.set_null_temperature(iSpecie, i, j, k);
                    }
                    else{
                        temp = this->m_emissivity.compute_dust_temp_energy(cumulated_energy, iSpecie);
                        this->m_dust.set_specie_temperature_at_position(iSpecie, i, j, k, temp);
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
        if(i%100000 == 0) {
            std::cout << "Launching photon: " << i + 1 << std::endl;
        }
        this->move(this->m_photons[i], generator, uniform_zero_one_distribution);
        while(this -> m_photons[i].get_on_grid_condition()){
            //TODO : nombre de la función
            //TODO : logica de la función, generar eventos (definir en el foton)
            if(this -> m_photons[i].get_is_scattering_condition()) {
                if(scattering_mode == 1){
                    //this->m_photons[i].get_random_direction(generator, uniform_zero_one_distribution);
                    this -> get_random_direction(m_photons[i],generator,uniform_zero_one_distribution);
                }
                else if(scattering_mode == 2){
                    int scattering_specie = this -> m_photons[i].find_specie_to_scattering(generator, uniform_zero_one_distribution, this -> m_dust.get_number_of_dust_species());
                    //this -> m_photons[i].get_henvey_greenstein_direction(generator, uniform_zero_one_distribution, this -> m_dust.get_dust_species()[scattering_specie].get_g_interpoled());
                    this -> get_henvey_greenstein_direction(this -> m_photons[i], generator, uniform_zero_one_distribution,this -> m_dust.get_dust_species()[scattering_specie].get_g_interpoled());
                }
            }
            else{
                this -> do_absorption_event(generator, uniform_zero_one_distribution, this -> m_photons[i], number_of_temperatures);
                //this -> m_photons[i].get_random_direction(generator, uniform_zero_one_distribution);
                this -> get_random_direction(m_photons[i],generator,uniform_zero_one_distribution);
            }
            this->m_photons[i].calculate_tau_path(generator, uniform_zero_one_distribution);
            this->move(this->m_photons[i], generator, uniform_zero_one_distribution);
        }
    }
}

void monte_carlo::get_henvey_greenstein_direction(photon &photon_i,std::mt19937 &generator,
                                             std::uniform_real_distribution<> &uniform_zero_one_distribution,
                                             const std::vector<double> &g) {
    //TODO : cambiar los pow cuando sean cuadrados
    int valuesOrientations[3] = {0,1,1};
    double g2,xi,dir_x,l2,dir_y,dir_z,linv;
    double g_value = g[photon_i.get_frequency_index()];
    xi = uniform_zero_one_distribution(generator);
    if (g_value != 0.0){
        g2 = std::pow(g_value,2);
        while (xi == 0.0 and g_value == 1.0){
            xi = uniform_zero_one_distribution(generator);
        }
        dir_x = (0.5/g_value) * (1.0 + g2 - std::pow(((1.0-g2)/(1.0 - g_value + 2 * g_value * xi)),2));
    }
    else{
        dir_x = 2.0 * xi - 1.0;
    }
    l2 = 2.0;
    while (l2 > 1.0){
        dir_y = 2.0 * uniform_zero_one_distribution(generator) - 1.0;
        dir_z = 2.0 * uniform_zero_one_distribution(generator) - 1.0;
        l2 = std::pow(dir_y,2) + std::pow(dir_z,2);
        if (l2 < 1.0e-4){
            l2 = 2.0;
        }
    }

    linv = std::sqrt((1.0-std::pow(dir_x,2))/l2);
    // Now normalize to sqrt(1-mu^2)  where mu = dirx
    dir_y = dir_y * linv;
    dir_z = dir_z * linv;

    if(std::abs(dir_x*dir_x+dir_y*dir_y+dir_z*dir_z - 1.0) > 1.0e-6){
        std::cerr<< "ERROR : Henyey-Greenstein direction vector not OK " << std::endl;
        exit(0);
    }

    // Now we rotate the vector

    double oldx = photon_i.get_direction()[0];
    double oldy = photon_i.get_direction()[1];
    double oldz = photon_i.get_direction()[2];

    //rotateVector
    double l = std::sqrt((oldx * oldx) + (oldy * oldy));
    double vx = l * dir_x - oldz * dir_z;
    double vy = dir_y;
    double vz = oldz * dir_x + l * dir_z;
    //dir_x = vx;
    //dir_z = vz;
    if (l > 1e-10) {
        double dx = oldx / l;
        double dy = oldy / l;
        dir_x = dx * vx - dy * vy;
        dir_y = dy * vx + dx * vy;
        dir_z = vz;
    }
    else{
        dir_x = vx;
        dir_y = vy;
        dir_z = vz;
    }

    common::chek_unity_vector(dir_x, dir_y, dir_z);
    //get orientations
    //obtain orientations. It is 0 (left,down) or 1 (right, up)
    int ix = std::floor(dir_x) + 1.0;
    int iy = std::floor(dir_y) + 1.0;
    int iz = std::floor(dir_z) + 1.0;

    std::vector<int> orientation(3);
    std::vector<double> direction(3);

    orientation[0] = valuesOrientations[ix];
    orientation[1] = valuesOrientations[iy];
    orientation[2] = valuesOrientations[iz];

    direction[0] = dir_x;
    direction[1] = dir_y;
    direction[2] = dir_z;

    photon_i.set_orientation(orientation);
    photon_i.set_direction(direction);
}

void monte_carlo::get_random_direction(photon& photon_i, std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution) {
    std::vector<int> values_orientations = {0, 1, 1};
    //We intialize this values just to avoid a warning
    double l2, dirx=0.0, diry=0.0, dirz=0.0, linv;
    int ix, iy, iz;
    bool equal_zero = true;
    l2 = 2.0;
    while (equal_zero) {
        while (l2 > 1.0) {
            dirx = 2.0 * uniform_zero_one_distribution(generator) - 1.0;
            diry = 2.0 * uniform_zero_one_distribution(generator) - 1.0;
            dirz = 2.0 * uniform_zero_one_distribution(generator) - 1.0;
            l2 = dirx * dirx + diry * diry + dirz * dirz;
            if (l2 < 1e-4) {
                l2 = 2.0;
            }
        }
        equal_zero = (dirx == 0.0) || (diry == 0.0)|| (dirz == 0.0);
    }

    //TODO : Verificar si el l2 = 0 afecta (ver codigo esteban)
    //TODO : Esteban hace que dirx diry y dirz no puedan ser 0
    linv = 1.0 / std::sqrt(l2);

    dirx = dirx * linv;
    diry = diry * linv;
    dirz = dirz * linv;

    common::chek_unity_vector(dirx, diry, dirz);

    ix = std::floor(dirx) + 1.0;
    iy = std::floor(diry) + 1.0;
    iz = std::floor(dirz) + 1.0;

    std::vector<int> orientation(3);
    std::vector<double> direction(3);

    orientation[0] = values_orientations[ix];
    orientation[1] = values_orientations[iy];
    orientation[2] = values_orientations[iz];

    direction[0] = dirx;
    direction[1] = diry;
    direction[2] = dirz;

    photon_i.set_orientation(orientation);
    photon_i.set_direction(direction);
}

void monte_carlo::do_absorption_event(std::mt19937& generator, std::uniform_real_distribution<>& uniform_real_distribution, photon& photon_i, int number_of_temperatures){
    int ix,iy,iz,ray_inu;
    std::vector<double> temp_local(m_dust.get_number_of_dust_species());
    ix = photon_i.get_grid_position()[0];
    iy = photon_i.get_grid_position()[1];
    iz = photon_i.get_grid_position()[2];
    this -> divide_absorved_energy(photon_i);
    this -> add_temperature_decoupled(photon_i);
    for (int i = 0; i < m_dust.get_number_of_dust_species(); ++i) {
        //this -> dust_specie_temperature[i] = dust_species_information[i].get_temperature()[iz][iy][ix];
        temp_local[i] = m_dust.get_dust_species()[i].get_temperature()[iz][iy][ix];
        //this->dust_specie_temperature[i] = temperatures[i][iz][iy][ix];
    }
    photon_i.set_dust_specie_temperature(temp_local);
    ray_inu = this -> pickRandomFreqDb(generator,uniform_real_distribution,photon_i,number_of_temperatures);
    photon_i.set_frequency_index(ray_inu);
}

int monte_carlo::pickRandomFreqDb(std::mt19937& generator, std::uniform_real_distribution<>& uniform_real_distribution, photon& photon_i, int number_of_temperatures){
    std::vector<double> enerCum(m_dust.get_number_of_dust_species() + 1);
    int numCumul = m_frequencies.get_number_frequency_points() + 1;
    int intplt = 1;
    double rn;
    int iSpec = 0;
    int iTemp = 0;
    int inuPick = 0;
    enerCum[0] = 0;
    for (int i = 1; i < (m_dust.get_number_of_dust_species() + 1); ++i) {
        enerCum[i] = enerCum[i - 1] + photon_i.get_dust_specie_energy()[i - 1];
    }
    if (m_dust.get_number_of_dust_species() > 1) {
        rn = uniform_real_distribution(generator);
        iSpec = common::hunt(enerCum, m_dust.get_number_of_dust_species() + 1, rn, m_dust.get_number_of_dust_species());
        if ((iSpec < 0 || (iSpec > m_dust.get_number_of_dust_species() - 1))) {
            std::cerr << "ERROR : Specie found out of range ..." << std::endl;
            exit(0);
        }
    }
    iTemp = common::hunt(m_emissivity.get_db_temp(), number_of_temperatures, photon_i.get_temp_local()[iSpec], number_of_temperatures);
    double eps = 0.0;

    if (iTemp >= number_of_temperatures - 1) {
        std::cerr << "ERROR : Too high temperature discovered" << std::endl;
        exit(0);
    }
    if (iTemp <= -1) {
        iTemp = 0;
    } else {
        eps = (photon_i.get_temp_local()[iSpec] - m_emissivity.get_db_temp()[iTemp]) / (m_emissivity.get_db_temp()[iTemp + 1] - m_emissivity.get_db_temp()[iTemp]);
        if ((eps > 1) || (eps < 0)) {
            std::cerr << "ERROR : In picking new random frequency, eps out of range" << std::endl;
            exit(0);
        }
    }
    std::vector<double> dbCumul(m_frequencies.get_number_frequency_points() + 1);
    //TODO : intplt es siempre 1
    if (intplt == 1) {
        for (int inu = 0; inu < numCumul; ++inu) {
            dbCumul[inu] = (1.0 - eps) * m_emissivity.get_db_cumulnorm()[iSpec][iTemp][inu] + eps * m_emissivity.get_db_cumulnorm()[iSpec][iTemp + 1][inu];
        }
        rn = uniform_real_distribution(generator);
        inuPick = common::hunt(dbCumul, m_frequencies.get_number_frequency_points(), (double) rn, m_frequencies.get_number_frequency_points());
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
        inuPick = common::hunt(m_emissivity.get_db_cumulnorm()[iSpec][iTemp - 1], m_frequencies.get_number_frequency_points(),
                               rn, m_frequencies.get_number_frequency_points());
    }
    return inuPick;
}

void monte_carlo::add_temperature_decoupled(photon& photon_i){
    int ix,iy,iz;
    ix = photon_i.get_grid_position()[0];
    iy = photon_i.get_grid_position()[1];
    iz = photon_i.get_grid_position()[2];
    double cumen, temperature;
    for (int iSpec = 0; iSpec < m_dust.get_number_of_dust_species(); ++iSpec) {
        cumen = m_dust.get_dust_species()[iSpec].get_cumulative_energy()[iz][iy][ix] / (m_dust.get_dust_species()[iSpec].get_densities()[iz][iy][ix] * m_grid.get_cell_volume());
        temperature = m_emissivity.compute_dust_temp_energy(cumen, iSpec);
        m_dust.set_specie_temperature_at_position(iSpec, ix, iy, iz, temperature);
    }
}

void monte_carlo::divide_absorved_energy(photon& photon_i){
    int ix,iy,iz;
    double alpha_A = 0;
    std::vector<double> ener_part(m_dust.get_number_of_dust_species());
    ix = photon_i.get_grid_position()[0];
    iy = photon_i.get_grid_position()[1];
    iz = photon_i.get_grid_position()[2];
    if (m_dust.get_number_of_dust_species() == 1) {
        ener_part[0] = m_stars.get_stars_information()[photon_i.get_star_source()].get_energy();
    }
    else {
        for (int i = 0; i < m_dust.get_number_of_dust_species(); ++i) {
            ener_part[i] = m_dust.get_dust_species()[i].get_densities()[iz][iy][ix] * m_dust.get_dust_species()[i].get_kappa_absorption_interpoled()[photon_i.get_frequency_index()];
            alpha_A += ener_part[i];

        }
        for (int i = 0; i < m_dust.get_number_of_dust_species(); ++i) {
            ener_part[i] = m_stars.get_stars_information()[photon_i.get_star_source()].get_energy() * ener_part[i] / alpha_A;
        }
    }
    photon_i.set_dust_specie_energy(ener_part);
}

void monte_carlo::move(photon& photon_i, std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution){
    double minor_distance, fraction, add_tmp, dum;
    bool carry_on = true;
    while (carry_on) {
        //First we obtain the actual ray and grid position of the photon
        photon_i.set_prev_ray_position();
        photon_i.set_prev_grid_position();
        //Then we calculate the minor distance to move to the next cell of the grid
        minor_distance = photon_i.advance_next_position(this -> m_grid.get_number_of_points_X(),
                                                        this -> m_grid.get_number_of_points_Y(),
                                                        this -> m_grid.get_number_of_points_Z(),
                                                        this -> m_grid.get_x_points(),
                                                        this -> m_grid.get_y_points(),
                                                        this -> m_grid.get_z_points());
        photon_i.calculate_opacity_coefficients(minor_distance,
                                                      this -> m_dust.get_number_of_dust_species(),
                                                      this -> m_dust.get_dust_species());
        if (photon_i.get_tau_path_gone() + photon_i.get_dtau() > photon_i.get_tau_path_total()) {
            fraction = (photon_i.get_tau_path_total() - photon_i.get_tau_path_gone()) / photon_i.get_dtau();
            photon_i.update_ray_position(fraction);
            //TODO : Cambiar nombre variable dum
            dum = (1.0 - photon_i.get_albedo()) * (photon_i.get_tau_path_total() - photon_i.get_tau_path_gone()) *
                  m_stars.get_stars_information()[photon_i.get_star_source()].get_energy() / photon_i.get_alpha_A_total();
            for (int i = 0; i < m_dust.get_number_of_dust_species(); ++i) {
                //add_tmp = dum * this->alpha_A_specie[i];
                add_tmp = dum * photon_i.get_alpha_A_specie()[i];
                m_dust.add_energy_specie(i, photon_i.get_grid_position(), add_tmp);
            }
            carry_on = false;
        } else {
            dum = (1.0 - photon_i.get_albedo()) * photon_i.get_dtau() * m_stars.get_stars_information()[photon_i.get_star_source()].get_energy() / photon_i.get_alpha_A_total();
            for (int i = 0; i < m_dust.get_number_of_dust_species(); ++i) {
                add_tmp = dum * photon_i.get_alpha_A_specie()[i];
                m_dust.add_energy_specie(i, photon_i.get_prev_grid_position(), add_tmp);
            }
            photon_i.update_tau_path_gone();
            carry_on = photon_i.get_on_grid_condition();
        }
    }
    double rn = uniform_zero_one_distribution(generator);
    //this->is_scattering = rn < this->albedo;
    photon_i.set_scattering_state(rn);
}

void monte_carlo::initialize_cartesian_regular_photons(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution, int number_of_photons){
    this -> m_photons.resize(number_of_photons);
    int star, ray_frequency;
    for (unsigned int i = 0; i < this -> m_photons.size(); ++i) {
        //We identify the star that is going to give us the photon
        star = m_stars.identify_star(generator, uniform_zero_one_distribution);
        //We obtain a random frequency for the photon
        ray_frequency = m_frequencies.get_random_frequency(generator, uniform_zero_one_distribution, m_stars.get_stars_information()[star].get_cumulative_spectrum());
        //We initialize a new photon
        this -> m_photons[i] = photon(generator,
                                      uniform_zero_one_distribution,
                                      m_dust.get_number_of_dust_species(),
                                      m_frequencies.get_number_frequency_points(),
                                      star,
                                      ray_frequency,
                                      m_stars.get_stars_information()[star].get_star_position());
        //We set the grid position of the photon
        this->m_photons[i].set_grid_position(this->m_grid.found_point_cartesian_regular(
                this->m_photons[i].get_ray_position()));
        //We get a random direction for the photon
        this ->get_random_direction(this -> m_photons[i],generator,uniform_zero_one_distribution);
    }
}


std::map<std::string,double> monte_carlo::read_main_file(void){
    //First we create a map with the default values
    std::map<std::string,double> simulation_parameters;
    //###PHOTONS###
    //Number of m_photons packages to be used in the Monte Carlo simulation
    simulation_parameters["nphot"] = 100000.0;
    //Number of m_photons packages to be used for the scattering Monte Carlo simulation
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

/*
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
    this -> m_grid.initialize_cartesian_regular();
    //We read the frequencies from the "wavelength_micron.inp" file, and we calculate
    //  the mean intesity.
    //Note that freq_nu = m_frequencies.get_frequencies()
    //  and freq_dnu = m_frequencies.get_mean_intensity()
    this -> m_frequencies.read_frequencies();
    this -> m_frequencies.calculate_mean_intensity();
    //We read the star sources from the "stars.inp" and we calculate the spectrum and its cumulative function,
    //  the luminosities, the energy, and we invoke the function "jitter_stars" in order to avoid the star being
    //  exactly in the position of a wall of the grid
    this -> m_stars.read_stars();
    this -> m_stars.calculate_spectrum(this -> m_frequencies.get_frequencies());
    this -> m_stars.calculate_cumulative_spectrum(this -> m_frequencies.get_mean_intensity());
    this -> m_stars.calculate_total_luminosities(this -> m_frequencies.get_mean_intensity());
    this -> m_stars.calculate_energy(simulation_parameters["nphot"]);
    this -> m_stars.fix_luminosities();
    this -> m_stars.jitter_stars(this -> m_grid.get_x_points(),this -> m_grid.get_y_points(),this -> m_grid.get_z_points());
    //We read the dust information.
    //First, we read the "dust_density.inp" file, to obtain the density of the species in the grid
    this -> m_dust.read_dust_species_density(this -> m_grid.get_number_of_points_X(),this -> m_grid.get_number_of_points_Y(),this -> m_grid.get_number_of_points_Z());
    //Then, we read the opacities meta file.
    //  This function also read each dustkappa_* file for each specie name.
    //  It's going to remap and interpolate the readed values according to the frequencies domain
    this -> m_dust.read_opacities_meta(this -> m_frequencies.get_frequencies());
    //The last process for the dust is to initialize the temperatures of each specie in the grid
    this -> m_dust.initialize_specie_temperature(this -> m_grid.get_number_of_points_X(),this -> m_grid.get_number_of_points_Y(),this -> m_grid.get_number_of_points_Z());
    //We calculate the temperatures DB. These values are precalculated temperatures that we are going to use in the simulation.
    this -> m_emissivity.generate_emissivity_table(simulation_parameters,this -> m_dust.get_number_of_dust_species(),this -> m_dust.get_dust_species(),this -> m_frequencies.get_frequencies(),this->m_frequencies.get_mean_intensity());
    this -> m_emissivity.compute_derivate(this -> m_dust.get_number_of_dust_species(),simulation_parameters["ntemp"],this -> m_frequencies.get_mean_intensity());
    //With this, we end the setup process, and we proceed to run the m_photons
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

    //Now we set up the m_photons
    //TODO : Cambiar la logica del programa
    //TODO : El foton sólo debería tener atributos, no métodos (o al menos no todos los métodos considerados hasta ahora)
    this -> m_photons.resize(simulation_parameters["nphot"]);
    for (unsigned int i = 0; i < this -> m_photons.size(); ++i) {
    //for (int i = 0; i < 10; ++i) {
        std::cout << "Launching photon: " << i+1 << std::endl;
        //TODO : New lleva el objeto al heap. Acorde a lo leído, el heap tiene más espacio, pero es necesario borrar
        this -> m_photons[i] = photon(generator,
                                    uniform_zero_one_distribution,
                                    this -> m_dust.get_number_of_dust_species(),
                                    this -> m_stars.get_number_of_stars(),
                                    this -> m_frequencies.get_number_frequency_points(),
                                    this -> m_stars.get_cumulative_luminosity(),
                                    this -> m_stars.get_stars_information());
        //TODO : calcular la posición del foton respecto a la pos del rasho
        this -> m_photons[i].set_grid_position(this -> m_grid.found_point_cartesian_regular(this -> m_photons[i].get_ray_position()));
        //TODO : Acá inicia el ciclo del camino completo, es saltar al sgte evento y verificar si estamos dentro de la grilla. Si está, hacer los eventos relacionados
        std::cout << "walking to the next event" << std::endl;
        this -> m_photons[i].walk_next_event(generator,
                                           uniform_zero_one_distribution,
                                           this -> m_dust.get_number_of_dust_species(),
                                           this -> m_dust.get_dust_species_to_change(),
                                           this -> m_stars.get_stars_information(),
                                           this -> m_grid.get_number_of_points_X(),
                                           this -> m_grid.get_number_of_points_Y(),
                                           this -> m_grid.get_number_of_points_Z(),
                                           this -> m_grid.get_x_points(),
                                           this -> m_grid.get_y_points(),
                                           this -> m_grid.get_z_points());
        //While the photon is on the grid
        while(this -> m_photons[i].get_on_grid_condition()){
            std::cout << "The photon " << i + 1 << " is on the grid" << std::endl;
            //We need to do a scattering event...
            if(this -> m_photons[i].get_is_scattering_condition()){
                std::cout << "The photon " << i + 1 << " is going to scatter" << std::endl;
                //If we do a scattering event, we need to know the type of scattering
                if(scattering_mode == 1){
                    this -> m_photons[i].get_random_direction(generator,uniform_zero_one_distribution);
                }
                if(scattering_mode == 2){
                    int scattering_specie = this -> m_photons[i].find_specie_to_scattering(generator,uniform_zero_one_distribution,this -> m_dust.get_number_of_dust_species());
                    this -> m_photons[i].get_henvey_greenstein_direction(generator,uniform_zero_one_distribution,this -> m_dust.get_dust_species()[scattering_specie].get_g_interpoled());
                }
            }
            //Or an absorption event
            else{
                std::cout << "The photon " << i + 1 << " is going to do absorption" << std::endl;
                this -> m_photons[i].do_absorption_event(generator,
                                                       uniform_zero_one_distribution,
                                                       this -> m_dust.get_number_of_dust_species(),
                                                       this -> m_dust.get_dust_species_to_change(),
                                                       this -> m_stars.get_stars_information(),
                                                       this -> m_emissivity.get_db_temp(),
                                                       this -> m_emissivity.get_db_logenertemp(),
                                                       this -> m_emissivity.get_db_enertemp(),
                                                       simulation_parameters["ntemp"],
                                                       this -> m_frequencies.get_number_frequency_points(),
                                                       this -> m_grid.get_cell_volume(),
                                                       this -> m_emissivity.get_db_cumulnorm());
                this -> m_photons[i].get_random_direction(generator,uniform_zero_one_distribution);
            }
            std::cout << "The photon " << i + 1 << " is generating new tau path" << std::endl;
            this -> m_photons[i].calculate_tau_path(generator,uniform_zero_one_distribution);
            std::cout << "The photon " << i + 1 << " is walking to the next event" << std::endl;
            this -> m_photons[i].move(generator,
                                               uniform_zero_one_distribution,
                                               this -> m_dust.get_number_of_dust_species(),
                                               this -> m_dust.get_dust_species_to_change(),
                                               this -> m_stars.get_stars_information(),
                                               this -> m_grid.get_number_of_points_X(),
                                               this -> m_grid.get_number_of_points_Y(),
                                               this -> m_grid.get_number_of_points_Z(),
                                               this -> m_grid.get_x_points(),
                                               this -> m_grid.get_y_points(),
                                               this -> m_grid.get_z_points());
        }
    }
}
*/
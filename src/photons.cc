#include <photons.hh>

photons::photons(void) {
    ;
}



void photons::initialize_photons(int number_of_photons,int number_of_species, int number_of_stars, int number_of_frequencies, const std::vector<double>& luminosities_cum, const std::vector<double>& star_energy, const std::vector<double>& star_position, const std::vector<int>& number_of_grid_points, const std::vector<double>& difference, const std::vector<double>& star_cumulative_spectrum){
    this -> number_of_photons = number_of_photons;
    this -> cumulative_energy_specie.resize(number_of_species,std::vector<std::vector<std::vector<double>>>(number_of_grid_points[2],std::vector<std::vector<double>>(number_of_grid_points[1],std::vector<double>(number_of_grid_points[0]))));
    this -> pothon_information.assign(number_of_photons, photon(number_of_species, number_of_stars, number_of_frequencies, luminosities_cum, star_energy, star_position, number_of_grid_points, difference, star_cumulative_spectrum));
    //Inicializar cada foton (internamente cada foton adquiere la info necesaria)
    //Inicializar vectores de calor de las especies
    ;
}

void photons::walk_full_path(int number_of_species, const std::vector<double>& star_energy,const std::vector<std::vector<std::vector<std::vector<double>>>>& densities, const std::vector<std::vector<double>>& kappa_A, const std::vector<std::vector<double>>& kappa_S, const std::vector<double>& star_energies, const std::vector<int>& number_of_grid_points, const std::vector<double>& grid_cell_walls_x,const std::vector<double>& grid_cell_walls_y,const std::vector<double>& grid_cell_walls_z, int scattering_mode, const std::vector<std::vector<double>>& species_g, double cell_volumes,const std::vector<double>& dbTemp, const std::vector<std::vector<double>>& dbLogEnerTemp, const std::vector<std::vector<double>>& dbEnerTemp, int number_of_temperatures,int number_of_frequencies, const std::vector<std::vector<std::vector<double>>>& dbCumulNorm){
    int specie_to_scattering;
    for (int i = 0; i < this -> number_of_photons; ++i) {
        this -> pothon_information[i].walk_next_event(number_of_species,densities,kappa_A,kappa_S, star_energies,number_of_grid_points,grid_cell_walls_x,grid_cell_walls_y,grid_cell_walls_z,this->cumulative_energy_specie);
        while (this -> pothon_information[i].get_on_grid_condition()){
            if (this -> pothon_information[i].get_is_scattering_condition()){
                if (scattering_mode == 1){
                    this -> pothon_information[i].get_random_direction();
                }
                else if (scattering_mode == 2){
                    specie_to_scattering = this -> pothon_information[i].find_specie_to_scattering(number_of_species);
                    const std::vector<double>& g = species_g[specie_to_scattering];
                    this -> pothon_information[i].getHenveyGreensteinDirection(g);
                }
            }
            else{
                this -> pothon_information[i].do_absorption_event(number_of_species,this->temperatures,this->cumulative_energy_specie,densities,cell_volumes,star_energies,kappa_A,dbTemp,dbLogEnerTemp,dbEnerTemp,number_of_temperatures,number_of_frequencies,dbCumulNorm);
                this -> pothon_information[i].get_random_direction();
            }
            this -> pothon_information[i].get_tau_path();
            this -> pothon_information[i].walk_next_event(number_of_species,densities,kappa_A,kappa_S, star_energies,number_of_grid_points,grid_cell_walls_x,grid_cell_walls_y,grid_cell_walls_z,this->cumulative_energy_specie);
        }
    }
}

photons::~photons() {
    ;
}



#include <photon.hh>

photon::photon(int number_of_species, int number_of_stars, int number_of_frequencies,
               const std::vector<double> &luminosities_cum, const std::vector<double> &star_energy,
               const std::vector<double> &star_position, const std::vector<int> &number_of_grid_points,
               const std::vector<double> &difference, const std::vector<double> &star_cumulative_spectrum) {
    //At first, we resize each vector in order to store relevant information for the photon
    this->orientation.resize(3);
    this->direction.resize(3);
    this->distance.resize(3);
    this->alpha_A_specie.resize(number_of_species);
    this->alpha_S_specie.resize(number_of_species);
    this->cumulative_alpha.resize(number_of_species + 1);
    this->enerPart.resize(number_of_species);
    this->enerCum.resize(number_of_species + 1);
    this->dbCumul.resize(number_of_species + 1);
    //TODO : Crear código para generar los valores aleatorios, conversar con rannou
    //Then, we identify the star source of the photon
    this->identify_star(number_of_stars, luminosities_cum);
    //TODO : Que es esta energía xd?
    //TODO : El código original la saca pero acorde al code de esteban no hace nada
    double photon_energy = star_energy[this->star_source];
    //We obtain the ray position of the photon in AU, and then we obtain the position of the photon in the grid
    this->ray_position[0] = star_position[0];
    this->ray_position[1] = star_position[1];
    this->ray_position[2] = star_position[2];
    this->grid_position[0] = this->found_point(ray_position[0], "cartesian", number_of_grid_points[0], difference[0]);
    this->grid_position[1] = this->found_point(ray_position[1], "cartesian", number_of_grid_points[1], difference[1]);
    this->grid_position[2] = this->found_point(ray_position[2], "cartesian", number_of_grid_points[2], difference[2]);
    //We get a random direction for the photon
    this->get_random_direction();
    //TODO : Crear código para generar los valores aleatorios, conversar con rannou
    //We get a random frequency for the photon
    //At this frequency, we know relevant information about the scattering
    this->get_random_frequency_inu(star_cumulative_spectrum, number_of_frequencies);
    //TODO : Qué es el tau path
    this->get_tau_path();
    //At first, the photon is on the grid
    this->on_grid = true;
}

void photon::identify_star(int number_of_stars, const std::vector<double> &luminosities_cum) {
    //TODO : Dejar los generadores aleatorios en la clase, se usan varias veces...o bien! generar una funcion en common no?
    //TDOO : podría generalizar el método.
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    int star_lum = dis(gen);
    //TODO : Obtener luminosities cum
    int star = common::hunt(luminosities_cum, number_of_stars + 1, star_lum, number_of_stars / 2);
    this->star_source = star;
}

int photon::found_point(double x, std::string type, int number_of_points, double difference) {
    if (type == "cartesian") {
        return std::floor(x / difference) + (number_of_points / 2);
    }
    if (type == "spherical") {
        return std::floor(x / difference);
    }
}

void photon::get_random_direction() {
    //TODO : Entender este vector
    //TODO : dejar este vector en el objeto
    std::vector<int> values_orientations = {0, 1, 1};
    //TODO : Generalizar
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double l2, dirx, diry, dirz, linv;
    int ix, iy, iz;
    bool equal_zero = true;
    l2 = 2.0;
    while (equal_zero) {
        while (l2 > 1.0) {
            dirx = 2.0 * dis(gen) - 1.0;
            diry = 2.0 * dis(gen) - 1, 0;
            dirz = 2.0 * dis(gen) - 1.0;
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

    //TODO : Por qué revisar si es un vector unitario?
    this->chek_unity_vector(dirx, diry, dirz);

    ix = std::floor(dirx) + 1.0;
    iy = std::floor(diry) + 1.0;
    iz = std::floor(dirz) + 1.0;

    this->orientation[0] = values_orientations[ix];
    this->orientation[1] = values_orientations[iy];
    this->orientation[2] = values_orientations[iz];

    this->direction[0] = dirx;
    this->direction[1] = diry;
    this->direction[2] = dirz;
}

void photon::chek_unity_vector(double x, double y, double z) {
    double module = std::sqrt(x * x + y * y + z * z);
    if (std::fabs(module - 1.0) > 1e-6) {
        std::cerr << "ERROR : Error unity vector " << x << " " << y << " " << z << std::endl;
        exit(0);
    }
}

void photon::get_random_frequency_inu(const std::vector<double> &star_cumulative_spectrum, int number_of_frequencies) {
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double rn = dis(gen);
    int ray_inu = common::hunt(star_cumulative_spectrum, number_of_frequencies + 1, rn, number_of_frequencies / 2);
    this->ray_inu = ray_inu;
}

void photon::get_tau_path() {
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double rn = dis(gen);
    this->tau_path_total = -1.0 * std::log(1.0 - rn);
    this->tau_path_gone = 0.0;
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

void photon::walk_next_event(int number_of_species,
                             const std::vector<std::vector<std::vector<std::vector<double>>>> &densities,
                             const std::vector<std::vector<double>> &kappa_A,
                             const std::vector<std::vector<double>> &kappa_S, const std::vector<double> &star_energies,
                             const std::vector<int> &number_of_points, const std::vector<double> &grid_cell_walls_x,
                             const std::vector<double> &grid_cell_walls_y, const std::vector<double> &grid_cell_walls_z,
                             std::vector<std::vector<std::vector<std::vector<double>>>> &cumulative_energy_specie) {
    //TODO : Generalizar
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double minor_distance, fraction, add_tmp;
    bool carry_on = true;
    while (carry_on) {
        //First we obtain the actual ray and grid position of the photon
        this->prev_ray_position[0] = this->ray_position[0];
        this->prev_ray_position[1] = this->ray_position[1];
        this->prev_ray_position[2] = this->ray_position[2];

        this->prev_grid_position[0] = this->grid_position[0];
        this->prev_grid_position[1] = this->grid_position[1];
        this->prev_grid_position[2] = this->grid_position[2];

        //Then we calculate the minor distance to move to the next cell of the grid
        minor_distance = this->advance_next_position(number_of_points, grid_cell_walls_x, grid_cell_walls_y,
                                                     grid_cell_walls_z);

        this->calculate_opacity_coefficients(minor_distance, number_of_species, densities, kappa_A, kappa_S);

        if (this->tau_path_gone + this->dtau > this->tau_path_total) {
            fraction = (this->tau_path_total - tau_path_gone) / this->dtau;

            this->ray_position[0] =
                    this->prev_ray_position[0] + fraction * (this->ray_position[0] - prev_ray_position[0]);
            this->ray_position[1] =
                    this->prev_ray_position[1] + fraction * (this->ray_position[0] - prev_ray_position[1]);
            this->ray_position[2] =
                    this->prev_ray_position[2] + fraction * (this->ray_position[0] - prev_ray_position[2]);

            this->grid_position[0] = this->prev_grid_position[0];
            this->grid_position[1] = this->prev_grid_position[1];
            this->grid_position[2] = this->prev_grid_position[2];

            double dum = (1.0 - this->albedo) * (this->tau_path_total - this->tau_path_gone) *
                         star_energies[this->star_source] / this->alpha_A_total;
            for (int i = 0; i < number_of_species; ++i) {
                add_tmp = dum * this->alpha_A_specie[i];
                cumulative_energy_specie[i][this->grid_position[2]][this->grid_position[1]][this->grid_position[0]] += add_tmp;
            }
            carry_on = false;
        } else {
            double dum = (1.0 - this->albedo) * this->dtau * star_energies[this->star_source] / this->alpha_A_total;
            for (int i = 0; i < number_of_species; ++i) {
                add_tmp = dum * this->alpha_A_specie[i];
                cumulative_energy_specie[i][this->grid_position[2]][this->grid_position[1]][this->grid_position[0]] += add_tmp;
            }
            this->tau_path_gone = this->tau_path_gone + this->dtau;
            carry_on = this->on_grid;
        }
    }
    double rn = dis(gen);
    this->is_scattering = rn < this->albedo;
}

double
photon::advance_next_position(const std::vector<int> &number_of_points, const std::vector<double> &grid_cell_walls_x,
                              const std::vector<double> &grid_cell_walls_y,
                              const std::vector<double> &grid_cell_walls_z) {

    std::vector<int> signs = {-1, 1};
    this->get_cell_walls(grid_cell_walls_x, grid_cell_walls_y, grid_cell_walls_z);
    this->distance[0] = (this->cell_walls[0] - this->ray_position[0]) / this->direction[0];
    this->distance[1] = (this->cell_walls[1] - this->ray_position[1]) / this->direction[1];
    this->distance[2] = (this->cell_walls[2] - this->ray_position[2]) / this->direction[2];

    double min_distance = std::min(std::min(distance[0], distance[1]), distance[2]);
    int count = 0;
    std::vector<int> indexes = {-1, -1, -1};
    for (int i = 0; i < 3; ++i) {
        if (this->distance[i] == min_distance) {
            indexes[count] = i;
            count++;
        }
    }

    this->ray_position[0] = this->ray_position[0] + min_distance * this->direction[0];
    this->ray_position[1] = this->ray_position[1] + min_distance * this->direction[1];
    this->ray_position[2] = this->ray_position[2] + min_distance * this->direction[2];

    //avoid bug assign cellWall to ray position
    //update grid position with signs
    for (int i = 0; i < count; i++) {
        this->ray_position[indexes[i]] = this->cell_walls[indexes[i]];
        this->grid_position[indexes[i]] = this->grid_position[indexes[i]] + signs[this->orientation[indexes[i]]];
    }

    this->is_on_grid(number_of_points[0], number_of_points[1], number_of_points[2]);
    return min_distance;
}

void photon::get_cell_walls(const std::vector<double> &grid_cell_walls_x, const std::vector<double> &grid_cell_walls_y,
                            const std::vector<double> &grid_cell_walls_z) {
    this->cell_walls[0] = grid_cell_walls_x[this->grid_position[0] + this->orientation[0]];
    this->cell_walls[1] = grid_cell_walls_y[this->grid_position[1] + this->orientation[1]];
    this->cell_walls[2] = grid_cell_walls_z[this->grid_position[2] + this->orientation[2]];
}

void photon::is_on_grid(int number_of_points_x, int number_of_points_y, int number_of_points_z) {
    bool on_x = (this->grid_position[0] >= 0) && (this->grid_position[0] < number_of_points_x);
    bool on_y = (this->grid_position[1] >= 0) && (this->grid_position[1] < number_of_points_y);
    bool on_z = (this->grid_position[2] >= 0) && (this->grid_position[2] < number_of_points_z);
    this->on_grid = on_x && on_y && on_z;
}

void photon::calculate_opacity_coefficients(double minor_distance, int number_of_species,
                                            std::vector<std::vector<std::vector<std::vector<double>>>> densities,
                                            std::vector<std::vector<double>> kappa_A,
                                            std::vector<std::vector<double>> kappa_S) {
    //We obatin the previos position of the photon before the movement
    int ix = this->prev_grid_position[0];
    int iy = this->prev_grid_position[0];
    int iz = this->prev_grid_position[0];

    //In order to obtain the total opacity of the cell, we need to considerate the densities of the species as the
    //  position of the photon
    double opacity_coefficient_alpha_A_total = 0;
    double opacity_coefficient_alpha_S_total = 0;

    for (int i = 0; i < number_of_species; ++i) {
        //We obtain the absorption opacity and the scattering opacity and then we do that value times the
        //  density and we accumulate the result
        this->alpha_A_specie[i] = densities[i][iz][iy][ix] * kappa_A[i][this->ray_inu];
        this->alpha_S_specie[i] = densities[i][iz][iy][ix] * kappa_S[i][this->ray_inu];
        opacity_coefficient_alpha_A_total = opacity_coefficient_alpha_A_total + this->alpha_A_specie[i];
        opacity_coefficient_alpha_S_total = opacity_coefficient_alpha_S_total + this->alpha_S_specie[i];
    }
    //We store the calculated values for the photon
    this->alpha_A_total = opacity_coefficient_alpha_A_total;
    this->alpha_S_total = opacity_coefficient_alpha_S_total;
    this->alpha_total = opacity_coefficient_alpha_A_total + opacity_coefficient_alpha_S_total;
    //The albedo is the radiation percentage that the specie reflex respect the radiation that affects it
    this->albedo = opacity_coefficient_alpha_S_total / this->alpha_total;
    this->dtau = this->alpha_total * minor_distance;

}

void photon::do_absorption_event(int number_of_species,
                                 std::vector<std::vector<std::vector<std::vector<double>>>> &temperatures,
                                 const std::vector<std::vector<std::vector<std::vector<double>>>> &cumulEner,
                                 const std::vector<std::vector<std::vector<std::vector<double>>>> &densities,
                                 double cellVolumes, const std::vector<double> &star_energies,
                                 const std::vector<std::vector<double>> &kappa_A, const std::vector<double> &dbTemp,
                                 const std::vector<std::vector<double>> &dbLogEnerTemp,
                                 const std::vector<std::vector<double>> &dbEnerTemp, int number_of_temperatures,
                                 int number_of_frequencies,
                                 const std::vector<std::vector<std::vector<double>>> &dbCumulNorm) {
    int ix = this->grid_position[0];
    int iy = this->grid_position[1];
    int iz = this->grid_position[2];
    this->divideAbsorvedEnergy(number_of_species, star_energies, densities, kappa_A);
    this->addTemperatureDecoupled(number_of_species, cumulEner, densities, cellVolumes, temperatures, dbTemp,
                                  dbLogEnerTemp, dbEnerTemp, number_of_temperatures);
    for (int i = 0; i < number_of_species; ++i) {
        this->tempLocal[i] = temperatures[i][iz][iy][ix];
    }
    this->pickRandomFreqDb(number_of_frequencies, number_of_species, number_of_temperatures, dbTemp, dbCumulNorm);
}

void photon::pickRandomFreqDb(int number_of_frequencies, int number_of_species, int number_of_temperatures,
                              const std::vector<double> &dbTemp,
                              const std::vector<std::vector<std::vector<double>>> &dbCumulNorm) {
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    //float sumTemp1, sumTemp2;
    int numCumul = number_of_frequencies + 1;
    int intplt = 1;
    double rn;
    int iSpec = 0;
    int iTemp = 0;
    int inuPick = 0;
    this->enerCum[0] = 0;
    for (int i = 1; i < (number_of_species + 1); ++i) {
        this->enerCum[i] = this->enerCum[i - 1] + enerPart[i - 1];
    }
    if (number_of_species > 1) {
        rn = dis(gen);
        iSpec = common::hunt(this->enerCum, number_of_species + 1, rn, enerCum[number_of_species / 2]);
        if ((iSpec < 0 || (iSpec > number_of_species - 1))) {
            std::cerr << "ERROR : Specie found out of range ..." << std::endl;
            exit(0);
            //printf("exit\n");
            //exit(1);
        }
    }
    iTemp = common::hunt(dbTemp, number_of_temperatures, this->tempLocal[iSpec], dbTemp[number_of_temperatures / 2]);
    double eps = 0.0;

    if (iTemp >= number_of_temperatures - 1) {
        //printf("ERROR: Too high temperature discovered\n");
        std::cerr << "ERROR : Too high temperature discovered" << std::endl;
        exit(0);
    }
    if (iTemp <= -1) {
        iTemp = 0;
    } else {
        eps = (tempLocal[iSpec] - dbTemp[iTemp]) / (dbTemp[iTemp + 1] - dbTemp[iTemp]);
        if ((eps > 1) || (eps < 0)) {
            std::cerr << "ERROR : In picking new random frequency, eps out of range" << std::endl;
            exit(0);
            //printf("ERROR: in pick randomFreq, eps<0 or eps>1\n");
            //printf("exit\n");
            //exit(1);
        }
    }

    if (intplt == 1) {
        for (int inu = 0; inu < numCumul; ++inu) {
            this->dbCumul[inu] =
                    (1.0 - eps) * dbCumulNorm[iSpec][iTemp][inu] + eps * dbCumulNorm[iSpec][iTemp + 1][inu];
        }
        rn = dis(gen);
        //rn = 0.20160651049169737;
        inuPick = common::hunt(this->dbCumul, number_of_frequencies, (double) rn, dbCumul[number_of_frequencies / 2]);
        //free(dbCumul);
    } else {
        //Now, within this species/size, find the frequency
        if (eps > 0.5) {
            if (iTemp < number_of_temperatures - 2) {
                iTemp++;
            }
        }
        rn = dis(gen);
        //verify dbCumulNorm
        inuPick = common::hunt(dbCumulNorm[iSpec][iTemp - 1], number_of_frequencies, rn,
                               dbCumulNorm[iSpec][iTemp - 1][number_of_frequencies / 2]);
    }
    this->ray_inu = inuPick;
}

void photon::addTemperatureDecoupled(int number_of_species,
                                     const std::vector<std::vector<std::vector<std::vector<double>>>> &cumulEner,
                                     const std::vector<std::vector<std::vector<std::vector<double>>>> &densities,
                                     double cellVolumes,
                                     std::vector<std::vector<std::vector<std::vector<double>>>> &temperatures,
                                     const std::vector<double> &dbTemp,
                                     const std::vector<std::vector<double>> &dbLogEnerTemp,
                                     const std::vector<std::vector<double>> &dbEnerTemp, int number_of_temperatures) {
    int ix = this->grid_position[0];
    int iy = this->grid_position[1];
    int iz = this->grid_position[2];
    double cumen;
    for (int iSpec = 0; iSpec < number_of_species; ++iSpec) {
        //TODO : Obtener cumulEner de la especie ispec y densidad en ispec y cell
        cumen = cumulEner[iSpec][iz][iy][ix] / (densities[iSpec][iz][iy][ix] * cellVolumes);
        //TODO : Este vector inicia en 0 para cada especie en el objeto de polvo...
        //TODO : Hay que cachar si puedo cargar en cada protón y luego sumar todo o depende de las temps anteriores
        //const std::vector<double>& dbTemp, std::vector<std::vector<double>>& dbLogEnerTemp, const std::vector<std::vector<double>>& dbEnerTemp, int number_of_temperatures, double energy, int iSpec
        temperatures[iSpec][iz][iy][ix] = this->computeDusttempEnergyBd(dbTemp, dbLogEnerTemp, dbEnerTemp,
                                                                        number_of_temperatures, cumen, iSpec);
        //dustTemperature->temperatures[iSpec][iz][iy][ix] = computeDusttempEnergyBd(emissivityDb, cumen, iSpec);
    }
}

double photon::computeDusttempEnergyBd(const std::vector<double> &dbTemp,
                                       const std::vector<std::vector<double>> &dbLogEnerTemp,
                                       const std::vector<std::vector<double>> &dbEnerTemp, int number_of_temperatures,
                                       double energy, int iSpec) {
    //printf("in computeDusttempEnergyBd\n");
    double tempReturn = 0.0;
    int itemp = 0;
    double logEner = std::log(energy);
    double eps;
    //TODO : Obtener dblogenergtemp de la especie
    itemp = common::hunt(dbLogEnerTemp[iSpec], number_of_temperatures, (double) logEner, number_of_temperatures / 2);
    //printf("itemp=%d\n", *itemp);
    if (itemp >= number_of_temperatures - 1) {
        std::cerr << "ERROR : Too high temperature discovered" << std::endl;
        exit(0);
    }

    if (itemp <= -1) {
        //Temperature presumably below lowest temp in dbase
        //TODO : Obtener enertemp de la especie
        eps = energy / dbEnerTemp[iSpec][0];
        if (eps >= 1) {
            //printf("exit\n");
            //exit(1);
            std::cerr << "ERROR : Too high temperature discovered" << std::endl;
            exit(0);
        }
        tempReturn = eps * dbTemp[0];
    } else {
        //TODO : Obtener enertemp de la especie
        eps = (energy - dbEnerTemp[iSpec][itemp]) / (dbEnerTemp[iSpec][itemp + 1] - dbEnerTemp[iSpec][itemp]);
        if ((eps > 1) || (eps < 0)) {
            std::cerr << "ERROR : Temperature found out of range..." << std::endl;
            printf("exit\n");
            //exit(1);
        }
        //TODO : Obtener dbTemp
        tempReturn = (1.0 - eps) * dbTemp[itemp] + eps * dbTemp[itemp + 1];
    }
    return tempReturn;
}

void photon::divideAbsorvedEnergy(int number_of_species, const std::vector<double> &star_energies,
                                  const std::vector<std::vector<std::vector<std::vector<double>>>> &density,
                                  const std::vector<std::vector<double>> &kappa_A) {
    //printf("in divideAbsorvedEnergy\n");
    int ix = this->grid_position[0];
    int iy = this->grid_position[1];
    int iz = this->grid_position[2];
    double alphaA = 0;
    if (number_of_species == 1) {
        //TODO : Obtener energia de la estrella
        this->enerPart[0] = star_energies[this->star_source];
    } else {
        for (int i = 0; i < number_of_species; ++i) {
            //TODO : Obtener densidad de la especie i
            this->enerPart[i] = density[i][iz][iy][ix] * kappa_A[i][this->ray_inu];
            //photon->enerPart[i] = dustDensity->densities[i][iz][iy][ix] * dustOpacity->kappaA[i][photon->iFrequency];
            alphaA += this->enerPart[i];
        }
        for (int i = 0; i < number_of_species; ++i) {
            this->enerPart[i] = star_energies[this->star_source] * this->enerPart[i] / alphaA;
        }
    }
}

void photon::getHenveyGreensteinDirection(const std::vector<double> &g) {
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    int valuesOrientations[3] = {0, 1, 1};
    double newx, newy, newz;
    double temp, g2, sign;
    bool equalZero = true;
    //printf("dustOpacity->g[iSpec][photon->iFrequency]=%10.10lg\n",dustOpacity->g[iSpec][photon->iFrequency]);
    //TODO : Obtener g de la especie obtenida
    double g_value = g[this->ray_inu];

    //get random numbers != 0.5 , and 2*random-1 != 0
    double rnX = 0.5;
    double rnY = 0.5;
    double rnZ = 0.5;

    while (rnX != 0.5) {
        rnX = dis(gen);
    }

    newx = 2 * rnX - 1.0;
    if (g_value > 0) {
        g2 = g_value * g_value;
        temp = (1.0 - g2) / (1.0 + g_value * newx);
        newx = (1.0 + g2 - temp * temp) / (2.0 * g_value);
        newx = std::max(newx, -1.0);
        newx = std::min(newx, 1.0);
    }

    float l2 = 2.0;
    while (l2 > 1.0) {
        rnY = dis(gen);//curand_uniform(&photon->state);
        rnZ = dis(gen);//curand_uniform(&photon->state);
        newy = 2 * rnY - 1;
        newz = 2 * rnZ - 1;
        l2 = newy * newy + newz * newz;
        if (l2 < 0.0001) {
            l2 = 2.0;
        }
        //l2 = (l2 < 0.0001 ? 2.0 : l2);
    }

    double linv = std::sqrt((1.0 - newx * newx) / l2);
    //float linv = sqrtf((1.0-newx*newx)/l2);
    newy = newy * linv;
    newz = newz * linv;

    double oldx = this->direction[0];
    double oldy = this->direction[1];
    double oldz = this->direction[2];

    //rotateVector
    double l = std::sqrt((oldx * oldx) + (oldx * oldy));
    //float l = sqrtf(oldx*oldx+oldy*oldy);
    double vx = l * newx - oldz * newz;
    double vy = newy;
    double vz = oldz * newx + l * newz;
    newx = vx;
    newz = vz;
    if (l > 1e-6) {
        double dx = oldx / l;
        double dy = oldy / l;
        newx = dx * vx - dy * vy;
        newy = dy * vx + dx * vy;
    }
    this->chek_unity_vector(newx, newy, newz);
    //get orientations
    //obtain orientations. It is 0 (left,down) or 1 (right, up)
    int ix = std::floor(newx) + 1.0;
    int iy = std::floor(newy) + 1.0;
    int iz = std::floor(newz) + 1.0;

    this->orientation[0] = valuesOrientations[ix];
    this->orientation[1] = valuesOrientations[iy];
    this->orientation[2] = valuesOrientations[iz];
    //printf("floor: %d, %d, %d\n",ix,iy,iz);
    /*if (ix==2 || iy==2 ||iz==2 ){
      printf("actual: dirx=%lf diry=%lf dirz=%lf\nnew: dirx=%lf diry=%lf dirz=%lf\n",photon->direction[0],photon->direction[1],photon->direction[2],newx,newy,newz);
    }

    if (ix==2 || iy==2 ||iz==2 ){
      checkUnitVector(newx,newy,newz);
      //printf("newx=%2.8lg newy=%2.8lg newz=%2.8lg\n",newx,newy,newz);
    }*/
    this->direction[0] = newx;
    this->direction[1] = newy;
    this->direction[2] = newz;
}

int photon::find_specie_to_scattering(int number_of_species) {

    //TODO : Esto pareciera ser general a todos los fotones, si es el caso, dejarlo en la clase superior
    this->cumulative_alpha[0] = 0.0;

    for (int i = 0; i < number_of_species; ++i) {
        this->cumulative_alpha[i + 1] = this->cumulative_alpha[i] + this->alpha_S_specie[i];
    }
    for (int i = 0; i < number_of_species; ++i) {
        this->cumulative_alpha[i] = this->cumulative_alpha[i] / this->cumulative_alpha[number_of_species];
    }

    this->cumulative_alpha[number_of_species] = 1.0;

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    double rn = dis(gen);
    int iSpec = common::hunt(this->cumulative_alpha, number_of_species + 1, (double) rn, number_of_species / 2);
    return iSpec;

}

int photon::get_star_source(void) {
    return this->star_source;
}

void photon::set_star_source(int star_source) {
    this->star_source = star_source;
}

bool photon::get_on_grid_condition() {
    return this->on_grid;
}

bool photon::get_is_scattering_condition() {
    return this->is_scattering;
}

photon::~photon(void) {
    ;
}

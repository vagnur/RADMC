#include <photon.hh>

photon::photon(void){
    ;
}

photon::photon(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
               int number_of_species, int number_of_frequencies, int star_source, int frequency_index, const std::vector<double>& star_ray_position,
               int number_of_points_x, int number_of_points_y, int number_of_points_z) {
    //At first, we resize each vector in order to store relevant information for the photon
    this->orientation.resize(3);
    this->direction.resize(3);
    this->distance.resize(3);
    this->ray_position.resize(3);
    this->grid_position.resize(3);
    this -> prev_grid_position.resize(3);
    this -> prev_ray_position.resize(3);
    this->alpha_A_specie.resize(number_of_species);
    this->alpha_S_specie.resize(number_of_species);
    this->cumulative_alpha.resize(number_of_species + 1);
    this->dust_specie_energy.resize(number_of_species);
    this->cumulative_dust_specie_energy.resize(number_of_species + 1);
    this -> cell_walls.resize(3);
    this -> dust_specie_temperature.resize(number_of_species);
    this -> dbCumul.resize(number_of_frequencies+1);
    //Then, we identify the star source of the photon
    this -> star_source = star_source;
    //We set the ray position according to the star source
    this -> set_ray_position(star_ray_position);
    //We get a random frequency for the photon
    //At this frequency, we know relevant information about the scattering properties of the dust specie
    this -> frequency_index = frequency_index;
    //We get a tau path for the photon
    this->calculate_tau_path(generator, uniform_zero_one_distribution);
    //We calculate if the photon is inside or outside the grid
    //this ->is_on_grid(number_of_points_x, number_of_points_y, number_of_points_z);
    //std::cout << this -> on_grid << std::endl;
}

void photon::calculate_tau_path(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution) {
    double rn = uniform_zero_one_distribution(generator);
    this->tau_path_total = -1.0 * std::log(1.0 - rn);
    this->tau_path_gone = 0.0;
}

void photon::is_on_grid(int number_of_points_x, int number_of_points_y, int number_of_points_z) {
    bool on_x = (this->grid_position[0] >= 0) && (this->grid_position[0] < number_of_points_x);
    bool on_y = (this->grid_position[1] >= 0) && (this->grid_position[1] < number_of_points_y);
    bool on_z = (this->grid_position[2] >= 0) && (this->grid_position[2] < number_of_points_z);
    this->on_grid = on_x && on_y && on_z;
}

void photon::calculate_opacity_coefficients(int number_of_species, const std::vector<dust_species>& dust_species_information) {
    //We obtain the previous position of the photon before the movement
    int ix = this->prev_grid_position[0];
    int iy = this->prev_grid_position[1];
    int iz = this->prev_grid_position[2];

    //In order to obtain the total opacity of the cell, we need to considerate the densities of the species as the
    //  position of the photon
    double opacity_coefficient_alpha_A_total = 0;
    double opacity_coefficient_alpha_S_total = 0;

    for (int i = 0; i < number_of_species; ++i) {
        //We obtain the absorption opacity and the scattering opacity, then we do that value times the
        //  density, and we accumulate the result

        //this->alpha_A_specie[i] = densities[i][iz][iy][ix] * kappa_A[i][this->frequency_index];
        //this->alpha_S_specie[i] = densities[i][iz][iy][ix] * kappa_S[i][this->frequency_index];
        this -> alpha_A_specie[i] = dust_species_information[i].get_densities()[iz][ix][iy] * dust_species_information[i].get_kappa_absorption_interpoled()[this->frequency_index];
        this -> alpha_S_specie[i] = dust_species_information[i].get_densities()[iz][ix][iy] * dust_species_information[i].get_kappa_scattering_interpoled()[this->frequency_index];
        opacity_coefficient_alpha_A_total = opacity_coefficient_alpha_A_total + this->alpha_A_specie[i];
        opacity_coefficient_alpha_S_total = opacity_coefficient_alpha_S_total + this->alpha_S_specie[i];
    }
    //We store the calculated values for the photon
    this->alpha_A_total = opacity_coefficient_alpha_A_total;
    this->alpha_S_total = opacity_coefficient_alpha_S_total;
    this->alpha_total = opacity_coefficient_alpha_A_total + opacity_coefficient_alpha_S_total;
    //The albedo is the radiation percentage that the specie reflex respect the radiation that affects it
    this->albedo = opacity_coefficient_alpha_S_total / this->alpha_total;
    this->dtau = this->alpha_total * this -> min_distance;

}

int photon::find_specie_to_scattering(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
                                      int number_of_species) {
    this->cumulative_alpha[0] = 0.0;
    for (int i = 0; i < number_of_species; ++i) {
        this->cumulative_alpha[i + 1] = this->cumulative_alpha[i] + this->alpha_S_specie[i];
    }
    for (int i = 0; i < number_of_species; ++i) {
        this->cumulative_alpha[i] = this->cumulative_alpha[i] / this->cumulative_alpha[number_of_species];
    }
    this->cumulative_alpha[number_of_species] = 1.0;

    double rn = uniform_zero_one_distribution(generator);
    int iSpec = common::hunt(this->cumulative_alpha, number_of_species + 1, rn, number_of_species);
    return iSpec;
}

void photon::update_ray_position(double fraction){
    this->ray_position[0] =
            this->prev_ray_position[0] + fraction * (this->ray_position[0] - prev_ray_position[0]);
    this->ray_position[1] =
            this->prev_ray_position[1] + fraction * (this->ray_position[1] - prev_ray_position[1]);
    this->ray_position[2] =
            this->prev_ray_position[2] + fraction * (this->ray_position[2] - prev_ray_position[2]);

    this->grid_position[0] = this->prev_grid_position[0];
    this->grid_position[1] = this->prev_grid_position[1];
    this->grid_position[2] = this->prev_grid_position[2];
}

void photon::update_tau_path_gone(){
    this->tau_path_gone = this->tau_path_gone + this->dtau;
}

int photon::get_star_source(void) const {
    return this->star_source;
}

bool photon::get_on_grid_condition() const {
    return this->on_grid;
}

std::string photon::get_next_event() const {
    if(this -> is_scattering){
        return "scattering";
    }
    else{
        return "absorption";
    }
}

const std::vector<double> &photon::get_ray_position() const {
    return ray_position;
}

double photon::get_dtau() const{
    return this -> dtau;
}

double photon::get_tau_path_total() const{
    return this -> tau_path_total;
}

double photon::get_tau_path_gone() const{
    return this -> tau_path_gone;
}

double photon::get_albedo() const{
    return this -> albedo;
}

double photon::get_alpha_A_total() const{
    return this -> alpha_A_total;
}

const std::vector<double>& photon::get_alpha_A_specie() const{
    return this -> alpha_A_specie;
}

const std::vector<int>& photon::get_grid_position() const{
    return this -> grid_position;
}

int photon::get_frequency_index() const{
    return this -> frequency_index;
}

const std::vector<double>& photon::get_dust_specie_energy() const{
    return this -> dust_specie_energy;
}

const std::vector<double>& photon::get_cumulative_dust_specie_energy() const{
    return this -> cumulative_dust_specie_energy;
}

const std::vector<double>& photon::get_temp_local() const{
    return this -> dust_specie_temperature;
}

const std::vector<int>& photon::get_prev_grid_position() const{
    return this -> prev_grid_position;
}

const std::vector<double>& photon::get_direction() const{
    return this -> direction;
}

const std::vector<int>& photon::get_orientation() const{
    return this -> orientation;
}

void photon::set_orientation(const std::vector<int>& orientation){
    this -> orientation[0] = orientation[0];
    this -> orientation[1] = orientation[1];
    this -> orientation[2] = orientation[2];
}

void photon::set_direction(const std::vector<double>& direction){
    this -> direction[0] = direction[0];
    this -> direction[1] = direction[1];
    this -> direction[2] = direction[2];
}

void photon::set_dust_specie_energy(const std::vector<double>& dust_specie_energy){
    this -> dust_specie_energy = dust_specie_energy;
}

void photon::set_dust_specie_temperature(const std::vector<double>& dust_specie_temperature){
    this -> dust_specie_temperature = dust_specie_temperature;
}

void photon::set_grid_position(const std::vector<int>& grid_position) {
    this -> grid_position = grid_position;
}

void photon::set_prev_ray_position(){
    this->prev_ray_position[0] = this->ray_position[0];
    this->prev_ray_position[1] = this->ray_position[1];
    this->prev_ray_position[2] = this->ray_position[2];
}

void photon::set_prev_grid_position() {
    this->prev_grid_position[0] = this->grid_position[0];
    this->prev_grid_position[1] = this->grid_position[1];
    this->prev_grid_position[2] = this->grid_position[2];
}

void photon::set_ray_position(std::vector<double> ray_position){
    this -> ray_position[0] = ray_position[0];
    this -> ray_position[1] = ray_position[1];
    this -> ray_position[2] = ray_position[2];
}

void photon::set_frequency_index(int frequency_index){
    this -> frequency_index = frequency_index;
}

void photon::set_scattering_state(double rn){
    this->is_scattering = rn < this->albedo;
}

photon::~photon(void) {
    ;
}

void photon::set_walls(const std::vector<double> &cell_walls) {
    this -> cell_walls = cell_walls;
}

const std::vector<double> &photon::get_cell_walls() const {
    return this -> cell_walls;
}

void photon::set_distance(const std::vector<double> distance) {
    this -> distance = distance;
}

void photon::set_min_distance(double min_distance) {
    this -> min_distance = min_distance;
}



//////////////////DEPRECIATED CODE////////////////////

/*
double
photon::advance_next_position(int number_of_points_X,int number_of_points_Y,int number_of_points_Z, const std::vector<double> &grid_cell_walls_x,
                              const std::vector<double> &grid_cell_walls_y,
                              const std::vector<double> &grid_cell_walls_z) {

    std::vector<int> signs = {-1, 1};
    this->obtain_cell_walls(grid_cell_walls_x, grid_cell_walls_y, grid_cell_walls_z);
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
        //TODO : Revisar signo +
        this->grid_position[indexes[i]] = this->grid_position[indexes[i]] + signs[this->orientation[indexes[i]]];
    }

    this->is_on_grid(number_of_points_X, number_of_points_Y, number_of_points_Z);
    return min_distance;
}
 */

/*
void photon::obtain_cell_walls(const std::vector<double> &grid_cell_walls_x, const std::vector<double> &grid_cell_walls_y,
                               const std::vector<double> &grid_cell_walls_z) {
    this->cell_walls[0] = grid_cell_walls_x[this->grid_position[0] + this->orientation[0]];
    this->cell_walls[1] = grid_cell_walls_y[this->grid_position[1] + this->orientation[1]];
    this->cell_walls[2] = grid_cell_walls_z[this->grid_position[2] + this->orientation[2]];
}
*/

/*
photon::photon(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
               int number_of_species, int number_of_stars, int number_of_frequencies,
               const std::vector<double> &luminosities_cum, const std::vector<star>& star_information) {
    //At first, we resize each vector in order to store relevant information for the photon
    this->orientation.resize(3);
    this->direction.resize(3);
    this->distance.resize(3);
    this->ray_position.resize(3);
    this->grid_position.resize(3);
    this -> prev_grid_position.resize(3);
    this -> prev_ray_position.resize(3);
    this->alpha_A_specie.resize(number_of_species);
    this->alpha_S_specie.resize(number_of_species);
    this->cumulative_alpha.resize(number_of_species + 1);
    this->dust_specie_energy.resize(number_of_species);
    this->cumulative_dust_specie_energy.resize(number_of_species + 1);
    this -> cell_walls.resize(3);
    this -> dust_specie_temperature.resize(number_of_species);
    this -> dbCumul.resize(number_of_frequencies+1);
    //TODO : Crear código para generar los valores aleatorios, conversar con rannou
    //Then, we identify the star source of the photon
    this->identify_star(generator,uniform_zero_one_distribution,number_of_stars, luminosities_cum);
    //TODO : Que es esta energía xd?
    //TODO : El código original la saca pero acorde al code de esteban no hace nada
    //TODO : Simplemente es la energía de la estrella. A priori no es necesario dejarlo en el fotón.
    //double photon_energy = star_energy[this->star_source];
    //double photon_energy = star_information[this -> star_source].get_energy();
    //We obtain the ray position of the photon in AU, and then we obtain the position of the photon in the grid
    this->ray_position[0] = star_information[this -> star_source].get_star_position()[0];
    this->ray_position[1] = star_information[this -> star_source].get_star_position()[1];
    this->ray_position[2] = star_information[this -> star_source].get_star_position()[2];
    //this->ray_position[0] = star_position[0];
    //this->ray_position[1] = star_position[1];
    //this->ray_position[2] = star_position[2];
    //this->grid_position[0] = this->found_point(ray_position[0], "cartesian", number_of_grid_points[0], difference[0]);
    //this->grid_position[1] = this->found_point(ray_position[1], "cartesian", number_of_grid_points[1], difference[1]);
    //this->grid_position[2] = this->found_ray_position_in_grid(ray_position[2], "cartesian", number_of_grid_points[2], difference[2]);
    //We get a random direction for the photon
    this->get_random_direction(generator, uniform_zero_one_distribution);
    //TODO : Crear código para generar los valores aleatorios, conversar con rannou
    //We get a random frequency for the photon
    //At this frequency, we know relevant information about the scattering properties of the dust specie
    this->get_random_frequency_inu(generator, uniform_zero_one_distribution,star_information[this -> star_source].get_cumulative_spectrum(), number_of_frequencies);
    //TODO : Qué es el tau path
    this->calculate_tau_path(generator,uniform_zero_one_distribution);
    //At first, the photon is on the grid
    this->on_grid = true;
}
 */

/*
void photon::identify_star(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
                           int number_of_stars, const std::vector<double> &luminosities_cum) {
    double star_lum = uniform_zero_one_distribution(generator);
    //TODO : Obtener luminosities cum
    int star = common::hunt(luminosities_cum, number_of_stars + 1, star_lum, number_of_stars);
    this->star_source = star;
}
 */

/*
int photon::found_ray_position_in_grid(double x, std::string type, int number_of_points, double difference) {
    return std::floor(x / difference) + (number_of_points / 2);
}
 */

/*
void photon::check_unity_vector(double x, double y, double z) {
    double module = std::sqrt((x * x) + (y * y) + (z * z));
    if (std::fabs(module - 1.0) > 1e-6) {
        std::cerr << "ERROR : Error unity vector " << x << " " << y << " " << z << std::endl;
        exit(0);
    }
}
 */

/*
void photon::get_random_frequency_inu(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
                                      const std::vector<double> &star_cumulative_spectrum, int number_of_frequencies) {
    double rn = uniform_zero_one_distribution(generator);
    int ray_inu = common::hunt(star_cumulative_spectrum, number_of_frequencies + 1, rn, number_of_frequencies);
    this->frequency_index = ray_inu;
}
 */

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
        this -> calculate_tau_path();
        this -> move_photon();
    }
}
*/

/*
void photon::walk_next_event(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,int number_of_species,
                             std::vector<dust_species>& dust_specie_information, const std::vector<star>& stars_information,
                             int number_of_points_X, int number_of_points_Y, int number_of_points_Z, const std::vector<double> &grid_cell_walls_x,
                             const std::vector<double> &grid_cell_walls_y, const std::vector<double> &grid_cell_walls_z) {
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

        //Then we calculate the minor distance to move_photon to the next cell of the grid
        minor_distance = this->advance_next_position(number_of_points_X, number_of_points_Y, number_of_points_Z, grid_cell_walls_x, grid_cell_walls_y,
                                                     grid_cell_walls_z);
        this->calculate_opacity_coefficients(minor_distance, number_of_species, dust_specie_information);
        if (this->tau_path_gone + this->dtau > this->tau_path_total) {
            fraction = (this->tau_path_total - tau_path_gone) / this->dtau;

            this->ray_position[0] =
                    this->prev_ray_position[0] + fraction * (this->ray_position[0] - prev_ray_position[0]);
            this->ray_position[1] =
                    this->prev_ray_position[1] + fraction * (this->ray_position[1] - prev_ray_position[1]);
            this->ray_position[2] =
                    this->prev_ray_position[2] + fraction * (this->ray_position[2] - prev_ray_position[2]);

            this->grid_position[0] = this->prev_grid_position[0];
            this->grid_position[1] = this->prev_grid_position[1];
            this->grid_position[2] = this->prev_grid_position[2];

            //TODO : Acá está la energia
            double dum = (1.0 - this->albedo) * (this->tau_path_total - this->tau_path_gone) *
                         stars_information[this->star_source].get_energy() / this->alpha_A_total;
            for (int i = 0; i < number_of_species; ++i) {
                add_tmp = dum * this->alpha_A_specie[i];
                //TODO : como recibimos y cambiamos el vector?
                dust_specie_information[i].add_energy(this -> grid_position[0],this -> grid_position[1],this -> grid_position[2],add_tmp);
                //cumulative_energy_specie[i][this->grid_position[2]][this->grid_position[1]][this->grid_position[0]] += add_tmp;
            }
            carry_on = false;
        } else {
            double dum = (1.0 - this->albedo) * this->dtau * stars_information[this->star_source].get_energy() / this->alpha_A_total;
            for (int i = 0; i < number_of_species; ++i) {
                add_tmp = dum * this->alpha_A_specie[i];
                dust_specie_information[i].add_energy(this -> prev_grid_position[0],this -> prev_grid_position[1],this -> prev_grid_position[2],add_tmp);
                //cumulative_energy_specie[i][this->grid_position[2]][this->grid_position[1]][this->grid_position[0]] += add_tmp;
            }
            this->tau_path_gone = this->tau_path_gone + this->dtau;
            carry_on = this->on_grid;
        }
    }
    double rn = uniform_zero_one_distribution(generator);
    this->is_scattering = rn < this->albedo;
}
 */

/*
void photon::do_absorption_event(std::mt19937& generator, std::uniform_real_distribution<>& uniform_real_distribution,
                                 int number_of_species,
                                 std::vector<dust_species>& dust_species_information,
                                 const std::vector<star>& star_information,
                                 const std::vector<double> &dbTemp,
                                 const std::vector<std::vector<double>> &dbLogEnerTemp,
                                 const std::vector<std::vector<double>> &dbEnerTemp, int number_of_temperatures,
                                 int number_of_frequencies, double cellVolumes,
                                 const std::vector<std::vector<std::vector<double>>> &dbCumulNorm) {
    int ix = this->grid_position[0];
    int iy = this->grid_position[1];
    int iz = this->grid_position[2];
    this -> divideAbsorvedEnergy(number_of_species, star_information, dust_species_information);
    this -> addTemperatureDecoupled(number_of_species,cellVolumes,dust_species_information,dbTemp,dbLogEnerTemp,dbEnerTemp,number_of_temperatures);
    for (int i = 0; i < number_of_species; ++i) {
        this -> dust_specie_temperature[i] = dust_species_information[i].get_temperature()[iz][iy][ix];
        //this->dust_specie_temperature[i] = temperatures[i][iz][iy][ix];
    }
    this->generate_random_frequency(generator, uniform_real_distribution,number_of_frequencies, number_of_species, number_of_temperatures, dbTemp, dbCumulNorm);

}
 */

/*
void photon::generate_random_frequency(std::mt19937& generator, std::uniform_real_distribution<>& uniform_real_distribution,
                              int number_of_frequencies, int number_of_species, int number_of_temperatures,
                              const std::vector<double> &dbTemp,
                              const std::vector<std::vector<std::vector<double>>> &dbCumulNorm) {
    //float sumTemp1, sumTemp2;
    int numCumul = number_of_frequencies + 1;
    int intplt = 1;
    double rn;
    int iSpec = 0;
    int iTemp = 0;
    int inuPick = 0;
    this->cumulative_dust_specie_energy[0] = 0;
    for (int i = 1; i < (number_of_species + 1); ++i) {
        this->cumulative_dust_specie_energy[i] = this->cumulative_dust_specie_energy[i - 1] + this -> dust_specie_energy[i - 1];
    }
    if (number_of_species > 1) {
        rn = uniform_real_distribution(generator);
        iSpec = common::hunt(this->cumulative_dust_specie_energy, number_of_species + 1, rn, number_of_species);
        if ((iSpec < 0 || (iSpec > number_of_species - 1))) {
            std::cerr << "ERROR : Specie found out of range ..." << std::endl;
            //exit(0);
            //printf("exit\n");
            //exit(1);
        }
    }
    iTemp = common::hunt(dbTemp, number_of_temperatures, this->dust_specie_temperature[iSpec], number_of_temperatures);
    double eps = 0.0;

    if (iTemp >= number_of_temperatures - 1) {
        //printf("ERROR: Too high temperature discovered\n");
        std::cerr << "ERROR : Too high temperature discovered" << std::endl;
        //exit(0);
    }
    if (iTemp <= -1) {
        iTemp = 0;
    } else {
        eps = (dust_specie_temperature[iSpec] - dbTemp[iTemp]) / (dbTemp[iTemp + 1] - dbTemp[iTemp]);
        if ((eps > 1) || (eps < 0)) {
            std::cerr << "ERROR : In picking new random frequency, eps out of range" << std::endl;
            //exit(0);
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
        rn = uniform_real_distribution(generator);
        //rn = 0.20160651049169737;
        inuPick = common::hunt(this->dbCumul, number_of_frequencies, (double) rn, number_of_frequencies);
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
        inuPick = common::hunt(dbCumulNorm[iSpec][iTemp - 1], number_of_frequencies, rn,
                               number_of_frequencies);
    }
    this->frequency_index = inuPick;
}
 */

/*
//void photon::addTemperatureDecoupled(int number_of_species,
//                                     const std::vector<std::vector<std::vector<std::vector<double>>>> &cumulEner,
//                                     const std::vector<std::vector<std::vector<std::vector<double>>>> &densities,
//                                     double cellVolumes,
//                                     std::vector<std::vector<std::vector<std::vector<double>>>> &temperatures,
//                                     const std::vector<double> &dbTemp,
//                                     const std::vector<std::vector<double>> &dbLogEnerTemp,
//                                     const std::vector<std::vector<double>> &dbEnerTemp, int number_of_temperatures) {
void photon::addTemperatureDecoupled(int number_of_species, double cellVolumes,
                                     std::vector<dust_species>& dust_species_information,
                                     const std::vector<double> &dbTemp,
                                     const std::vector<std::vector<double>> &dbLogEnerTemp,
                                     const std::vector<std::vector<double>> &dbEnerTemp, int number_of_temperatures) {
    int ix = this->grid_position[0];
    int iy = this->grid_position[1];
    int iz = this->grid_position[2];
    double cumen, temperature;
    for (int iSpec = 0; iSpec < number_of_species; ++iSpec) {
        //TODO : Obtener cumulEner de la especie ispec y densidad en ispec y cell
        //cumen = cumulEner[iSpec][iz][iy][ix] / (densities[iSpec][iz][iy][ix] * cellVolumes);
        cumen = dust_species_information[iSpec].get_cumulative_energy()[iz][iy][ix] / (dust_species_information[iSpec].get_densities()[iz][iy][ix]*cellVolumes);
        //TODO : Este vector inicia en 0 para cada especie en el objeto de polvo...
        //TODO : Hay que cachar si puedo cargar en cada protón y luego sumar todo o depende de las temps anteriores
        //const std::vector<double>& dbTemp, std::vector<std::vector<double>>& dbLogEnerTemp, const std::vector<std::vector<double>>& dbEnerTemp, int number_of_temperatures, double energy, int iSpec
        temperature = this ->computeDusttempEnergyBd(dbTemp, dbLogEnerTemp, dbEnerTemp, number_of_temperatures, cumen,iSpec);
        dust_species_information[iSpec].set_temperature_at_position(ix,iy,iz,temperature);
        //temperatures[iSpec][iz][iy][ix] = this->computeDusttempEnergyBd(dbTemp, dbLogEnerTemp, dbEnerTemp,
        //                                                               number_of_temperatures, cumen, iSpec);
        //dustTemperature->temperatures[iSpec][iz][iy][ix] = computeDusttempEnergyBd(emissivityDb, cumen, iSpec);
    }
}
 */

/*
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
    itemp = common::hunt(dbLogEnerTemp[iSpec], number_of_temperatures, logEner, number_of_temperatures);
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
            //printf("exit\n");
            exit(0);
        }
        //TODO : Obtener dbTemp
        tempReturn = (1.0 - eps) * dbTemp[itemp] + eps * dbTemp[itemp + 1];
    }
    return tempReturn;
}
 */

/*
//void photon::divideAbsorvedEnergy(int number_of_species, const std::vector<double> &star_energies,
//                                  const std::vector<std::vector<std::vector<std::vector<double>>>> &density,
//                                  const std::vector<std::vector<double>> &kappa_A) {
void photon::divideAbsorvedEnergy(int number_of_species, const std::vector<star>& star_information,
                                  const std::vector<dust_species>& dust_specie_information){
    //printf("in divideAbsorvedEnergy\n");
    int ix = this->grid_position[0];
    int iy = this->grid_position[1];
    int iz = this->grid_position[2];
    double alphaA = 0;
    if (number_of_species == 1) {
        //TODO : Obtener energia de la estrella
        this -> dust_specie_energy[0] = star_information[this -> star_source].get_energy();
        //this->dust_specie_energy[0] = star_energies[this->star_source];
    }
    else {
        for (int i = 0; i < number_of_species; ++i) {
            //TODO : Obtener densidad de la especie i
            this -> dust_specie_energy[i] = dust_specie_information[i].get_densities()[iz][iy][ix] * dust_specie_information[i].get_kappa_absorption_interpoled()[this -> frequency_index];
            //this->dust_specie_energy[i] = density[i][iz][iy][ix] * kappa_A[i][this->frequency_index];
            //photon->dust_specie_energy[i] = dustDensity->densities[i][iz][iy][ix] * dustOpacity->kappaA[i][photon->iFrequency];
            alphaA += this->dust_specie_energy[i];
        }
        for (int i = 0; i < number_of_species; ++i) {
            this -> dust_specie_energy[i] = star_information[this -> star_source].get_energy() * this -> dust_specie_energy[i] / alphaA;
            //this->dust_specie_energy[i] = star_energies[this->star_source] * this->dust_specie_energy[i] / alphaA;
        }
    }
}
 */

/*
void photon::get_henvey_greenstein_direction(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
                                          const std::vector<double> &g) {
    int valuesOrientations[3] = {0, 1, 1};
    double newx, newy, newz;
    double temp, g2;
    //printf("dustOpacity->g[iSpec][photon->iFrequency]=%10.10lg\n",dustOpacity->g[iSpec][photon->iFrequency]);
    //TODO : Obtener g de la especie obtenida
    double g_value = g[this->frequency_index];
    std::cout << this -> frequency_index << std::endl;
    std::cout << g_value << std::endl;

    //get random numbers != 0.5 , and 2*random-1 != 0
    double rnX = 0.5;
    double rnY = 0.5;
    double rnZ = 0.5;

    while (rnX != 0.5) {
        rnX = uniform_zero_one_distribution(generator);
    }

    newx = 2 * rnX - 1.0;
    std::cout << newx << std::endl;
    if (g_value > 0) {
        g2 = g_value * g_value;
        temp = (1.0 - g2) / (1.0 + g_value * newx);
        newx = (1.0 + g2 - temp * temp) / (2.0 * g_value);
        newx = std::max(newx, -1.0);
        newx = std::min(newx, 1.0);
    }
    std::cout << newx << std::endl;

    double l2 = 2.0;
    while (l2 > 1.0) {
        rnY = uniform_zero_one_distribution(generator);//curand_uniform(&photon->state);
        rnZ = uniform_zero_one_distribution(generator);//curand_uniform(&photon->state);
        newy = 2 * rnY - 1.0;
        newz = 2 * rnZ - 1.0;
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
    common::check_unity_vector(newx, newy, newz);
    //get orientations
    //obtain orientations. It is 0 (left,down) or 1 (right, up)
    int ix = std::floor(newx) + 1.0;
    int iy = std::floor(newy) + 1.0;
    int iz = std::floor(newz) + 1.0;

    this->orientation[0] = valuesOrientations[ix];
    this->orientation[1] = valuesOrientations[iy];
    this->orientation[2] = valuesOrientations[iz];
    //printf("floor: %d, %d, %d\n",ix,iy,iz);
    if (ix==2 || iy==2 ||iz==2 ){
      printf("actual: dirx=%lf diry=%lf dirz=%lf\nnew: dirx=%lf diry=%lf dirz=%lf\n",photon->direction[0],photon->direction[1],photon->direction[2],newx,newy,newz);
    }

    if (ix==2 || iy==2 ||iz==2 ){
      checkUnitVector(newx,newy,newz);
      //printf("newx=%2.8lg newy=%2.8lg newz=%2.8lg\n",newx,newy,newz);
    }
    this->direction[0] = newx;
    this->direction[1] = newy;
    this->direction[2] = newz;
}
*/

/*
void photon::get_random_direction(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution) {
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

    common::check_unity_vector(dirx, diry, dirz);

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
*/

/*
void photon::get_henvey_greenstein_direction(std::mt19937 &generator,
                                             std::uniform_real_distribution<> &uniform_zero_one_distribution,
                                             const std::vector<double> &g) {
    //TODO : cambiar los pow cuando sean cuadrados
    int valuesOrientations[3] = {0,1,1};
    double g2,xi,dir_x,l2,dir_y,dir_z,linv;
    double g_value = g[this -> frequency_index];
    g_value = 0.12445256783491307;
    xi = uniform_zero_one_distribution(generator);
    if (g_value != 0.0){
        g2 = std::pow(g_value,2);
        while (xi == 0.0 and g_value == 1.0){
            xi = uniform_zero_one_distribution(generator);
        }
        xi = 0.17741459161887413;
        dir_x = (0.5/g_value) * (1.0 + g2 - std::pow(((1.0-g2)/(1.0 - g_value + 2 * g_value * xi)),2));
    }
    else{
        dir_x = 2.0 * xi - 1.0;
    }
    l2 = 2.0;
    while (l2 > 1.0){
        dir_y = 0.74997921613727536;
        //dir_y = 2.0 * uniform_zero_one_distribution(generator) - 1.0;
        dir_y = 2.0 * dir_y - 1.0;
        dir_z = 0.30334996416120863;
        //dir_z = 2.0 * uniform_zero_one_distribution(generator) - 1.0;
        dir_z = 2.0 * dir_z - 1.0;
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
    }

    //TODO : Todo ok hasta aca

    // Rotacion vector
    //std::cout << dir_x << " " << dir_y << " " << dir_z << std::endl;
    //std::cout << std::sqrt((dir_x*dir_x)+(dir_y*dir_y)+(dir_z*dir_z)) << std::endl;

    double oldx = this->direction[0];
    double oldy = this->direction[1];
    double oldz = this->direction[2];

    oldx = -0.15734385689173966;
    oldy = -1.3453880951885377e-002;
    oldz = -0.98745222860944748;

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

    common::check_unity_vector(dir_x, dir_y, dir_z);
    //get orientations
    //obtain orientations. It is 0 (left,down) or 1 (right, up)
    int ix = std::floor(dir_x) + 1.0;
    int iy = std::floor(dir_y) + 1.0;
    int iz = std::floor(dir_z) + 1.0;

    this->orientation[0] = valuesOrientations[ix];
    this->orientation[1] = valuesOrientations[iy];
    this->orientation[2] = valuesOrientations[iz];

    this -> direction[0] = dir_x;
    this -> direction[1] = dir_y;
    this -> direction[2] = dir_z;
}
*/
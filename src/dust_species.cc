#include <dust_species.hh>

dust_species::dust_species(void){
    ;
}

void dust_species::initialize_temperature(int number_of_points_x,int number_of_points_y,int number_of_points_z) {
    std::vector<std::vector<std::vector<double>>> temperatures(number_of_points_z,std::vector<std::vector<double>>(number_of_points_y,std::vector<double>(number_of_points_x,0)));
    this -> temperatures = temperatures;
    std::vector<std::vector<std::vector<double>>> cumulative_energy(number_of_points_z,std::vector<std::vector<double>>(number_of_points_y,std::vector<double>(number_of_points_x,0)));
    this -> cumulative_energy = cumulative_energy;
}

void dust_species::set_density(const std::vector<std::vector<std::vector<double>>>& densities){
    this -> densities = densities;
}

void dust_species::set_lambda(const std::vector<double>& lambda){
    this -> lambda = lambda;
}

void dust_species::set_kappa_absorption(const std::vector<double>& kappa_absorption){
    this -> kappa_absorption = kappa_absorption;
}

void dust_species::set_kappa_scattering(const std::vector<double>& kappa_scattering){
    this -> kappa_scattering = kappa_scattering;
}

void dust_species::set_g(const std::vector<double>& g){
    this -> g = g;
}

void dust_species::set_frequency(const std::vector<double>& frequency){
    this -> frequency = frequency;
}

const std::vector<std::vector<std::vector<double>>>& dust_species::get_densities() const{
    return this -> densities;
}

const std::vector<double>& dust_species::get_frequency() const{
    return this -> frequency;
}

const std::vector<double>& dust_species::get_absoprtion() const{
    return this -> kappa_absorption;
}

const std::vector<double>& dust_species::get_scattering() const{
    return this -> kappa_scattering;
}

const std::vector<double>& dust_species::get_g() const{
    return this -> g;
}

const std::vector<double>& dust_species::get_lambda() const{
    return this -> lambda;
}

dust_species::~dust_species(void){
    ;
}

const std::vector<double>& dust_species::get_kappa_absorption_remapped() const {
    return this -> kappa_absorption_remapped;
}

void dust_species::set_kappa_absorption_remapped(const std::vector<double>& kappa_absorption_remapped) {
    this -> kappa_absorption_remapped = kappa_absorption_remapped;
}

const std::vector<double>& dust_species::get_kappa_scattering_remapped() const {
    return this -> kappa_scattering_remapped;
}

void dust_species::set_kappa_scattering_remapped(const std::vector<double>& kappa_scattering_remapped) {
    this -> kappa_scattering_remapped = kappa_scattering_remapped;
}

const std::vector<double>& dust_species::get_g_remapped() const {
    return this -> g_remapped;
}

void dust_species::set_g_remapped(const std::vector<double>& g_remapped) {
    this -> g_remapped = g_remapped;
}

const std::vector<double>& dust_species::get_kappa_absorption_interpoled() const {
    return this -> kappa_absorption_interpoled;
}

void dust_species::set_kappa_absorption_interpoled(const std::vector<double>& kappa_absorption_interpoled) {
    this -> kappa_absorption_interpoled = kappa_absorption_interpoled;
}

const std::vector<double>& dust_species::get_kappa_scattering_interpoled() const {
    return this -> kappa_scattering_interpoled;
}

void dust_species::set_kappa_scattering_interpoled(const std::vector<double>& kappa_scattering_interpoled) {
    this -> kappa_scattering_interpoled = kappa_scattering_interpoled;
}

const std::vector<double>& dust_species::get_g_interpoled() const {
    return this -> g_interpoled;
}

void dust_species::set_g_interpoled(const std::vector<double>& g_interpoled) {
    this -> g_interpoled = g_interpoled;
}

void dust_species::add_energy(int pos_X, int pos_Y, int pos_Z, double add_tmp) {
    this -> cumulative_energy[pos_Z][pos_Y][pos_X] += add_tmp;
}

const std::vector<std::vector<std::vector<double>>> &dust_species::get_cumulative_energy() const {
    return this -> cumulative_energy;
}

const std::vector<std::vector<std::vector<double>>> &dust_species::get_temperature() const{
    return this -> temperatures;
};

void dust_species::set_temperature_at_position(int ix, int iy, int iz, double temperature){
    this -> temperatures[iz][iy][ix] = temperature;
}

void dust_species::set_null_temperature(int ix, int iy, int iz){
    this -> temperatures[iz][iy][ix] = 0;
}
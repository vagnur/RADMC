#include <dust_species.hh>

dust_species::dust_species(void){
    ;
}

void dust_species::set_density(std::vector<std::vector<std::vector<double>>> densities){
    this -> densities = densities;
}

void dust_species::set_lambda(std::vector<double> lambda){
    this -> lambda = lambda;
}

void dust_species::set_kappa_absorption(std::vector<double> kappa_absorption){
    this -> kappa_absorption = kappa_absorption;
}

void dust_species::set_kappa_scattering(std::vector<double> kappa_scattering){
    this -> kappa_scattering = kappa_scattering;
}

void dust_species::set_g(std::vector<double> g){
    this -> g = g;
}

void dust_species::set_frequency(std::vector<double> frequency){
    this -> frequency = frequency;
}

std::vector<std::vector<std::vector<double>>> dust_species::get_densities(){
    return this -> densities;
}

std::vector<double> dust_species::get_frequency(){
    return this -> frequency;
}

std::vector<double> dust_species::get_absoprtion(){
    return this -> kappa_absorption;
}

std::vector<double> dust_species::get_scattering(){
    return this -> kappa_scattering;
}

std::vector<double> dust_species::get_g(){
    return this -> g;
}

dust_species::~dust_species(void){
    ;
}
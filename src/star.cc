#include <star.hh>

star::star(void){
    ;
}

star::star(double radio, double mass, double x_position, double y_position, double z_position){
    this -> radio = radio;
    this -> mass = mass;
    this -> positions.push_back(x_position);
    this -> positions.push_back(y_position);
    this -> positions.push_back(z_position);
}

void star::set_flux(const std::vector<double>& flux){
    this -> flux = flux;
}

void star::set_spectrum(const std::vector<double>& spectrum){
    this -> spectrum = spectrum;
}

void star::set_cumulative_spectrum(const std::vector<double>& cumulative_spectrum){
    this -> cumulative_spectrum = cumulative_spectrum;
}

const std::vector<double>& star::get_flux(void) const{
    return this -> flux;
}

const std::vector<double>& star::get_spectrum(void) const{
    return this -> spectrum;
}

double star::get_star_radio(void) const{
    return this -> radio;
}

double star::get_luminosity(void) const{
    return this -> luminosity;
}

void star::set_luminosity(double luminosity){
    this -> luminosity = luminosity;
}

void star::set_energy(double star_energy){
    this -> energy = star_energy;
}

const std::vector<double>& star::get_star_position(void) const{
    return this -> positions;
}

void star::set_star_position(const std::vector<double>& star_position){
    this -> positions = star_position;
}

double star::get_energy(void) const{
    return this -> energy;
}

star::~star(void){
    ;
}

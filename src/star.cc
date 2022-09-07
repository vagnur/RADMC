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

void star::set_flux(std::vector<double> flux){
    this -> flux = flux;
}

void star::set_spectrum(std::vector<double> spectrum){
    this -> spectrum = spectrum;
}

void star::set_cumulative_spectrum(std::vector<double> cumulative_spectrum){
    this -> cumulative_spectrum = cumulative_spectrum;
}

std::vector<double> star::get_flux(){
    return this -> flux;
}

std::vector<double> star::get_spectrum(){
    return this -> spectrum;
}

double star::get_star_radio(){
    return this -> radio;
}

double star::get_luminosity(void){
    return this -> luminosity;
}

void star::set_luminosity(double luminosity){
    this -> luminosity = luminosity;
}

star::~star(void){
    ;
}

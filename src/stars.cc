#include <stars.hh>

stars::stars(void){
    ;
}

void stars::read_stars(void){
    //Read the input file for the stars information
    std::ifstream input_file;
    input_file.open("inputs/stars.inp");
    //In case that the file couldn't be opened, an error message is displayed.
    if(!input_file){
        std::cerr << "Mandatory file \"stars.inp\" could not be opened. Make sure that the file exists" << std::endl;
        exit(0);
    }
    int number_of_stars,iformat,number_of_frequencies;
    double star_radio,star_mass,x_position,y_position,z_position,lambda_point;
    //First we read the format of the file.
    //TODO : Ahora mismo siempre se trabaja con el valor 2. El manual indica que si el valor es 1 entonces
    //TODO : los valores de la frecuencia serán valores directamente en hertz y no en micrones
    input_file >> iformat;
    this -> iformat = iformat;
    //Then we read the number of stars and the number of frequency points
    input_file >> number_of_stars;
    this -> number_of_stars = number_of_stars;
    //We resize the vector that is going to alocate each star
    this -> stars_information.resize(number_of_stars);
    input_file >> number_of_frequencies;
    this -> number_of_frequencies = number_of_frequencies;
    //Now we create a new star object with the information relevant to that star
    for (int i = 0; i < number_of_stars; ++i) {
        input_file >> star_radio;
        input_file >> star_mass;
        input_file >> x_position;
        input_file >> y_position;
        input_file >> z_position;
        star str(star_radio,star_mass,x_position,y_position,z_position);
        this -> stars_information[i] = str;
        //TODO : ¿Es necesario borrar al estrella antes de generar la nueva?
        //TODO : Creo que sería bueno borrar la memoria de la estrella en el destructor (estudiar y ver rendimiento)
    }
    //TODO : En el codigo original verifican varias cosas sobre la grilla y la estrella.
    //TODO : Esto implica que hay que leer la grilla antes
    //Now we read the frequency points from the file
    //TODO : estas frecuencias son las mismas que ya se leyeron xd
    //TODO : el programa original compara si son las mismas o no y se detiene si no lo son
    //TODO : evaluar si es necesario, a priori se considera que no y no se trabaja con los valores
    for (int i = 0; i < number_of_frequencies; ++i) {
        input_file >> lambda_point;
    }
    //We read the flux of each star
    std::vector<double> star_flux;
    double star_flux_value;
    for (int i = 0; i < number_of_stars; ++i) {
        //If the first flux present in the file is negative, then we need to read only that value
        input_file >> star_flux_value;
        if(star_flux_value < 0){
            star_flux.push_back(star_flux_value);
            this -> stars_information[i].set_flux(star_flux);
            star_flux.clear();
        }
        //If not, then we read the full spectrum
        else {
            star_flux.resize(number_of_frequencies);
            star_flux[0] = star_flux_value;
            for (int j = 1; j < number_of_frequencies; ++j) {
                input_file >> star_flux_value;
                star_flux[j] = star_flux_value;
            }
            this -> stars_information[i].set_flux(star_flux);
            star_flux.clear();
        }
    }
    //We close the file and the process end
    input_file.close();
}

void stars::calculate_spectrum(const std::vector<double>& frequencies){
    //First we declare the vectors that we are going to use in the process
    //Spectrum is the vector that is going to store the calculated spectrum of each star
    std::vector<double> spectrum;
    //This vector represent the flux that we have read from each star
    std::vector<double> star_flux;
    double black_body_spectrum, spectrum_surface, star_radio,parsec_square_over_pi;
    //For each star...
    for (int i = 0; i < this -> number_of_stars; ++i) {
        //First we make space in the vector
        spectrum.resize(this -> number_of_frequencies);
        //We get the flux of the star
        star_flux = this -> stars_information[i].get_flux();
        //If the first value{ of the flux is negative, then blackbody is taken with a temperature equal to the absolute value of the flux
        if(star_flux[0] < 0){
            //Then we calculate the blackbody spectrum for each intensity
            for (int j = 0; j < this -> number_of_frequencies; ++j) {
                black_body_spectrum = common::black_body_planck_function(std::abs(star_flux[0]),frequencies[j]);
                spectrum[j] = black_body_spectrum;
            }
            //Then we assing the calculated spectrum to the star
            this -> stars_information[i].set_spectrum(spectrum);
        }
        //In other case, we need to calculate the intensity at the stellar surface.
        else{
            //We obtain the radio of the star
            star_radio = this -> stars_information[i].get_star_radio();
            for (int j = 0; j < this -> number_of_frequencies; ++j) {
                //Note that the spectrum is given in the usual flux-at-one-parsec format, so we have to translate this here to intensity at stellar surface
                //For this we are going to use:
                //F_v(r) = L_v / 4*PI*r^2; with L_v = 4*PI*R^2*PI*B_v(T_*)
                //TODO : Creo que esta es la formula, pero me quedan dudas
                //Note that parsec^2/pi = 3.0308410e32 m^2, but we want it in cm^2 = 3.0308410e36
                //TODO : Acorde a google parsec^2/PI = 3.0307577e32 (verificar con seba)
                //parsec_square_over_pi = 3.0307577e36;
                parsec_square_over_pi = 3.0308410e36;
                spectrum_surface = (parsec_square_over_pi * star_flux[j]) / std::pow(star_radio,2);
                spectrum[j] = spectrum_surface;
            }
            //Then we assign the calculated spectrum to the star
            this -> stars_information[i].set_spectrum(spectrum);
        }
        //We delete the calculated values in order to pass to the next star
        spectrum.clear();
    }
}

void stars::calculate_cumulative_spectrum(const std::vector<double>& mean_intensity){
    //This vector is the spectrum of each star
    std::vector<double> star_spectrum;
    //In this vector we are going to save the acumulated spectrum of each star
    //std::vector<double> cumulative_spectrum(this->number_of_frequencies+1);
    std::vector<double> cumulative_spectrum;
    //For each star...
    for (int i = 0; i < this -> number_of_stars; ++i) {
        cumulative_spectrum.resize(this -> number_of_frequencies+1);
        //star_cumulative = 0;
        cumulative_spectrum[0] = 0.0;
        //We obtain the calculated spectrum of each star
        star_spectrum = stars_information[i].get_spectrum();
        //Then for each spectrum, we calculate the cumulative spectrum until point j
        for (int j = 1; j < this -> number_of_frequencies+1; ++j) {
            cumulative_spectrum[j] = cumulative_spectrum[j-1] + star_spectrum[j-1] * mean_intensity[j-1];
        }
        for (int j = 0; j < this -> number_of_frequencies + 1; ++j) {
            cumulative_spectrum[j] = cumulative_spectrum[j] / cumulative_spectrum[this -> number_of_frequencies];
        }
        //Then we assign the cumulative spectrum to the star
        stars_information[i].set_cumulative_spectrum(cumulative_spectrum);
        //We delete the calculated values in order to pass to the next star
        cumulative_spectrum.clear();
    }
}

void stars::calculate_total_luminosities(const std::vector<double>& mean_intensity) {
    //TODO : Esta funcion requiere verificar caleta de cosas, por ahora hacemos lo básico
    double total_luminosity = 0.0;
    double star_luminosity;
    std::vector<double> star_spec;
    double star_radio;
    double x;
    double PI = 3.14159265358979323846264338328;
    double FOURPI = 12.5663706143591729538505735331;
    for (int i = 0; i < this -> number_of_stars; ++i) {
        star_spec = this -> stars_information[i].get_spectrum();
        star_radio = this -> stars_information[i].get_star_radio();
        star_luminosity = 0.0;
        std::cout.precision(17);
        for (int j = 0; j < this -> number_of_frequencies; ++j) {
            x = PI * star_spec[j] * FOURPI * std::pow(star_radio,2);
            star_luminosity = star_luminosity + x * mean_intensity[j];
        }
        this -> stars_information[i].set_luminosity(star_luminosity);
        total_luminosity = total_luminosity + star_luminosity;
    }
    this -> total_luminosity = total_luminosity;
    this -> cumulative_luminosity.resize(number_of_stars+1);
    this -> cumulative_luminosity[0] = 0.0;
    for (int i = 1; i < this->number_of_stars+1; ++i) {
        star_luminosity = this -> stars_information[i-1].get_luminosity();
        this -> cumulative_luminosity[i] = star_luminosity / total_luminosity;
    }
    this -> cumulative_luminosity[this->number_of_stars] = 1.0;
}

const std::vector<star>& stars::get_stars_information(void) const{
    return this -> stars_information;
}

stars::~stars(void){
    ;
}

void stars::calculate_energy(int number_of_photons) {
    //TODO : Cambiar acorde a entrada
    int number_of_sources = this -> number_of_stars;
    double star_energy;
    for (int i = 0; i < this -> number_of_stars; ++i) {
        if(this -> stars_information[i].get_luminosity() < 0.0){
            std::cerr << "ERROR: At star " << i << std::endl;
            std::cerr << "ERROR: Cannot treat star with zero luminosity" << std::endl;
            exit(0);
        }
    }
    for (int i = 0; i < this -> number_of_stars; ++i) {
        star_energy = this -> stars_information[i].get_luminosity() * number_of_sources / number_of_photons;
        this -> stars_information[i].set_energy(star_energy);
    }
}

void stars::jitter_stars(std::vector<double> cell_walls_x, std::vector<double> cell_walls_y, std::vector<double> cell_walls_z){
    double small_x,small_y,small_z,szx,szy,szz;
    int number_of_points_x = cell_walls_x.size();
    int number_of_points_y = cell_walls_y.size();
    int number_of_points_z = cell_walls_z.size();
    small_x = 0.48957949203816064943e-12;
    small_y = 0.38160649492048957394e-12;
    small_z = 0.64943484920957938160e-12;
    szx = cell_walls_x[number_of_points_x] - cell_walls_x[0];
    szy = cell_walls_y[number_of_points_y] - cell_walls_y[0];
    szz = cell_walls_z[number_of_points_z] - cell_walls_z[0];
    std::vector<double> star_position;
    for (int i = 0; i < this -> number_of_stars; ++i) {
        star_position = this -> stars_information[i].get_star_position();
        star_position[0] = star_position[0] + szx * small_x;
        star_position[1] = star_position[1] + szy * small_y;
        star_position[2] = star_position[2] + szz * small_z;
        this -> stars_information[i].set_star_position(star_position);
    }
}

void stars::fix_luminosities() {
    this -> cumulative_luminosity[0] = 0.0;
    if(this -> stars_information[0].get_luminosity() > 0.0){
        this -> cumulative_luminosity[1] = 1.0 / this -> number_of_stars;
    }
    else{
        this -> cumulative_luminosity[1] = 0.0;
    }
    for (int i = 1; i < this->number_of_stars; ++i) {
        if(this -> stars_information[i].get_luminosity() > 0.0){
            this -> cumulative_luminosity[i+1] = this -> cumulative_luminosity[i] + 1.0 / this -> number_of_stars;
        }
    }
    this -> cumulative_luminosity[this -> number_of_stars] = 1.0;
}

int stars::get_number_of_stars() const {
    return number_of_stars;
}

const std::vector<double> &stars::get_cumulative_luminosity() const {
    return cumulative_luminosity;
}


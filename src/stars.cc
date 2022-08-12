#include <stars.hh>

void stars::read_stars(void){
    //Read the input file for the stars information
    std::ifstream input_file;
    input_file.open("inputs/stars.inp");
    //In case that the file couldn't be opened, an error message is displayed.
    if(!input_file){
        std::cerr << "Mandatory file \"stars.inp\" could not be opened. Make sure that the file exists" << std::endl;
        exit(0);
    }
    int number_of_stars,iformat,number_of_lambdas;
    double star_radio,star_mass,x_position,y_position,z_position,lambda_point;
    //First we read the format of the file.
    //TODO : Ahora mismo siempre se trabaja con el valor 2. El manual indica que si el valor es 1 entonces
    //TODO : los valores de la frecuencia serán valores directamente en hertz y no en micrones
    //TODO : HAY QUE VERIFICAR SI ESTO AFECTA A NIVEL DE SIMULACION
    input_file >> iformat;
    this -> iformat = iformat;
    //Then we read the number of stars and the number of frequency points
    input_file >> number_of_stars;
    this -> number_of_stars = number_of_stars;
    //We resize the vector that is going to alocate each star
    this -> stars_information.resize(number_of_stars);
    input_file >> number_of_lambdas;
    this -> number_of_lambdas = number_of_lambdas;
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
        //TODO : Creo que sería bueno borrar la memoria de la estrella en el destructor
    }
    //TODO : En el codigo original verifican varias cosas sobre la grilla y la estrella.
    //TODO : Esto implica que hay que leer la grilla antes
    //Now we read the frequency points from the file
    //TODO : estas frecuencias son las mismas que ya se leyeron xd
    //TODO : el programa original compara si son las mismas o no y se detiene si no lo son
    //TODO : evaluar si es necesario, a priori se considera que no y no se trabaja con los valores
    for (int i = 0; i < number_of_lambdas; ++i) {
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
            star_flux.resize(number_of_lambdas);
            star_flux[0] = star_flux_value;
            for (int j = 1; j < number_of_lambdas; ++j) {
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

void stars::calculate_spectrum(std::vector<double> mean_intensity){
    //First we declare the vectors that we are going to use in the process
    //Spectrum is a the vector that is going to store the calculated spectrum of each star
    std::vector<double> spectrum;
    //This vector represent the flux that we have read from each star
    std::vector<double> star_flux;
    double black_body_spectrum, spectrum_surface, star_radio,parsec_square_over_pi;
    //For each star...
    for (int i = 0; i < this -> number_of_stars; ++i) {
        //First we make space in the vector
        spectrum.resize(this -> number_of_lambdas);
        //We get the flux of the star
        star_flux = this -> stars_information[i].get_flux();
        //If the first value{ of the flux is negative, then blackbody is taken with a temperature equal to the absolute value of the flux
        if(star_flux[0] < 0){
            //Then we calculate the blackbody spectrum for each intensity
            for (int j = 0; j < this -> number_of_lambdas; ++j) {
                black_body_spectrum = this -> black_body_planck_funcion(std::abs(star_flux[0]),mean_intensity[j]);
                spectrum[j] = black_body_spectrum;
            }
            //Then we assing the calculated spectrum to the star
            this -> stars_information[i].set_spectrum(spectrum);
        }
        //In other case, we need to calculate the intensity at the stellar surface.
        else{
            //We obtain the radio of the star
            star_radio = this -> stars_information[i].get_star_radio();
            for (int j = 0; j < this -> number_of_lambdas; ++j) {
                //Note that the spectrum is given in the usual flux-at-one-parsec format, so we have to translate this here to intensity at stellar surface
                //For this we are going to use:
                //F_v(r) = L_v / 4*PI*r^2; with L_v = 4*PI*R^2*PI*B_v(T_*)
                //TODO : Creo que esta es la formula, pero me quedan dudas
                //Note that parsec^2/pi = 3.0308410e32 m^2, but we want it in cm^2 = 3.0308410e36
                //TODO : Acorde a google parsec^2/PI = 3.0307577e32 (verificar con seba)
                parsec_square_over_pi = 3.0307577e36;
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

void stars::calculate_cumulative_spectrum(std::vector<double> mean_intensity){
    //This vector is the spectrum of each star
    std::vector<double> star_spectrum;
    //In this vector we are going to save the acumulated spectrum of each star
    std::vector<double> cumulative_spectrum(this->number_of_lambdas);
    double star_cumulative;
    //For each star...
    for (int i = 0; i < this -> number_of_stars; ++i) {
        star_cumulative = 0;
        //We make space in the vector to store the calculated cumulative spectrums
        cumulative_spectrum.resize(this -> number_of_lambdas);
        //We obtain the calculated spectrum of each star
        star_spectrum = stars_information[i].get_spectrum();
        //Then for each spectrum, we calculate the cumulative spectrum until point i
        for (int j = 0; j < this -> number_of_lambdas; ++j) {
            star_cumulative = star_cumulative + star_spectrum[j] * mean_intensity[j];
            cumulative_spectrum[i] = star_cumulative;
        }
        //The we assign the cumulative spectrum to the star
        stars_information[i].set_cumulative_spectrum(cumulative_spectrum);
        //We delete the calculated values in order to pass to the next star
        cumulative_spectrum.clear();
    }
}

double stars::black_body_planck_funcion(double temperature,double frequency){
    //This function calculates the Blackbody thermal radiation with the planck function
    //The planck function is definided as:
    //B_v(T) = [2 * h * v^3 / c^2] / [e^(h*v/K_b*T)-1]
    //Where h is the planck constant with value = 6.62607015e−34 J⋅Hz^−1, but we want it in erg⋅Hz^-1 = 6.6260701e-27 erg⋅Hz^-1
    //      v is the frequency in Hz
    //      c is the speed of light with value = 299792458 m / s^2, but we want it in cm / s^2 = 29979245800 cm / s^2
    //      K_b is the Boltzmann constant with value = 1.380649e−23 J⋅K^−1, but we want it in erg⋅K^-1 = 1.380649e-16 erg⋅K^-1
    double h = 6.6260701e-27;
    double c = 29979245800;
    double K_b = 380649e-16;
    return ((2*h*std::pow(frequency,3))/std::pow(c,2))/(std::exp(h*frequency/K_b*temperature)-1);
    //TODO : El codigo original usa otros valores para las constantes, evaluar con seba
}
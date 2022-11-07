#include <frequencies.hh>

frequencies::frequencies(void){
    ;
}

void frequencies::read_frequencies(void){
    //Light speed
    //double light_speed = 29979245800 // in cm / s^2
    //double light_speed = 2.9979245800000e10;
    //double light_speed = 2.9979245800000e10;
    double light_speed_microns = common::get_light_speed_microns();
    //Variable to store the number of frequency points for the stellar spectra
    int number_of_frequency_points;
    //Variable to store each frequency before transformation and after transformation
    double frequency, lambda;
    //We read the input file for the frequencies
    // TODO : Generalizar el método para leer waveleght o frequencies.inp
    // TODO : La unica diferencia radica en la transformación de los datos
    std::ifstream input_file;
    input_file.open("inputs/wavelength_micron.inp");
    //In case that the file couldn't be opened, an error message is displayed.
    if(!input_file){
        std::cerr << "Mandatory file \"wavelength_micron.inp\" could not be opened. Make sure that the file exists" << std::endl;
        exit(0);
    }
    //At first, we read the number of points
    input_file >> number_of_frequency_points;
    this -> number_of_frequencies = number_of_frequency_points;
    //We make space for each frequency point in the vectors
    this -> frequency_values.resize(number_of_frequency_points);
    this -> mean_intensity_values.resize(number_of_frequency_points);
    //Then we read each point
    for (int i = 0; i < number_of_frequency_points; ++i) {
        input_file >> lambda;
        //We are reading the wavelenght (lambda), but we want the frequency.
        //We can calculate it from f = v / lambda. Since the frequency is in micron, we need
        //to do C times 1e4 (that's going to be v in the equation) in order to get the lambda value.
        frequency = light_speed_microns / lambda;
        //TODO : En otra parte del código simplemente hacen que Hz = 2.99792458E+14 / wavelenght in micrones
        //TODO : Preguntar seba
        //frequency = 2.99792458e14 / lambda;
        this->frequency_values[i] = frequency;
    }
    // TODO : Solo funciona para funciones crecientes o decrecientes, por lo que hay que verificar que los puntos siempre suban o bajen
    //We close the file and the process end
    input_file.close();
}

void frequencies::calculate_mean_intensity(void){
    // TODO : The zeroth moment J ν is called the mean intensity and is indeed the angular average
    // TODO : of I ν (n). If we are in a homogeneous and isotropic radiation field, then J ν = I ν .
    // TODO : Afecta esto? o sea claro, afecta, pero no sé si tiene sentido para el problema a resolver
    //In order to calculate the mean intensity we need to calculate the next integral
    //1/2 ∫[-1,1] I_v(n,u)du
    // TODO : La integral es clara pero por qué I_v(n,u) du es la frecuencia???
    double mean_intensity;
    //We are going to iterate over each frequency in order to calculate the mean intensity
    for (int i = 0; i < this->number_of_frequencies; ++i) {
        //If we are in the first value, then we calculate the integral from [0,1]
        if(i == 0){
            mean_intensity = 0.5 * std::abs(this->frequency_values[1] - this->frequency_values[0]);
        }
        //If we are in the last value, then we calculate the integral from [-1,0]
        else if(i == this->number_of_frequencies-1){
            mean_intensity = 0.5 * std::abs(this->frequency_values[this->number_of_frequencies-1] - this->frequency_values[this->number_of_frequencies-2]);
        }
            //in any other case, we calculate the integral from [-1,1]
        else{
            mean_intensity = 0.5 * std::abs(this->frequency_values[i+1] - this->frequency_values[i-1]);
        }
        //We store the calculated mean intensity
        this -> mean_intensity_values[i] = mean_intensity;
    }
}

int frequencies::get_random_frequency(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
                                      const std::vector<double> &star_cumulative_spectrum) {
    //We generate a random number between 0 and 1
    double rn = uniform_zero_one_distribution(generator);
    //We search that value in the vector via hunt, and we get the index of the value in the vector
    int ray_inu = common::hunt(star_cumulative_spectrum, number_of_frequencies + 1, rn, number_of_frequencies);
    return ray_inu;
}

int frequencies::get_number_frequency_points() const{
    return this -> number_of_frequencies;
}

const std::vector<double>& frequencies::get_frequencies() const{
    return this -> frequency_values;
}

const std::vector<double>& frequencies::get_mean_intensity() const{
    return this -> mean_intensity_values;
}

frequencies::~frequencies(void){
    ;
}
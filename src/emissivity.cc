#include <emissivity.hh>

emissivity::emissivity(void){
    ;
}

void emissivity::generate_emissitivy_table(std::map<std::string,double> simulation_parameters,int number_of_species, int number_of_frequencies, std::vector<double> kappa_absorption, std::vector<double> freq_nu, std::vector<double> freq_dnu){
    //First we get the number of temperatures, temp0 and temp1 from the main input
    int ntemp = simulation_parameters["ntemp"];
    double temp0 = simulation_parameters["temp0"];
    double temp1 = simulation_parameters["temp1"];
    //Variable to store the precalculated temperature of each iteration
    double precalculated_temperature;
    //Then we create a vector to contain the precalculated temperatures
    std::vector<double> db_temp(ntemp);
    double x,y;
    for (int i = 0; i < ntemp; ++i) {
        //TODO : Preguntar formula de esto
        x = temp1/temp0;
        y = i/(ntemp-1.0);
        precalculated_temperature = temp0 * std::pow(x,y);
        db_temp[i] = precalculated_temperature;
    }

    double demis;
    std::vector<double> fnu_diff;
    //This vector its going to store the precaculated temperatures
    this -> db_enertemp.resize(number_of_species, std::vector<double> (ntemp, 0));
    //This vector its going to store the logaritmic temperature
    this -> db_logenertemp.resize(number_of_species, std::vector<double> (ntemp,0));
    //This vector its going to store the emissivity of each specie at each temperature
    this -> db_emiss.resize(number_of_species, std::vector<std::vector<double>> (ntemp,std::vector<double>(ntemp,0)));
    //We iterate over each dust specie
    for (int i = 0; i < number_of_species; ++i) {
        //And over each temperature
        for (int j = 0; j < ntemp; ++j) {
            //TODO : Preguntar por el evento de absorcion
            fnu_diff = absoprtion_event(db_temp[j],number_of_frequencies, kappa_absorption, freq_nu, freq_dnu);
            //We stored the last value in the fnu_diff vector, so now we get it
            demis = fnu_diff.back();
            fnu_diff.pop_back();
            //We store the relevant values in each vector
            db_enertemp[i][j] = demis;
            db_logenertemp[i][j] = std::log(demis);
            db_emiss[i][j] = fnu_diff;
        }
    }
}

//TODO : Verificar vectores que se usan en la simulacion de los calculados para almacenar en el objeto de la clase
void emissivity::compute_derivate(int number_of_species, int number_of_temperatures, int number_of_frequencies, std::vector<double> freq_dnu) {
    //This vector its going to represent the emisitivity of each specie at each temperature
    std::vector<double> diffemis(number_of_frequencies);
    //This vector its going to store the cumulative value of the calculated derivative
    std::vector<double> db_cumul(number_of_frequencies + 1);
    //In this cube, we are going to store the normal cumulative value of the derivative for each specie at each temperature
    std::vector<std::vector<std::vector<double>>> db_cumulnorm(number_of_species,std::vector<std::vector<double>>(number_of_temperatures,std::vector<double>(number_of_frequencies+1,0)));
    //This vector its going to store the cumulative normal value for a specific specie at a specific temperature
    std::vector<double> db_cumulnorm_values(number_of_frequencies+1);
    //We iterate over each specie
    for (int i = 0; i < number_of_species; ++i) {
        //For the first temeprature, we simply the emissivity function itself
        diffemis = this -> db_emiss[i][0];
        //Now we accumulate over each frequency
        db_cumul[0] = 0.0;
        for (int j = 1; j < number_of_frequencies + 1; ++j) {
            db_cumul[j] = db_cumul[j-1]  + diffemis[j] * freq_dnu[j];
        }
        //Now we obtain the normal value of the calculated accumulation
        for (int j = 0; j < number_of_frequencies + 1; ++j) {
            db_cumulnorm_values[j] = db_cumul[j] / db_cumul[number_of_frequencies];
            //db_cumulnorm[i][0][j] = db_cumul[j] / db_cumul[number_of_frequencies];
        }
        //We store the normal cumulative value
        db_cumulnorm[i][0] = db_cumulnorm_values;
        //For the rest of the temperatures, we take the difference
        for (int j = 1; j < number_of_temperatures; ++j) {
            for (int k = 0; k < number_of_frequencies + 1; ++k) {
                //We take the diference of the emissivity of specie i, at temperature j, at frequency k
                diffemis[k] = this -> db_emiss[i][j][k] - this -> db_emiss[i][j-1][k];
                //For the fist frequency we store 0
                if(k == 0){
                    db_cumul[0] = 0.0;
                }
                //For any other, we accumulate with the last calculated value
                else{
                    db_cumul[k] = db_cumul[k-1] + diffemis[k-1] * freq_dnu[k-1];
                }
            }
            //Now we obtain the normal value of the calculated accumulation
            for (int k = 0; k < number_of_frequencies + 1; ++k) {
                db_cumulnorm_values[k] = db_cumul[k] /db_cumul[number_of_frequencies];
            }
            //We store the normal cumulative value
            db_cumulnorm[i][j] = db_cumulnorm_values;
        }
    }
}

std::vector<double> emissivity::absoprtion_event(double temp, int number_of_frequencies, std::vector<double> cell_alpha, std::vector<double> freq_nu, std::vector<double> freq_dnu){
    std::vector<double> fnu_diff(number_of_frequencies+1);
    double fourpi = 12.5663706143591729538505735331;
    double demis = 0.0;
    for (int i = 0; i < number_of_frequencies; ++i) {
        fnu_diff[i] = cell_alpha[i] * common::black_body_planck_function(temp,freq_nu[i]);
        demis = demis + (fnu_diff[i] * freq_dnu[i]);
    }
    demis = demis * fourpi;
    fnu_diff[number_of_frequencies] = demis;
    return fnu_diff;
}

emissivity::~emissivity(void){
    ;
}
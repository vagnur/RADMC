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
    std::vector<std::vector<double>> db_enertemp(number_of_species, std::vector<double> (ntemp, 0));
    std::vector<std::vector<std::vector<double>>> db_emiss(number_of_species, std::vector<std::vector<double>> (ntemp,std::vector<double>(ntemp,0)));
    std::vector<std::vector<double>> db_logenertemp(number_of_species, std::vector<double> (ntemp,0));
    //We iterate over each dust specie
    for (int i = 0; i < number_of_species; ++i) {
        //And over each temperature
        for (int j = 0; j < ntemp; ++j) {
            //cell_alpha = kappa_absorption
            fnu_diff = absoprtion_event(db_temp[j],number_of_frequencies, kappa_absorption, freq_nu, freq_dnu);
            demis = fnu_diff.back();
            fnu_diff.pop_back();
            db_enertemp[i][j] = demis;
            db_logenertemp[i][j] = std::log(demis);
            db_emiss[i][j] = fnu_diff;
        }
    }
}

std::vector<double> emissivity::absoprtion_event(double temp, int number_of_frequencies, std::vector<double> cell_alpha, std::vector<double> freq_nu, std::vector<double> freq_dnu){
    std::vector<double> fnu_diff(number_of_frequencies+1);
    double fourpi = 12.5663706143591729538505735331;
    double demis = 0.0;
    for (int i = 0; i < number_of_frequencies; ++i) {
        fnu_diff[i] = cell_alpha[i] * common::black_body_planck_funcion(temp,freq_nu[i]);
        demis = demis + (fnu_diff[i] * freq_dnu[i]);
    }
    demis = demis * fourpi;
    fnu_diff[number_of_frequencies] = demis;
    return fnu_diff;
}

emissivity::~emissivity(void){
    ;
}
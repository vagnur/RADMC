#include <emissivity.hh>

emissivity::emissivity(void){
    ;
}

void emissivity::generate_emissivity_table(std::map<std::string,double>& simulation_parameters, int number_of_species, const std::vector<dust_species>& dust_species_information, const std::vector<double>& freq_nu, const std::vector<double>& freq_dnu){
    int number_of_frequencies = freq_nu.size();
    //First we get the number of temperatures, temp0 and temp1 from the main input
    int ntemp = simulation_parameters["ntemp"];
    double temp0 = simulation_parameters["temp0"];
    double temp1 = simulation_parameters["temp1"];
    this -> number_of_temperatures = ntemp;
    this -> temp0 = temp0;
    this -> temp1 = temp1;
    //Variable to store the precalculated temperature of each iteration
    double precalculated_temperature;
    //Then we create a vector to contain the precalculated temperatures
    this -> db_temp.resize(ntemp);
    //std::vector<double> db_temp(ntemp);
    double x,y;
    for (int i = 0; i < ntemp; ++i) {
        //TODO : Preguntar formula de esto
        x = temp1/temp0;
        y = i/(ntemp-1.0);
        precalculated_temperature = temp0 * std::pow(x,y);
        this -> db_temp[i] = precalculated_temperature;
    }
    double demis;
    std::vector<double> fnu_diff;
    //This vector is going to store the precalculated temperatures for each specie
    this -> db_enertemp.resize(number_of_species, std::vector<double> (ntemp, 0));
    //This vector is going to store the logarithmic temperature
    this -> db_logenertemp.resize(number_of_species, std::vector<double> (ntemp,0));
    //This vector is going to store the emissivity of each specie, at each temperature and at each frequency
    this -> db_emiss.resize(number_of_species, std::vector<std::vector<double>> (ntemp,std::vector<double>(number_of_frequencies,0)));
    //We iterate over each dust specie
    for (int i = 0; i < number_of_species; ++i) {
        //And over each temperature
        for (int j = 0; j < ntemp; ++j) {
            //TODO : Preguntar por este evento
            fnu_diff = absorption_event(this -> db_temp[j], number_of_frequencies, dust_species_information[i].get_kappa_absorption_interpoled(), freq_nu, freq_dnu);
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
void emissivity::compute_derivate(int number_of_species, int number_of_temperatures, const std::vector<double>& freq_dnu) {
    int number_of_frequencies = freq_dnu.size();
    //This vector is going to represent the emissivity of each specie at each temperature
    std::vector<double> diffemis(number_of_frequencies);
    //This vector is going to store the cumulative value of the calculated derivative
    std::vector<double> db_cumul(number_of_frequencies + 1);
    //In this cube, we are going to store the normal cumulative value of the derivative for each specie, at each temperature and at each frequency
    this -> db_cumulnorm.resize(number_of_species,std::vector<std::vector<double>>(number_of_temperatures,std::vector<double>(number_of_frequencies+1,0)));
    //This vector is going to store the cumulative normal value for a specific specie at a specific temperature
    std::vector<double> db_cumulnorm_values(number_of_frequencies+1);
    //We iterate over each specie
    for (int i = 0; i < number_of_species; ++i) {
        //For the first temperature, we take simply the emissivity function itself
        diffemis = this -> db_emiss[i][0];
        //Now we accumulate over each frequency
        db_cumul[0] = 0.0;
        for (int j = 1; j < number_of_frequencies + 1; ++j) {
            db_cumul[j] = db_cumul[j-1]  + diffemis[j-1] * freq_dnu[j-1];
        }
        //Now we obtain the normal value of the calculated accumulation
        for (int j = 0; j < number_of_frequencies + 1; ++j) {
            db_cumulnorm_values[j] = db_cumul[j] / db_cumul[number_of_frequencies];
            //db_cumulnorm[i][0][j] = db_cumul[j] / db_cumul[number_of_frequencies];
        }
        //We store the normal cumulative value
        this -> db_cumulnorm[i][0] = db_cumulnorm_values;
        //For the rest of the temperatures, we take the difference
        for (int j = 1; j < number_of_temperatures; ++j) {
            for (int k = 0; k < number_of_frequencies + 1; ++k) {
                //We take the difference of the emissivity of specie i, at temperature j, at frequency k
                diffemis[k] = this -> db_emiss[i][j][k] - this -> db_emiss[i][j-1][k];
                //For the fist frequency we store 0 ad db_cumul
                if(k == 0){
                    db_cumul[0] = 0.0;
                }
                //For any other case, we accumulate with the last calculated value
                else{
                    db_cumul[k] = db_cumul[k-1] + diffemis[k-1] * freq_dnu[k-1];
                }
            }
            //Now we obtain the normal value of the calculated accumulation
            for (int k = 0; k < number_of_frequencies + 1; ++k) {
                db_cumulnorm_values[k] = db_cumul[k] /db_cumul[number_of_frequencies];
            }
            //We store the normal cumulative value
            this -> db_cumulnorm[i][j] = db_cumulnorm_values;
        }
    }
}

std::vector<double> emissivity::absorption_event(double temp, int number_of_frequencies, const std::vector<double>& cell_alpha, const std::vector<double>& freq_nu, const std::vector<double>& freq_dnu){
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

double emissivity::compute_dust_temp_energy(double energy, int iSpec){
    double tempReturn = 0.0;
    int itemp = 0;
    double logEner = std::log(energy);
    double eps;
    //TODO : Obtener dblogenergtemp de la especie
    itemp = common::hunt(this -> db_logenertemp[iSpec], this -> number_of_temperatures, logEner, this -> number_of_temperatures);
    //itemp = common::hunt(dbLogEnerTemp[iSpec], number_of_temperatures, logEner, number_of_temperatures);
    //printf("itemp=%d\n", *itemp);
    if (itemp >= this -> number_of_temperatures - 1) {
        std::cout << "aca 1" << std::endl;
        std::cerr << "ERROR : Too high temperature discovered" << std::endl;
        exit(0);
    }

    if (itemp <= -1) {
        //Temperature presumably below lowest temp in dbase
        //TODO : Obtener enertemp de la especie
        //eps = energy / dbEnerTemp[iSpec][0];
        eps = energy / db_enertemp[iSpec][0];
        if (eps >= 1) {
            //printf("exit\n");
            //exit(1);
            std::cout << "aca 2" << std::endl;
            std::cerr << "ERROR : Too high temperature discovered" << std::endl;
            exit(0);
        }
        //tempReturn = eps * dbTemp[0];
        tempReturn = eps * this -> db_temp[0];
    } else {
        //TODO : Obtener enertemp de la especie
        //eps = (energy - dbEnerTemp[iSpec][itemp]) / (dbEnerTemp[iSpec][itemp + 1] - dbEnerTemp[iSpec][itemp]);
        eps = (energy - this -> db_enertemp[iSpec][itemp]) / (this -> db_enertemp[iSpec][itemp + 1] - this -> db_enertemp[iSpec][itemp]);
        if ((eps > 1) || (eps < 0)) {
            std::cerr << "ERROR : Temperature found out of range..." << std::endl;
            //printf("exit\n");
            exit(0);
        }
        //TODO : Obtener dbTemp
        //tempReturn = (1.0 - eps) * dbTemp[itemp] + eps * dbTemp[itemp + 1];
        tempReturn = (1.0 - eps) * this -> db_temp[itemp] + eps * db_temp[itemp + 1];
    }
    return tempReturn;
}

emissivity::~emissivity(void){
    ;
}

const std::vector<double>& emissivity::get_db_temp() const{
    return this -> db_temp;
}

const std::vector<std::vector<double>>& emissivity::get_db_enertemp() const{
    return this -> db_enertemp;
}

const std::vector<std::vector<double>>& emissivity::get_db_logenertemp() const{
    return this -> db_logenertemp;
}

const std::vector<std::vector<std::vector<double>>>& emissivity::get_db_cumulnorm() const{
    return this -> db_cumulnorm;
}
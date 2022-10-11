#ifndef RADMC_EMISSIVITY_HH
#define RADMC_EMISSIVITY_HH

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <common.hh>
#include <dust.hh>

class emissivity {

private:

    std::vector<double> db_temp;
    std::vector<std::vector<double>> db_enertemp;
    std::vector<std::vector<double>> db_logenertemp;
    std::vector<std::vector<std::vector<double>>> db_emiss;
    std::vector<std::vector<std::vector<double>>> db_cumulnorm;

public:

    //Empty constructor
    emissivity(void);

    //
    void generate_emissivity_table(std::map<std::string,double>& simulation_parameters, int number_of_species, const std::vector<dust_species>& dust_species_information, const std::vector<double>& freq_nu, const std::vector<double>& freq_dnu);

    //
    void compute_derivate(int number_of_species, int number_of_temperatures, const std::vector<double>& freq_dnu);

    //
    std::vector<double> absorption_event(double temp, int number_of_frequencies, const std::vector<double>& cell_alpha, const std::vector<double>& freq_nu, const std::vector<double>& freq_dnu);

    const std::vector<double> &get_db_temp() const;

    const std::vector<std::vector<double>> & get_db_enertemp() const;

    const std::vector<std::vector<double>> &get_db_logenertemp() const;


    //Empty destructor
    ~emissivity(void);

    const std::vector<std::vector<std::vector<double>>> &get_db_cumulnorm() const;
};

#endif //RADMC_EMISSIVITY_HH

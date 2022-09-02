#ifndef RADMC_EMISSIVITY_HH
#define RADMC_EMISSIVITY_HH

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <common.hh>

class emissivity {

private:



public:

    //Empty constructor
    emissivity(void);

    void generate_emissitivy_table(std::map<std::string,double> simulation_parameters, int number_of_species, int number_of_frequencies, std::vector<double> kappa_absorption, std::vector<double> freq_nu, std::vector<double> freq_dnu);

    std::vector<double> absoprtion_event(double temp, int number_of_frequencies, std::vector<double> cell_alpha, std::vector<double> freq_nu, std::vector<double> freq_dnu);

    //Empty destructor
    ~emissivity(void);

};


#endif //RADMC_EMISSIVITY_HH
#ifndef RADMC_EMISSIVITY_HH
#define RADMC_EMISSIVITY_HH

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <common.hh>
#include <dust.hh>


/*
 * This class is going to store the precalculated temperatures of each specie.
 * The precalculated temperatures are used to save computation time.
 */

class emissivity {

private:

    int number_of_temperatures;
    double temp0;
    double temp1;
    //This vector stores the precalculated temperatures
    std::vector<double> db_temp;
    //This vector is going to store the precalculated temperatures of each specie
    std::vector<std::vector<double>> db_enertemp;
    //Same as above, but in logarithmic domain
    std::vector<std::vector<double>> db_logenertemp;
    //This vector is going to store the emissivity of each specie, at each temperature and at each frequency
    std::vector<std::vector<std::vector<double>>> db_emiss;
    //This vector is going to store the normal cumulative value of the derivative for each specie, at each temperature and at each frequency
    std::vector<std::vector<std::vector<double>>> db_cumulnorm;

public:

    //Empty constructor
    emissivity(void);

    //Function to generate the precalculated temperatures of the simulation
    //Input: map simulation_parameters -> this map contains the simulation parameters, for this function, we use the number of temperatures and the limits of the precalculated temperatures
    //       int number_of_species -> number of dust species
    //       vector dist_species_information -> vector with the information of each dust specie
    //       vector freq_nu -> Vector with the frequencies
    //       vector freq_dnu -> Vector with the mean frequencies
    //Output: It has no output, but each class vector is going to store the relevant information
    void generate_emissivity_table(std::map<std::string,double>& simulation_parameters, int number_of_species, const std::vector<dust_species>& dust_species_information, const std::vector<double>& freq_nu, const std::vector<double>& freq_dnu);

    //This function calculates the derivative of the temperatures in order to calculate the cumulative value of each specie at each temperature
    //Input: int number_of_species -> number of dust species
    //       int number_of_temperatures -> number of temperatures
    //       vector freq_dnu -> Vector with the mean frequencies
    void compute_derivate(int number_of_species, int number_of_temperatures, const std::vector<double>& freq_dnu);

    //This function generate the absorption event of each frequency. This function is called for each specie at each temperature
    //Input: double temp -> The temperature that is going to be used for the black body planck function
    //       int number_of_frequencies -> number of frequencies
    //       vector cell_alpha -> Value of the absorption opacity of the present cell of the grid
    //       vector freq_nu -> Vector with the frequencies values
    //       vector freq_dnu -> Vector with the mean frequencies
    //Output: Vector with the cumulative emissivity and the emissivity of each frequency
    std::vector<double> absorption_event(double temp, int number_of_frequencies, const std::vector<double>& cell_alpha, const std::vector<double>& freq_nu, const std::vector<double>& freq_dnu);


    //Getters
    const std::vector<double> &get_db_temp() const;

    const std::vector<std::vector<double>> & get_db_enertemp() const;

    const std::vector<std::vector<double>> &get_db_logenertemp() const;

    const std::vector<std::vector<std::vector<double>>> &get_db_cumulnorm() const;

    //Empty destructor
    ~emissivity(void);

    double compute_dust_temp_energy(double energy, int iSpec);
};

#endif //RADMC_EMISSIVITY_HH

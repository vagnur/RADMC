#ifndef RADMC_MONTE_CARLO_HH
#define RADMC_MONTE_CARLO_HH

#include <map>
#include <vector>
#include <random>
#include <iostream>

#include <common.hh>
#include <grid.hh>
#include <frequencies.hh>
#include <stars.hh>
#include <dust.hh>
#include <emissivity.hh>
#include <photon.hh>

class monte_carlo {

private:

    //Grid object
    grid m_grid;
    //Frequencies object
    frequencies m_frequencies;
    //Stars object
    stars m_stars;
    //Dust object
    dust m_dust;
    //Emissivity object
    emissivity m_emissivity;
    //Vector with the pothons that we are going to use in the simulation
    std::vector<photon> m_photons;

public:

    //Empty constuctor
    monte_carlo(void);

    //Method to read the main file of the simulator
    std::map<std::string,double> read_main_file(void);

    //Method that execute the monte carlo therm in a regular cartesian grid
    void do_monte_carlo_therm_regular_cartesian_grid(void);

    //Method that initialize the photons that we are going to launch in the grid
    //Input : generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    //        int number_of_photons -> Number of photons that we want to initialize
    //Output : It has no output, but the photons vector is going to store a new initialized vector in each index position
    void initialize_cartesian_regular_photons(std::mt19937 &generator,
                                              std::uniform_real_distribution<> &uniform_zero_one_distribution,
                                              int number_of_photons);

    //This method move the selected photon to the next position and accumulate energy in the dust species of the grid position
    //Input : photon photon_i -> Is the <i> photon that we are going to use in the method
    //        generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    //Output -> It has no output, but the pothon is going to a new grid position and the dust species absorve energy in the process
    void move(photon &photon_i, std::mt19937 &generator,
              std::uniform_real_distribution<> &uniform_zero_one_distribution);

    //This method add temperature to each dust specie in the photon grid position
    //Input : photon photon_i -> Is the <i> photon that we are going to use in the method
    //Output : It has no output, but is going to store temperature in each dust specie vector
    void add_temperature_decoupled(photon& photon_i);

    //This method divide the energy in the dust species accord to the dust density and the absorption properties of the dust
    //Input : photon photon_i -> Is the <i> photon that we are going to use in the method
    //Output : It has no output, but the vector is going to store the energy that give to each dust specie
    void divide_absorved_energy(photon& photon_i);

    //TODO : Evaluar dejar esto en emisividad
    int
    pickRandomFreqDb(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_real_distribution,
                     photon& photon_i,
                     int number_of_temperatures);

    //
    void do_absorption_event(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_real_distribution,
                             photon& photon_i, int number_of_temperatures);

    void launch_photons(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_zero_one_distribution,
                        int scattering_mode, int number_of_photons, int number_of_temperatures);

    void calculate_dust_temperature();

    void write_dust_temperatures_file();

    void get_random_direction(photon &photon_i, std::mt19937 &generator,
                              std::uniform_real_distribution<> &uniform_zero_one_distribution);

    void get_henvey_greenstein_direction(photon &photon_i, std::mt19937 &generator,
                                         std::uniform_real_distribution<> &uniform_zero_one_distribution,
                                         const std::vector<double> &g);

    ~monte_carlo(void);
};


#endif //RADMC_MONTE_CARLO_HH

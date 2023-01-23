#ifndef RADMC_MONTE_CARLO_HH
#define RADMC_MONTE_CARLO_HH

#include <map>
#include <vector>
#include <random>
#include <iostream>

#include <common.hh>
//#include <grid.hh>
#include <cartesian_regular_grid.hh>
#include <spherical_regular_grid.hh>
#include <frequencies.hh>
#include <stars.hh>
#include <dust.hh>
#include <emissivity.hh>
#include <photon.hh>

/*
 * An object of this class is going to initialize the necessary objects to start a simulation.
 * Is going to simulate each photon on the grid and is goint to generate an output file with the temperature of the dust
 *  in each grid position.
 */

class monte_carlo {

private:

    //Grid object
    general_grid *m_grid;
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

    //This method move_photon the selected photon to the next position and accumulate energy in the dust species of the grid position
    //Input : photon photon_i -> Is the <i> photon that we are going to use in the method
    //        generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    //Output -> It has no output, but the pothon is going to a new grid position and the dust species absorve energy in the process
    void move_photon(photon &photon_i, std::mt19937 &generator,
                     std::uniform_real_distribution<> &uniform_zero_one_distribution);

    //This method add temperature to each dust specie in the photon grid position
    //Input : photon photon_i -> Is the <i> photon that we are going to use in the method
    //Output : It has no output, but is going to store temperature in each dust specie vector
    void add_temperature_decoupled(photon& photon_i);

    //This method divide the energy in the dust species accord to the dust density and the absorption properties of the dust
    //Input : photon photon_i -> Is the <i> photon that we are going to use in the method
    //Output : It has no output, but the vector is going to store the energy that give to each dust specie
    void divide_absorbed_energy(photon& photon_i);

    //This method give a random frequency
    //Input : photon photon_i -> Is the <i> photon that we are going to use in the method
    //        generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    //        int number of temperatures -> number of temperatures in the temperatures DB
    //Output : The index of the new frequency
    int pickRandomFreqDb(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_real_distribution,
                     photon& photon_i, int number_of_temperatures);

    //This method simulate the absorption event between the photon and the dust were the photon is
    //Input : photon photon_i -> Is the <i> photon that we are going to use in the method
    //        generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    //        int number of temperatures -> number of temperatures in the temperatures DB
    //Output : It has no output, but the photon and the dust species are going to store new information
    void do_absorption_event(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_real_distribution,
                             photon& photon_i, int number_of_temperatures);

    //This method launch each photon on the grid, and supervise the movement of the photon until it leave the grid
    //Input : generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    //        int scattering mode -> Type of scattering that we are going to use in the simulation
    //        int number of temperatures -> number of temperatures in the temperatures DB
    //        int number of photons -> number of photons to use in the simulation
    //Output : It has no output
    void launch_photons(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_zero_one_distribution,
                        int scattering_mode, int number_of_photons, int number_of_temperatures);

    //Method to generate a random direction for the photon
    //Input : photon photon_i -> Is the <i> photon that we are going to use in the method
    //        generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    //Output : It has no output, but the photon is going to store a new direction
    void get_random_direction(photon &photon_i, std::mt19937 &generator,
                              std::uniform_real_distribution<> &uniform_zero_one_distribution);

    //Method to obtain a henvey greenstein direction for the scattering mode = 2.
    //Input : photon photon_i -> Is the <i> photon that we are going to use in the method
    //        generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    //        vector g : vector with the g values of the dust
    //Output : It has no output, but the photon is going to store a new direction
    void get_henvey_greenstein_direction(photon &photon_i, std::mt19937 &generator,
                                         std::uniform_real_distribution<> &uniform_zero_one_distribution,
                                         const std::vector<double> &g);

    //This method calculate the dust temperature of each dust specie in each cell of the grid
    //Input : It has no input
    //Output : It has no output, but each dust specie is going to store it temperature in each position of the grid
    void calculate_dust_temperature();

    //This method is going to write the dust temperature file
    //Input : It has no input
    //Output : It has no output, but each dust specie is going to store it temperature in each position of the grid
    void write_dust_temperatures_file();

    ~monte_carlo(void);
};


#endif //RADMC_MONTE_CARLO_HH

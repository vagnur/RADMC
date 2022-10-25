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

    grid grid_object;
    frequencies frequencies_object;
    stars stars_object;
    dust dust_object;
    emissivity emissivity_object;
    std::vector<photon> photons;

public:

    monte_carlo(void);

    void do_monte_carlo_therm_regular_cartesian_grid();

    std::map<std::string,double> read_main_file(void);

    ~monte_carlo(void);

    void do_monte_carlo_therm_regular_cartesian_grid_2();

    void intialize_cartesian_regular_photons(std::mt19937 &generator,
                                             std::uniform_real_distribution<> &uniform_zero_one_distribution,
                                             int number_of_photons);

    void walk_next_event(photon &photon_i, std::mt19937 &generator,
                         std::uniform_real_distribution<> &uniform_zero_one_distribution);

    void add_temperature_decoupled(photon& photon_i);

    void divide_absorved_energy(photon& photon_i);

    int
    pickRandomFreqDb(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_real_distribution,
                     photon& photon_i,
                     int number_of_temperatures);

    void do_absorption_event(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_real_distribution,
                             photon& photon_i, int number_of_temperatures);

    void launch_photons(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_zero_one_distribution,
                        int scattering_mode, int number_of_photons, int number_of_temperatures);

    void calculate_dust_temperature();

    void write_temperatures_file();
};


#endif //RADMC_MONTE_CARLO_HH

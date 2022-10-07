#ifndef RADMC_MONTE_CARLO_HH
#define RADMC_MONTE_CARLO_HH

#include <map>
#include <vector>
#include <random>

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

};


#endif //RADMC_MONTE_CARLO_HH

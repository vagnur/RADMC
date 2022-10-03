#ifndef RADMC_MONTECARLO_HH
#define RADMC_MONTECARLO_HH

#include <map>

#include <common.hh>
#include <grid.hh>
#include <frequencies.hh>
#include <stars.hh>
#include <dust.hh>

class monteCarlo {

private:

    grid grid_object;
    frequencies frequencies_object;
    stars stars_object;
    dust dust_object;


public:

    monteCarlo(void);

    void do_monte_carlo_therm_regular_cartesian_grid();

    std::map<std::string,double> read_main_file(void);

    ~monteCarlo(void);

};


#endif //RADMC_MONTECARLO_HH

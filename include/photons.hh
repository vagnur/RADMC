#ifndef RADMC_PHOTONS_HH
#define RADMC_PHOTONS_HH

#include <vector>
#include <photon.hh>

class photons {

private:

    int number_of_photons;
    std::vector<std::vector<std::vector<std::vector<double>>>> cumulative_energy_specie;
    std::vector<std::vector<std::vector<std::vector<double>>>> temperatures;
    std::vector<photon> pothon_information;

public:

    photons(void);

    ~photons();

};


#endif
#ifndef RADMC_PHOTONS_HH
#define RADMC_PHOTONS_HH

#include <vector>
#include <photon.hh>

class photons {

private:

    std::vector<photon> pothon_information;

public:

    photons(void);

    ~photons();

};


#endif
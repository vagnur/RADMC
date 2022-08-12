#ifndef RADMC_DUST_SPECIES_HH
#define RADMC_DUST_SPECIES_HH

#include <vector>
#include <string>

/*
 * An object of this class represent a dust specie present in the simulation
 * Is going to store the relevant information of each specie, as the name, the ... and the ...
 */
class dust_species {

private:

    //A 3D vector that store the densities of the dust in the grid
    std::vector<std::vector<std::vector<double>>> densities;
    //Metadata from dustocap.inp file
    //This parameter is 1 if we are going to read "dustkappa_*.inp" file.
    //If 10, then we read "dustkapscatmat_*.inp" file
    int input_style;
    //For now, this is always 0.
    int iquantum;
    //Name of the dust specie
    std::string specie_name;
    //Vector to store the wavelength points
    std::vector<double> lambda;
    //Vector to store the absorption opacity
    std::vector<double> kappa_absorption;
    //Vector to store the scattering opacity
    std::vector<double> kappa_scattering;
    //Vector to store the mean scattering angle
    std::vector<double> g;
    //Vector to store the frequency points
    std::vector<double> frequency;

public:

    //Empty constructor
    dust_species(void);

    //Setter for the density
    void set_density(std::vector<std::vector<std::vector<double>>> densities);

    //Setter for the wavelength points
    void set_lambda(std::vector<double> lambda);

    //Setter for the absorption opacity
    void set_kappa_absorption(std::vector<double> kappa_absorption);

    //Setter for the scattering opacity
    void set_kappa_scattering(std::vector<double> kappa_scattering);

    //Setter for the scattering angle
    void set_g(std::vector<double> g);

    //Setter for the frequency points
    void set_frequency(std::vector<double> frequency);

    //Empty destructor
    ~dust_species(void);

};

#endif //RADMC_DUST_SPECIES_HH

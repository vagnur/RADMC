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

    //A 3D vector that store the densities of the dust specie in the grid
    std::vector<std::vector<std::vector<double>>> densities;
    //A 3D vector that store the temperature of the dust specie in the grid
    std::vector<std::vector<std::vector<double>>> temperatures;
    //A 3D vector hat store the cumulative energy of the dust specie in the grid
    std::vector<std::vector<std::vector<double>>> cumulative_energy;
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
    //Vector to store the absorption opacity (remapped)
    std::vector<double> kappa_absorption_remapped;
    //Vector to store the scattering opacity (remapped)
    std::vector<double> kappa_scattering_remapped;
    //Vector to store the mean scattering angle (remapped)
    std::vector<double> g_remapped;
    //Vector to store the absorption opacity (interpoled)
    std::vector<double> kappa_absorption_interpoled;
    //Vector to store the scattering opacity (remapped)
    std::vector<double> kappa_scattering_interpoled;
    //Vector to store the mean scattering angle (remapped)
    std::vector<double> g_interpoled;
    //Vector to store the frequency points
    std::vector<double> frequency;

public:

    //Empty constructor
    dust_species(void);

    void initialize_temperature(int number_of_points_x,int number_of_points_y,int number_of_points_z);

    //Setter for the density
    void set_density(const std::vector<std::vector<std::vector<double>>>& densities);

    //Setter for the wavelength points
    void set_lambda(const std::vector<double>& lambda);

    //Setter for the absorption opacity
    void set_kappa_absorption(const std::vector<double>& kappa_absorption);

    //Setter for the scattering opacity
    void set_kappa_scattering(const std::vector<double>& kappa_scattering);

    //Setter for the scattering angle
    void set_g(const std::vector<double>& g);

    //Setter for the frequency points
    void set_frequency(const std::vector<double>& frequency);

    //Getter for the densities
    const std::vector<std::vector<std::vector<double>>>& get_densities() const;

    //Getter for frequency
    const std::vector<double>& get_frequency() const;

    //Getter for absorption
    const std::vector<double>& get_absoprtion() const;

    //Getter for scattering
    const std::vector<double>& get_scattering() const;

    //Getter for mean scattering angle
    const std::vector<double>& get_g() const;

    const std::vector<double>& get_lambda() const;

    const std::vector<double>& get_kappa_absorption_remapped() const;

    void set_kappa_absorption_remapped(const std::vector<double>& kappa_absorption_remapped);

    const std::vector<double>& get_kappa_scattering_remapped() const;

    void set_kappa_scattering_remapped(const std::vector<double>& kappa_scattering_remapped);

    const std::vector<double>& get_g_remapped() const;

    void set_g_remapped(const std::vector<double>& g_remapped);

    const std::vector<double>& get_kappa_absorption_interpoled() const;

    void set_kappa_absorption_interpoled(const std::vector<double>& kappa_absorption_interpoled);

    const std::vector<double>& get_kappa_scattering_interpoled() const;

    void set_kappa_scattering_interpoled(const std::vector<double>& kappa_scattering_interpoled);

    const std::vector<double>& get_g_interpoled() const;

    void set_g_interpoled(const std::vector<double>& g_interpoled);

    void add_energy(int pos_X, int pos_Y, int pos_Z, double add_tmp);

    const std::vector<std::vector<std::vector<double>>> &get_cumulative_energy() const;

    const std::vector<std::vector<std::vector<double>>> &get_temperature() const;

    //Empty destructor
    ~dust_species(void);

    void set_temperature_at_position(int ix, int iy, int iz, double temperature);
};

#endif //RADMC_DUST_SPECIES_HH

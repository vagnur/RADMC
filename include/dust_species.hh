#ifndef RADMC_DUST_SPECIES_HH
#define RADMC_DUST_SPECIES_HH

#include <vector>
#include <string>

/*
 * An object of this class represent a dust specie present in the simulation
 * Is going to store the relevant information of each specie, as the name, the frequencies and the scattering relevant information.
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

    //This method add the inputted energy in the inputted grid position
    //Input : int iD -> position in the grid in dimension iD (D = x, y or z)
    //        double add_energy -> energy that we are going to add in that position of the grid
    //Output : It has no output, just the modification of the energy vector
    void add_energy(int pos_X, int pos_Y, int pos_Z, double add_energy);

    //This method set the inputted temperature in the inputted grid position
    //Input : int iD -> position in the grid in dimension iD (D = x, y or z)
    //        double temperature -> energy that we are going to add in that position of the grid
    //Output : It has no output, just the modification of the temperature vector
    void set_temperature_at_position(int ix, int iy, int iz, double temperature);

    //This method set a null temperature in the inputted grid position
    //Input : int iD -> position in the grid in dimension iD (D = x, y or z)
    //Output : It has no output, just the modification of the temperature vector
    void set_null_temperature(int ix, int iy, int iz);

    //This method initialize the temperatures and the cumulative energy vectors in 0 for each position in the grid
    //Input : int number_of_points_D -> Number of points in the D dimension (D = x,y or z)
    //Output : It has no output but is going to modify the vectors mentioned before
    void initialize_temperature(int number_of_points_x,int number_of_points_y,int number_of_points_z);

    //Setters
    void set_density(const std::vector<std::vector<std::vector<double>>>& densities);
    void set_lambda(const std::vector<double>& lambda);
    void set_kappa_absorption(const std::vector<double>& kappa_absorption);
    void set_kappa_scattering(const std::vector<double>& kappa_scattering);
    void set_g(const std::vector<double>& g);
    void set_frequency(const std::vector<double>& frequency);
    void set_kappa_absorption_remapped(const std::vector<double>& kappa_absorption_remapped);
    void set_kappa_scattering_remapped(const std::vector<double>& kappa_scattering_remapped);
    void set_g_remapped(const std::vector<double>& g_remapped);
    void set_kappa_absorption_interpoled(const std::vector<double>& kappa_absorption_interpoled);
    void set_kappa_scattering_interpoled(const std::vector<double>& kappa_scattering_interpoled);
    void set_g_interpoled(const std::vector<double>& g_interpoled);

    //Getters
    const std::vector<std::vector<std::vector<double>>>& get_densities() const;
    const std::vector<double>& get_frequency() const;
    const std::vector<double>& get_absoprtion() const;
    const std::vector<double>& get_scattering() const;
    const std::vector<double>& get_g() const;
    const std::vector<double>& get_lambda() const;
    const std::vector<double>& get_kappa_absorption_remapped() const;
    const std::vector<double>& get_kappa_scattering_remapped() const;
    const std::vector<double>& get_g_remapped() const;
    const std::vector<double>& get_kappa_absorption_interpoled() const;
    const std::vector<double>& get_kappa_scattering_interpoled() const;
    const std::vector<double>& get_g_interpoled() const;
    const std::vector<std::vector<std::vector<double>>> &get_cumulative_energy() const;
    const std::vector<std::vector<std::vector<double>>> &get_temperature() const;

    //Empty destructor
    ~dust_species(void);
};

#endif //RADMC_DUST_SPECIES_HH

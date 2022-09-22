#ifndef RADMC_PHOTON_HH
#define RADMC_PHOTON_HH

#include <vector>
#include <random>
#include <common.hh>

class photon {

private:

    photon(int number_of_species, int number_of_stars, const std::vector<double> &star_energy,
           const std::vector<double> &star_position, const std::vector<int> &number_of_grid_points,
           const std::vector<double> &difference);

    int star_source;
    int ray_inu;
    double tau_path_total;
    double tau_path_gone;
    bool on_grid;
    std::vector<double> alpha_A_specie;
    double alpha_A_total;
    std::vector<double> alpha_S_specie;
    double alpha_S_total;
    double alpha_total;
    double albedo;
    double dtau;
    std::vector<double> cumulative_alpha;
    std::vector<double> ray_position;
    std::vector<int> grid_position;
    std::vector<double> prev_ray_position;
    std::vector<int> prev_grid_position;
    std::vector<double> direction;
    std::vector<int> orientation;
    std::vector<double> cell_walls;
    std::vector<double> distance;
    std::vector<double> enerPart;
    std::vector<double> enerCum;
    std::vector<double> tempLocal;
    std::vector<double> dbCumul;
    bool is_scattering;

public:

    photon(void);

    int get_star_source(void);

    void set_star_source(int star_source);

    ~photon(void);

    void get_tau_path();

    void get_random_frequency_inu(const std::vector<double> &star_cumulative_spectrum, int number_of_frequencies);

    void chek_unity_vector(double x, double y, double z);

    void get_random_direction();

    void identify_star(int number_of_stars, const std::vector<double> &luminosities_cum);

    void get_cell_walls(std::vector<double> grid_cell_walls_x, std::vector<double> grid_cell_walls_y,
                        std::vector<double> grid_cell_walls_z);

    void is_on_grid(int number_of_points_x, int number_of_points_y, int number_of_points_z);

    double advance_next_position();

    void calculate_opacity_coefficients(double minor_distance, int number_of_species,
                                        std::vector<std::vector<std::vector<std::vector<double>>>> densities,
                                        std::vector<std::vector<double>> kappa_A,
                                        std::vector<std::vector<double>> kappa_S);

    void walk_next_event(int number_of_species, const std::vector<double> &star_energies);

    int find_specie_to_scattering(int number_of_species);

    void getHenveyGreensteinDirection(const std::vector<double> &g);

    void divideAbsorvedEnergy(int number_of_species, const std::vector<double> &star_energies,
                              const std::vector<std::vector<std::vector<double>>> &density,
                              const std::vector<std::vector<double>> &kappa_A);

    double computeDusttempEnergyBd(const std::vector<double> &dbTemp, std::vector<std::vector<double>> &dbLogEnerTemp,
                                   const std::vector<std::vector<double>> &dbEnerTemp, int number_of_temperatures,
                                   double energy, int iSpec);

    void addTemperatureDecoupled(int number_of_species);

    void pickRandomFreqDb(int number_of_frequencies, int number_of_species, int number_of_temperatures,
                          std::vector<double> dbTemp, std::vector<std::vector<std::vector<double>>> dbCumulNorm);

    void
    do_absorption_event(int number_of_species, std::vector<std::vector<std::vector<std::vector<double>>>> temperatures);

    void walk_full_path();

    double found_point(double x, std::string type, int number_of_points, double diference);
};


#endif //RADMC_PHOTON_HH

#ifndef RADMC_PHOTON_HH
#define RADMC_PHOTON_HH

#include <vector>
#include <random>
#include <common.hh>
#include <star.hh>
#include <dust_species.hh>

class photon {

private:

    int star_source;
    int ray_inu;
    double tau_path_total;
    double tau_path_gone;
    bool on_grid;
    double alpha_A_total;
    double alpha_S_total;
    double alpha_total;
    double albedo;
    double dtau;
    std::vector<double> ray_position;
    std::vector<double> prev_ray_position;
    std::vector<int> grid_position;
    std::vector<int> prev_grid_position;
    std::vector<int> orientation;
    std::vector<double> direction;
    std::vector<double> cell_walls;
    std::vector<double> distance;
    std::vector<double> enerPart;
    std::vector<double> enerCum;
    std::vector<double> tempLocal;
    std::vector<double> dbCumul;
    std::vector<double> alpha_S_specie;
    std::vector<double> alpha_A_specie;
    std::vector<double> cumulative_alpha;
    bool is_scattering;

public:

    photon();

    photon(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
           int number_of_species, int number_of_stars, int number_of_frequencies,
           const std::vector<double> &luminosities_cum, const std::vector<star>& star_information);

    photon(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_zero_one_distribution,
           int number_of_species, int number_of_frequencies, int star_source, int ray_inu,
           const std::vector<double> &star_ray_position);

    int get_star_source(void);

    void set_star_source(int star_source);

    ~photon(void);

    void get_tau_path(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution);

    void get_random_frequency_inu(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
                                  const std::vector<double> &star_cumulative_spectrum, int number_of_frequencies);

    void chek_unity_vector(double x, double y, double z);

    void get_random_direction(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution);

    void identify_star(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
                       int number_of_stars, const std::vector<double> &luminosities_cum);

    void is_on_grid(int number_of_points_x, int number_of_points_y, int number_of_points_z);


    void calculate_opacity_coefficients(double minor_distance, int number_of_species,
                                        const std::vector<dust_species>& dust_species_information);


    int find_specie_to_scattering(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
                                  int number_of_species);

    void get_henvey_greenstein_direction(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
                                      const std::vector<double> &g);

    int found_point(double x, std::string type, int number_of_points, double difference);

    void get_cell_walls(const std::vector<double> &grid_cell_walls_x, const std::vector<double> &grid_cell_walls_y,
                        const std::vector<double> &grid_cell_walls_z);

    bool get_on_grid_condition();

    bool get_is_scattering_condition();

    void divideAbsorvedEnergy(int number_of_species, const std::vector<star>& star_information,
                              const std::vector<dust_species>& dust_specie_information);

    double advance_next_position(int number_of_points_X, int number_of_points_Y, int number_of_points_Z, const std::vector<double> &grid_cell_walls_x,
                                 const std::vector<double> &grid_cell_walls_y,
                                 const std::vector<double> &grid_cell_walls_z);


    void addTemperatureDecoupled(int number_of_species, double cellVolumes,
                                 std::vector<dust_species>& dust_species_information,
                                 const std::vector<double> &dbTemp,
                                 const std::vector<std::vector<double>> &dbLogEnerTemp,
                                 const std::vector<std::vector<double>> &dbEnerTemp, int number_of_temperatures);

    double
    computeDusttempEnergyBd(const std::vector<double> &dbTemp, const std::vector<std::vector<double>> &dbLogEnerTemp,
                            const std::vector<std::vector<double>> &dbEnerTemp, int number_of_temperatures,
                            double energy, int iSpec);

    void pickRandomFreqDb(std::mt19937& generator, std::uniform_real_distribution<>& uniform_real_distribution,
                          int number_of_frequencies, int number_of_species, int number_of_temperatures,
                          const std::vector<double> &dbTemp,
                          const std::vector<std::vector<std::vector<double>>> &dbCumulNorm);



    void
    walk_next_event(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distributionint, int number_of_species, std::vector<dust_species>& dust_species_information,
                    const std::vector<star>& stars_information, int number_of_points_X, int number_of_points_Y, int number_of_points_Z,
                    const std::vector<double> &grid_cell_walls_x, const std::vector<double> &grid_cell_walls_y,
                    const std::vector<double> &grid_cell_walls_z);

    const std::vector<double> &getRayPosition() const;

    void setGridPosition(const std::vector<int> &gridPosition);

    void do_absorption_event(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_real_distribution,
                             int number_of_species, std::vector<dust_species> &dust_species_information,
                             const std::vector<star> &star_information, const std::vector<double> &dbTemp,
                             const std::vector<std::vector<double>> &dbLogEnerTemp,
                             const std::vector<std::vector<double>> &dbEnerTemp, int number_of_temperatures,
                             int number_of_frequencies, double cellVolumes,
                             const std::vector<std::vector<std::vector<double>>> &dbCumulNorm);

    void set_ray_position(std::vector<double> star_ray_position);

    void set_prev_ray_position();

    void set_prev_grid_position();

    double get_dtau();

    double get_tau_path_total();

    double get_tau_path_gone();

    void update_ray_position(double fraction);

    double get_albedo();

    double get_alpha_A_total();

    const std::vector<double> &get_alpha_A_specie() const;

    const std::vector<int> &get_grid_position() const;

    void update_tau_path_gone();

    void set_scattering_state(double rn);

    int get_ray_inu() const;

    void set_ener_part(const std::vector<double> &ener_part);

    void set_temp_local(const std::vector<double> &temp_local);

    const std::vector<double> &get_ener_part() const;

    const std::vector<double> &get_ener_cum() const;

    const std::vector<double> &get_temp_local() const;

    void set_ray_inu(int ray_inu);

    const std::vector<int> &get_prev_grid_position() const;
};


#endif //RADMC_PHOTON_HH

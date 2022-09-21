#ifndef RADMC_PHOTON_HH
#define RADMC_PHOTON_HH

#include <vector>
#include <random>
#include <common.hh>

class photon {

private:

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
    std::vector<double> energy;
    bool is_scattering;

public:

    photon(void);

    int get_star_source(void);

    void set_star_source(int star_source);

    ~photon(void);

};


#endif //RADMC_PHOTON_HH

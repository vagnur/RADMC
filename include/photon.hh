#ifndef RADMC_PHOTON_HH
#define RADMC_PHOTON_HH

#include <vector>
#include <random>
#include <common.hh>
#include <star.hh>
#include <dust_species.hh>

class photon {

private:

    //Represent the index of the star source in the star information vector
    int star_source;
    //Represent the index of the frequency in the frequencies vector
    int frequency_index;
    //TODO : Entender los tau y comentar
    double tau_path_total;
    double tau_path_gone;
    double dtau;
    //Opacity coeficient information
    //TODO : Verificar dónde se usan los vectores. Considerar dejar en clase monte carlo si no es necesario que el foton cargue con esa info
    double alpha_A_total;
    double alpha_S_total;
    double alpha_total;
    double albedo;
    std::vector<double> dbCumul;
    std::vector<double> alpha_S_specie;
    std::vector<double> alpha_A_specie;
    std::vector<double> cumulative_alpha;
    //This bool indicates if the photon is in the grid or not
    bool on_grid;
    //This bool represent (if true) that the photon is going to do a scattering event
    bool is_scattering;
    //Ray position of the photon
    std::vector<double> ray_position;
    //When we are going to move_photon the photon, we need the previous ray position
    std::vector<double> prev_ray_position;
    //Photon position in the grid
    std::vector<int> grid_position;
    //When we are going to move_photon the photon, we need the previous grid position
    std::vector<int> prev_grid_position;
    //Orientation of the photon
    //TODO : dar ejemplo
    std::vector<int> orientation;
    //Direction of the photon
    //TODO : dar ejemplo
    std::vector<double> direction;
    //Cell walls (limits) of the photon's grid position
    std::vector<double> cell_walls;
    //Distance that the photon is going to "walk"
    std::vector<double> distance;
    //Information of the dust specie in the cell of the photn
    //TODO : Dejar en dust specie parece xd
    std::vector<double> dust_specie_energy;
    std::vector<double> cumulative_dust_specie_energy;
    std::vector<double> dust_specie_temperature;

public:

    //Empty constructor
    photon(void);

    //Photon constructor
    //Input : generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    //        int number_of_species -> Number of dust species
    //        int number_of_frequencies -> Number of frequencies
    //        int star_source -> Index of the star that launch the photon
    //        int frequency_index -> Index of the frequency for the photon
    //        vector star_ray_position -> Position of the ray that launch the star
    //Output : A new photon to use in the MC simulation
    photon(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_zero_one_distribution,
           int number_of_species, int number_of_frequencies, int star_source, int frequency_index,
           const std::vector<double> &star_ray_position);

    //
    //Input : generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    void calculate_tau_path(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution);

    //This method check if the photon is still on the grid or not
    //Input number_of_points_D -> Number of points in the D dimension
    //Output: It has no output, but the photon is going to change a bool varialbe to know if still on the grid or not
    void is_on_grid(int number_of_points_x, int number_of_points_y, int number_of_points_z);

    //This method calculate the opacity coefficients that the vector is going to use for scattering and absorption in a cell of the grid
    //Input : minor_distance -> The distance that the photon is going to walk
    //        number_of_species -> Number of dust species
    //        dust_species_information -> Vector with the information of each dust specie
    //Output : No output, but the photon is going to store relevant information in the opacity vectors
    void calculate_opacity_coefficients(double minor_distance, int number_of_species,
                                        const std::vector<dust_species>& dust_species_information);

    //This method found the index of the specie that we are going to use for the scattering
    //Input : generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    //        number_of_species -> Number of dust species
    //Output : Index of the specie that we are going to use in the scattering event
    int find_specie_to_scattering(std::mt19937& generator, std::uniform_real_distribution<>& uniform_zero_one_distribution,
                                  int number_of_species);

    //Method to find the walls of the cell grid where the photon is
    //Input : grid_cell_wals_D - > Vector with the walls of the grid in each D dimension
    //Output : It has no output, but the photon is going to store the walls in each dimension
    void obtain_cell_walls(const std::vector<double> &grid_cell_walls_x, const std::vector<double> &grid_cell_walls_y,
                           const std::vector<double> &grid_cell_walls_z);

    //Method to move_photon the photon to the next position
    //Input : number_of_points_D -> Number of points in each D dimension
    //        grid_cell_wals_D - > Vector with the walls of the grid in each D dimension
    //Output : It return the minor distance (double) and the photon update the relevant vector with the new position
    double advance_next_position(int number_of_points_X, int number_of_points_Y, int number_of_points_Z, const std::vector<double> &grid_cell_walls_x,
                                 const std::vector<double> &grid_cell_walls_y,
                                 const std::vector<double> &grid_cell_walls_z);

    //Method to update the ray position of the photon
    //Input : fraction -> Fraction of the movement
    //Output : It has no output, but the photon is going to update the positions vectors
    void update_ray_position(double fraction);

    //
    void update_tau_path_gone();

    //Getter
    int get_star_source(void) const;
    int get_frequency_index() const;
    bool get_on_grid_condition() const;
    bool get_is_scattering_condition() const;
    double get_dtau() const;
    double get_tau_path_total() const;
    double get_tau_path_gone() const;
    double get_albedo() const;
    double get_alpha_A_total() const;
    const std::vector<double> &get_ray_position() const;
    const std::vector<double> &get_alpha_A_specie() const;
    const std::vector<int> &get_grid_position() const;
    const std::vector<double> &get_dust_specie_energy() const;
    const std::vector<double> &get_cumulative_dust_specie_energy() const;
    const std::vector<double> &get_temp_local() const;
    const std::vector<int> &get_prev_grid_position() const;
    const std::vector<double> &get_direction() const;

    //Setter
    void set_grid_position(const std::vector<int> &grid_position);
    void set_ray_position(std::vector<double> ray_position);
    void set_scattering_state(double rn);
    void set_dust_specie_energy(const std::vector<double> &dust_specie_energy);
    void set_dust_specie_temperature(const std::vector<double> &dust_specie_temperature);
    void set_frequency_index(int frequency_index);
    void set_orientation(const std::vector<int> &orientation);
    void set_direction(const std::vector<double> &direction);
    void set_prev_ray_position();
    void set_prev_grid_position();

    ~photon(void);
};


#endif //RADMC_PHOTON_HH

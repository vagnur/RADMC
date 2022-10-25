#ifndef RADMC_DUST_HH
#define RADMC_DUST_HH

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <dust_species.hh>
#include <common.hh>

/*
 * This class is going to work with the files "dust_density.inp", "dustopac.inp" and "dustkappa_*.inp".
 * From the "dust_density.inp" is going to obtain the density of each specie present in the grid and is going
 * to store according to the type of grid present in the simulation.
 * With the "dustopac.inp" file is going to know the dust species (the names) that are present in the simulation.
 * Each specie present in the dustopac file need to have its own "dustkappa_*.inp" (where * is the name of the
 * dust specie) and is going to obtain information about the specie and their scattering properties.
*/
class dust {

private:

    //Metadata from the "dust_density.inp" file
    //The iformat right now is always 1 and there is no further information about this parameter
    int iformat_dust_density;
    //This is the number of cells in the grid. Basically, is the number of points present for each dust specie
    int number_of_cells;
    //The number of dust species present in the simulation
    int number_of_species;
    //Each dust specie present in the simulation is going to have several information, therefore, we are going to store it
    //in a vector from the class dust_species; so every dust_specie object is going to store the relevant information
    //from each specie.
    std::vector<dust_species> dust_species_information;
    //Metada from the "dustopac.inp" file
    //The iformat is always 2 and there is no further information about the parameter
    int iformat_dustopac;


public:

    //Empty constructor
    dust(void);

    //This method is going to read the density from each dust specie present in the "dust_density.inp" file
    //Input: int number_of_point_C -> Represent the amount of point in the grid in the C coordinate.
    //Ouput: It has no output, but it's going to store a dust_specie object with a 3D vector
    //       that represent the density in the grid for each specie.
    void read_dust_species_density(int number_of_point_x,int number_of_point_y,int number_of_point_z);

    //This method is going to read the "dustopac.inp" file. That file only gives us the names of the dust species,
    //after we read that file, we need to read each "dustkappa_<specie>.inp" file. We are going to do that with
    //the method read_opacities.
    //Input: It has no input.
    //Output: It has no output. But since we are going to call the read_opacities method, we are going to store info in the
    //        dust_species objets.
    void read_opacities_meta(const std::vector<double>& frequencies);

    //This method is going to read each "dustkappa_<specie>.inp" file, and it's going to store the relevant information
    //according to the data presented in the file.
    //Input: int input_style -> Can be 1 or 10 and indicates which type file to open.
    //       string specie_name -> Name of the dust specie.
    //Output: It has no output, but we are going to store vectors with the relevant information in each dust_specie object.
    void read_opacities(int specie_position, int input_style, std::string specie_name, const std::vector<double>& frequencies);

    //This method is going to call the common methods remap and interpolate of the vectors kappa absorption, kappa scattering and g,
    //  for the inputted specie.
    //Input : int specie -> Position of the specie in the main vector dust_species_information (basically, is the specie)
    //        vector frequencies ->
    //        int iformat -> This format indicates which vector are considered
    //Output : it has no output, but each specie is goin to store the remapped vector and the interpoled vector
    void remap_opacities_values(int specie, const std::vector<double>& frequencies, int iformat);

    //This method is going to call a method in each dust specie that is going to initialize the temperature and energy
    //  vectors in 0.
    //Input : int number_of_points_<dimension> -> The number of points in each dimension (x,y or z)
    //Output : it has no output, but each specie is going to store 0 in the temperature and energy vectors
    void initialize_specie_temperature(int number_of_points_x,int number_of_points_y,int number_of_points_z);

    //This method is going to convert the frequencies in the input to the corresponding lambda values
    //Input : vector<double> frequncy -> Represent the frequency points that we want to convert
    //        int number_of_frequecy -> Indicate the number of points that we want to convert
    //Output : It's going to return a vector<double> with the lambda points
    std::vector<double> convert_lambda_to_frequency(std::vector<double> frequency,int number_of_frequency);

    //Getter for the dust species
    const std::vector<dust_species>& get_dust_species(void) const;

    //Getter for the species but not in const type so we can change things in the vector
    std::vector<dust_species>& get_dust_species_to_change(void);

    //Getter for the number of dust species
    int get_number_of_dust_species(void) const;

    //Empty destructor
    ~dust(void);

    void add_energy_specie(int specie, std::vector<int> grid_position, double add_tmp);

    void set_specie_temperature_at_position(int iSpec, int ix, int iy, int iz, double temperature);

    void set_null_temperature(int iSpecie, int i, int j, int k);
};


#endif //RADMC_DUST_HH

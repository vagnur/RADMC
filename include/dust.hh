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
 * With the "dustopac.inp" and the "dustkappa_*.inp" is going to...
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
    void read_opacities_meta(void);

    //This method is going to read each "dustkappa_<specie>.inp" file, and it's going to store the relevant information
    //according to the data presented in the file.
    //Input: int input_style -> Can be 1 or 10 and indicates which type file to open.
    //       string specie_name -> Name of the dust specie.
    //Output: It has no output, but we are going to store vectors with the relevant information in each dust_specie object.
    void read_opacities(int specie_position, int input_style, std::string specie_name);

    //This method is going to convert the frequencies in the input to the corresponding lambda values
    //Input : vector<double> frequncy -> Represent the frequency points that we want to convert
    //        int number_of_frequecy -> Indicate the number of points that we want to convert
    //Output : Its going to return a vector<double> with the lambda points
    std::vector<double> convert_lambda_to_frequency(std::vector<double> frequency,int number_of_frequency);

    //Getter for the dust species
    std::vector<dust_species> get_dust_species(void);

    //Getter for the number of dust species
    int get_number_of_dust_species(void);

    //Empty destructor
    ~dust(void);

};


#endif //RADMC_DUST_HH

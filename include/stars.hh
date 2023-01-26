#ifndef RADMC_STARS_HH
#define RADMC_STARS_HH

//Include of STL libraries
#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <random>
//In order to make a vector to store the information of each star
#include <star.hh>
//In order to execute common functions in this class
#include <common.hh>

/*
 * An object of this class is going to store the information of each star present in the simulation
 * Each star it's a photon source, and we are going to store their relevant information present in the "stars.inp" file
 */
class stars{

private:

    //Number of frequency points present in each star
    int number_of_frequencies;
    //Number of stars present in the simulation
    int number_of_stars;
    //Medata from the star. If 1, the list of wavelengths will instead be a list of frequencies in Herz. In 2 is a list in micron.
    int iformat;
    //Total luminosity of the stars in the simulation
    double total_luminosity;
    //Vector with the cumulative luminosity of the stars in the simulation
    std::vector<double> cumulative_luminosity;
    //Vector with the information of each star
    std::vector<star> stars_information;

public:

    //Empty constructor
	stars(void);

    //Method that is going to read the information present in the "stars.inp" file.
    void read_stars(void);

    //Method that is going to calculate the spectrum of each star.
    //For this, we have 2 options : The star can be treated as a blackbody, and in that case we use the planck function
    //If the star is not a blackbody, then we calculate the spectrum in the surface of the star.
    //Input : vector frequencies -> Vector that contains the frequencies in the system
    //Output : It has no output, but in the star_i it's going to store the calculated spectrum in a vector
    void calculate_spectrum(const std::vector<double>& frequencies);

    //In order to do the simulations, we need the cumulative spectrum of each star
    //Input : vector frequencies -> Vector that contains the frequencies in the system
    //Output : It has no output, but the star_i it's going to store the calculated cumulative spectrum in a vector
    void calculate_cumulative_spectrum(const std::vector<double>& mean_intensity);

    //This method is going to calculate the luminosity of each star and the total luminosity of the system
    //Input : vector mean_intensity -> vector that contains the mean intensity of the frequency
    //Output : It has no output, but we are going to store the cumulative spectrum in the class atribute
    void calculate_total_luminosities(const std::vector<double>& mean_intensity);

    //This method calculate the energy of each star
    //Input : int number of photons -> Number of photons in the simulation
    //Output : It has no output, but each star is going to store it energy
    void calculate_energy(int number_of_photons);

    //
    void fix_luminosities();

    //This method move the star just a little in order to avoid having a star exactly in a wall of the grid
    //Input : vector cell_walls_D -> Wells of the grid in each D dimension
    //Output : It has no output, but each star is going to have a new position
    void jitter_stars(std::vector<double> cell_walls_x, std::vector<double> cell_walls_y, std::vector<double> cell_walls_z);

    //Method that identify the star source from where the photon is going to be launched
    //Input : generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    //Output : Int that is going to be the index of the star in the vector star_information
    int identify_star(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_zero_one_distribution);

    //Getter for the stars
    const std::vector<star>& get_stars_information(void) const;
    const std::vector<double>& get_cumulative_luminosity() const;
    int get_number_of_stars() const;

    //Empty destructor
	~stars(void);

    void fix_spherical_position(std::vector<double> theta_points, int number_of_theta_points);

};

#endif //RADMC_STARS_HH
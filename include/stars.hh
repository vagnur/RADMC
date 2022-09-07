#ifndef RADMC_STARS_HH
#define RADMC_STARS_HH

//Include of STL libraries
#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cmath>
//In order to make a vector to store the information of each star
#include <star.hh>
#include <common.hh>

/*
 * An object of this class is going to store the information of each star present in the simulation
 * Each star it's the photon source, and we are going to store their relevant information present in the "stars.inp" file
 */
class stars{

private:

    //Number of frequency points present in each star
    int number_of_frequencies;
    //Number of stars present in the simulation
    int number_of_stars;
    //Medata from the star. If 1m the list of wavelengths (see below) will instead be a list of frequencies in Herz. In 2 is a list in micron.
    int iformat;
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
    //Input : vector mean_intesity -> Vector that contains the mean intensity (frequency dnu) in the system
    //Output : It has no output, but in the star_i it's going to store the calculated spectrum in a vector
    void calculate_spectrum(std::vector<double> mean_intensity);

    //In order to do the simulations, we need the cumulative spectrum of each star
    //Input : vector mean_intesity -> Vector that contains the mean intensity (frequency dnu) in the system
    //Output : It has no output, but in the star_i it's going to store the calculated cumulative spectrum in a vector
    void calculate_cumulative_spectrum(std::vector<double> mean_intensity);

    void calculate_total_luminosities(std::vector<double> freq_dnu);

    //Getter for the stars
    std::vector<star> get_stars_information(void);

    //Empty destructor
	~stars(void);
};

#endif //RADMC_STARS_HH
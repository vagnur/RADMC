#ifndef RADMC_FREQUENCIES_HH
#define RADMC_FREQUENCIES_HH

//Include of STL libraries
#include <iostream>
#include <vector>
#include <fstream>
#include <common.hh>
#include <random>

/*
 * This class represent the frequencies present in the simulation that are read from the "wavelength_micron.inp"
 * It's going to store the frequencies and the methods are associated with the procurement of the information
 * and the data transformation
 */
class frequencies{

private:

    //Number of frequency points present in the file
    int number_of_frequencies;
    //Vector that is going to store the frequency points present in the file
    std::vector<double> frequency_values;
    //Vector that is going to store the mean intensity values (dnu)
    std::vector<double> mean_intensity_values;

public:

    //Empty constructor
    frequencies(void);

    //Method to read the frequencies file, that can be "wavelenght_micron.inp" or ... (only working for wavelenght for now)
    //Input : It has no input
    //Otput : It has no output, but it's going to store the frequency values in the vector with the same name
    void read_frequencies(void);

    //Method to calculate the mean intensity according to the frequencies that has been read with the read_frequencies method
    //Inpút : It has no input
    //Output : It has no output, but it's going to store the mean intensity values (dnu) in the vector with the same name
    void calculate_mean_intensity(void);

    //This method get a random frequency from the frequencies vector
    //Input : generator and uniform_zer_one_distribution are objects to generate random number (search the object type for more info)
    //        vector star_cumulative_spectrum -> In order to select the frequency, we need to know the cumulative spectrum of the star that
    //                                           is going to launch the photon
    //Output : A integer number that represent the position of the selected frequency in the frequency vector.
    int get_random_frequency(std::mt19937 &generator, std::uniform_real_distribution<> &uniform_zero_one_distribution,
                             const std::vector<double> &star_cumulative_spectrum);

    //Getters
    int get_number_frequency_points() const;
    const std::vector<double>& get_frequencies() const;
    const std::vector<double>& get_mean_intensity() const;

    //Destructor of the class
    ~frequencies(void);
};

#endif //RADMC_FREQUENCIES_HH
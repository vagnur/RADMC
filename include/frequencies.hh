#ifndef RADMC_FREQUENCIES_HH
#define RADMC_FREQUENCIES_HH

//Include of STL libraries
#include <iostream>
#include <vector>
#include <fstream>
#include <common.hh>


/*
 * This class represent the frequencies present in the simulation that are readed from the "wavelength_micron.inp"
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
    //Otput : It has no output, but its going to store the frequency values in the vector with the same name
    void read_frequencies(void);

    //Method to calculate the mean intensity according to the frequencies that has been read with the read_frequencies method
    //Inpút : It has no input
    //Output : It has no output, but its going to store the mean intesity values (dnu) in the vector with the same name
    void calculate_mean_intensity(void);

    //Getter for the number of frequency points
    int get_number_frequency_points();

    //Getter for the frequencies vector
    std::vector<double> get_frequencies();

    //Getter for the mean intensity vector
    std::vector<double> get_mean_intensity();

    //Destructor of the class
    ~frequencies(void);

};

#endif //RADMC_FREQUENCIES_HH
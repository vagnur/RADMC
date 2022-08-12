#ifndef RADMC_GRID_HH
#define RADMC_GRID_HH

//Include of STL libraries
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

class grid{

private:
    //The grid can be regular
    std::string type;
    std::vector<bool> present_points;
    //TODO : Renombrar variables para generelizar con coordenadas polares
    std::vector<double> x_points;
    std::vector<double> y_points;
    std::vector<double> z_points;
    double diference_x;
    double diference_y;
    double diference_z;

public:

    //Empty constructor
    grid(void);

    //Method to initialize a 3d cartesian regular grid
    //TODO : Generalizar el método para grillas regulares
    void initialize_cartesian_regular(void);
    void initialize_spherical_regular(void);

    double found_point(double);

    void calculate_points_delta(void);

    //TODO : Delete this test
    void test_vectors(void);

    //Destructor of the class
    ~grid(void);

};

#endif //RADMC_GRID_HH

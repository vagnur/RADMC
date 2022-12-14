#ifndef RADMC_STAR_HH
#define RADMC_STAR_HH

#include <vector>

/*
 * A star is the source of m_photons for the monte carlo therm simulation
 * Each star has relevant information that is going to be stored in each star object
 */
class star {

private:

    //Radio of the star
    double radio;
    //Mass of the star
    double mass;
    //Luminosity of the star
    double luminosity;
    //Energy of the star
    double energy;
    //Position of the star at the x,y and z component
    std::vector<double> positions;
    //Flux of the star
    std::vector<double> flux;
    //Spectrum of the star
    std::vector<double> spectrum;
    //Cumulative spectrum of the star
    std::vector<double> cumulative_spectrum;

public:

    //Empty constructor
    star(void);

    //Constructor of the class. We read this information in the stars class and then we create each star with this meta information
    star(double radio, double mass, double x_position, double y_position, double z_position);

    //Setters
    void set_luminosity(double luminosity);
    void set_energy(double star_energy);
    void set_cumulative_spectrum(const std::vector<double>& cumulative_spectrum);
    void set_star_position(const std::vector<double>& star_position);
    void set_spectrum(const std::vector<double>& spectrum);
    void set_flux(const std::vector<double>& flux);

    //Getters
    double get_luminosity(void) const;
    double get_energy(void) const;
    const std::vector<double>& get_cumulative_spectrum() const;
    const std::vector<double>& get_star_position() const;
    const std::vector<double>& get_spectrum() const;
    const std::vector<double>& get_flux() const;
    double get_star_radio() const;

    //Empty destructor
    ~star(void);

};


#endif //RADMC_STAR_HH

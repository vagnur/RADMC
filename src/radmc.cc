#include <iostream>
#include <map>
#include <vector>
#include <monte_carlo.hh>
#include <fstream>
#include <frequencies.hh>
#include <dust.hh>
#include <common.hh>
#include <stars.hh>
#include <emissivity.hh>

/*
std::map<std::string,double> read_main_file();
std::vector<std::string> tokenize(std::string s, std::string del = " ");
void f(std::vector<int>& entrada);
void g(std::vector<int> entrada);
 */

//TODO : Gabo del futuro, recuerda verificar el tema de los digitos significativos. Si bien los valores dan,
//TODO : hay que ver cómo afecta el resultado final la presicion de los valores respecto al radmc original

int main() {
    std::cout << "INICIANDO PRUEBAS" << std::endl;

    monte_carlo mc;
    //mc.do_monte_carlo_therm_regular_cartesian_grid();
    mc.do_monte_carlo_therm_regular_cartesian_grid();

    //std::vector<double> prueba = {0.0,0.3,0.8,1.0};
    //std::cout << common::hunt(prueba,4,0.5,2) << std::endl;

    /*
    //Probando emissitivy
    std::cout.precision(17);

    std::map<std::string,double> simulation_parameters = read_main_file();

    frequencies fr;
    fr.read_frequencies();
    fr.calculate_mean_intensity();
    //freq_dnu y freq_nu correctos (exactamente igual)
    std::vector<double> freq_dnu = fr.get_mean_intensity();
    std::vector<double> freq_nu = fr.get_frequencies();
    int number_of_frequencies = fr.get_number_frequency_points();

    dust ds;
    ds.read_dust_species_density(32,32,32);
    ds.read_opacities_meta();
    std::vector<dust_species> species = ds.get_dust_species();
    std::vector<double> lambda = species[0].get_lambda();
    std::vector<double> frequency = species[0].get_frequency();
    std::vector<double> abs = species[0].get_absoprtion();
    std::vector<double> scat = species[0].get_scattering();
    std::vector<double> g = species[0].get_g();
    int number_of_dust_species = ds.get_number_of_dust_species();

    //Testeados, exactamente igual (de paso quiere decir que los vectores iniciales estan bien tambien :D )
    std::vector<double> abs_remap = common::remap_function(abs.size(), frequency, abs, freq_nu.size(),freq_nu,2,1);

    //for (int i = 0; i < abs_remap.size(); ++i) {
    //    std::cout << abs_remap[i] << std::endl;
    //}
    //exit(0);

    std::vector<double> scat_remap = common::remap_function(scat.size(),frequency,scat,freq_nu.size(),freq_nu,2,1);
    std::vector<double> g_remap = common::remap_function(g.size(),frequency,g,freq_nu.size(),freq_nu,1,1);

    //void emissivity::generate_emissitivy_table(std::map<std::string,double> simulation_parameters,int number_of_species, int number_of_frequencies, std::vector<double> kappa_absorption, std::vector<double> freq_nu, std::vector<double> freq_dnu){
    emissivity emis;
    emis.generate_emissivity_table(simulation_parameters,number_of_dust_species,number_of_frequencies,abs_remap,freq_nu,freq_dnu);
    */

    /*
    //PROBANDO INTERPOLATION Y EMISIVIDAD
    std::map<std::string,double> simulation_parameters = read_main_file();

    frequencies fr;
    fr.read_frequencies();
    fr.calculate_mean_intensity();
    //freq_dnu y freq_nu correctos (exactamente igual)
    std::vector<double> freq_dnu = fr.get_mean_intensity();
    std::vector<double> freq_nu = fr.get_frequencies();
    int number_of_frequencies = fr.get_number_frequency_points();

    std::cout.precision(17);

    dust ds;
    ds.read_dust_species_density(32,32,32);
    ds.read_opacities_meta();
    std::vector<dust_species> species = ds.get_dust_species();
    std::vector<double> lambda = species[0].get_lambda();
    std::vector<double> frequency = species[0].get_frequency();
    std::vector<double> abs = species[0].get_absoprtion();
    std::vector<double> scat = species[0].get_scattering();
    std::vector<double> g = species[0].get_g();
    int number_of_dust_species = ds.get_number_of_dust_species();

    std::vector<double> abs_remap = common::remap_function(abs.size(), frequency, abs, freq_nu.size(),freq_nu,2,1);
    std::vector<double> abs_interpol = common::interpolation_function(abs,frequency,abs.size(),freq_nu,number_of_frequencies,abs_remap);

    emissivity emis;
    emis.generate_emissivity_table(simulation_parameters,number_of_dust_species,number_of_frequencies,abs_interpol,freq_nu,freq_dnu);
    emis.compute_derivate(number_of_dust_species,simulation_parameters["ntemp"],number_of_frequencies,freq_dnu);
    */

    /*
     // PROBANDO DIMENSIONES DEL VECTOR EN 3D
    std::vector<std::vector<std::vector<int>>> prueba;// = {{{1,2},{3,4}},{{5,6},{7,8}}};
    prueba.resize(2,std::vector<std::vector<int>>(2,std::vector<int>(2,1)));
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                std::cout << prueba[i][j][k];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    prueba.clear();
    prueba.resize(2,std::vector<std::vector<int>>(2,std::vector<int>(2,3)));
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                std::cout << prueba[i][j][k];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    */


    /*
     // PRUEBAS FRECUENCIA
    frequencies fr;
    fr.read_frequencies();
    fr.calculate_mean_intensity();
    //std::vector<double> prueba = fr.get_frequencies();
    std::vector<double> prueba = fr.get_mean_intensity();
    std::cout.precision(17);
    for (int i = 0; i < prueba.size(); ++i) {
        std::cout << prueba[i] << std::endl;
    }
    */

    /*
    //PRUEBAS DUST DENSITY
    dust ds;
    ds.read_dust_species_density(32,32,32);
    std::vector<dust_species> species = ds.get_dust_species();
    std::vector<std::vector<std::vector<double>>> densities = species[0].get_densities();
    for (int i = 0; i < 32; ++i) {
        for (int j = 0; j < 32; ++j) {
            for (int k = 0; k < 32; ++k) {
                std::cout << densities[i][j][k];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
     */

    //PRRUEBA DUST OTROS
    /*
    dust ds;
    ds.read_dust_species_density(32,32,32);
    ds.read_opacities_meta();
    std::vector<dust_species> species = ds.get_dust_species();
    std::vector<double> frequency = species[0].get_frequency();
    std::vector<double> abs = species[0].get_absoprtion();
    std::vector<double> scat = species[0].get_scattering();
    std::vector<double> g = species[0].get_g();

    frequencies fr;
    fr.read_frequencies();
    fr.calculate_mean_intensity();
    std::vector<double> freq_nu = fr.get_mean_intensity();
    std::vector<double> freq = fr.get_frequencies();

    std::vector<double> abs_remap = common::remap_function(abs.size(), frequency, abs, freq.size(),freq,2,1);
     */

    /*
    //PRUEBAS STARS
    std::map<std::string,double> simulation_parameters = read_main_file();
    frequencies fr;
    fr.read_frequencies();
    fr.calculate_mean_intensity();
    std::vector<double> freq_dnu = fr.get_mean_intensity();
    std::vector<double> freq_nu = fr.get_frequencies();
    stars st;
    st.read_stars();
    st.calculate_spectrum(freq_nu);
    st.calculate_cumulative_spectrum(freq_dnu);
    st.calculate_total_luminosities(freq_dnu);
    st.calculate_energy(simulation_parameters["nphot"]);
    st.fix_luminosities();
     */

    /*
    std::vector<int> prueba = {1,2,3};
    f(prueba);
    for (unsigned int i = 0; i < prueba.size(); ++i) {
        std::cout << prueba[i] << std::endl;
    }
     */
}


    /*
    std::vector<int> prueba = {1,2,3};
    g(prueba);
    return 0;
}

void f(std::vector<int>& entrada){
    entrada.push_back(10);
}

void g(std::vector<int> entrada){
    int number_of_data = entrada.size();
    for (int i = 0; i < number_of_data; ++i) {
        std::cout << entrada[i] << std::endl;
    }
}
     */
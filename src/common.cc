#include <common.hh>

double common::get_light_speed(){
    return light_speed;
}

double common::get_light_speed_microns(){
    return light_speed_microns;
}

std::vector<double> common::remap_function(int number_of_old_points, std::vector<double> old_points, std::vector<double> old_function, int number_of_new_points, std::vector<double> new_points, int elow, int eup){
    //This vector it's going to store the remapped values
    std::vector<double> new_function(number_of_new_points);
    //In order to iterate over the vectors, we need to know the step, the start and the end of each vector (old and new)
    int old_step,old_start,old_end;
    int new_step,new_start,new_end;
    int hunted_value;
    //These values are used in the logarithmic interpolation
    double x,y;
    //First, we need to know if the old points and the new points are ascending or descending.
    //Since this are monotonic functions, we need to compare only 2 points
    //Ascending case for old points
    if(old_points[1] >= old_points[0]){
        old_step = 1;
        old_start = 0;
        old_end = number_of_old_points - 1;
    }
        //Descending case for old points
    else{
        old_step = - 1;
        old_start = number_of_old_points - 1;
        old_end = 0;
    }
    //Ascending case for new points
    if(new_points[1] >= new_points[0]){
        new_step = 1;
        new_start = 0;
        //new_end = number_of_new_points - 1;
        new_end = number_of_new_points;
    }
        //Descending case for new points
    else{
        new_step = -1;
        new_start = number_of_new_points - 1;
        new_end = -1;
    }

    //Now we iterate over the new values, from the lower to the highest value
    //for (int i = new_start; i <= new_end; i = i + new_step) {
    //for(int i = 0; i < number_of_new_points; ++i){
    int i = new_start;
    while(i != new_end){
        //If the new point is lower than the lowest boundary...
        if(new_points[i] < old_points[old_start]){
            //We have several options
            if(elow == 0){
                //If elow = 0, means that the method doesn't allow a point out of bound
                std::cerr << "Out of bound point in remapping function (lower boundary)" << std::endl;
                exit(0);
            }
            else if(elow == 1){
                //If elow = 1, means that we take boundary value of old function
                new_function[i] = old_function[old_start];
            }
            else if(elow == 2){
                //If elow = 2, means that we are going to do a logarithmic extrapolation
                if(number_of_old_points > 1){
                    if(old_function[old_start+old_step] > 0 && old_function[old_start] > 0){
                        x = old_function[old_start+old_step]/old_function[old_start];
                        y = (std::log(new_points[i])-std::log(old_points[old_start])) / (std::log(old_points[old_start+old_step])-std::log(old_points[old_start]));
                        new_function[i] = old_function[old_start] * std::pow(x,y);
                    }
                    else{
                        new_function[i] = 0.0;
                    }
                }
                else{
                    new_function[i] = old_function[old_start];
                }
            }
            else if(elow == 3){
                //If elow = 3, means that if the value it's out of bound it's going to be 0
                new_function[i] = 0.0;
            }
            else if(elow == 4){
                //If elow = 4, we do a smoother version of elow = 3
                if(number_of_new_points > 1){
                    x = (new_points[i]-old_points[old_start]) / (old_points[old_start+old_step]-old_points[old_start]);
                    if(x > -1.0){
                        new_function[i] = (1-std::abs(x))*old_function[old_start];
                    }
                    else{
                        new_function[i] = 0.0;
                    }
                }
                else{
                    new_function[i] = old_function[old_start];
                }
            }
        }
        //If the new point is higher than the highest boundary...
        else if(new_points[i] > old_points[old_end]){
            //Same as before, we have the same options
            if(eup == 0){
                std::cerr << "Out of bound point in remapping function (higher boundary)" << std::endl;
                exit(0);
            }
            else if(eup == 1){
                new_function[i] = old_function[old_end];
            }
            else if(eup == 2){
                if(number_of_old_points > 1){
                    if(old_function[old_end-old_step] > 0.0 && old_function[old_end] > 0.0){
                        x = old_function[old_end-old_step] / old_function[old_end];
                        y = (std::log(new_points[i]) - std::log(old_points[old_end])) / (std::log(old_points[old_end-old_step])-std::log(old_points[old_end]));
                        new_function[i] = old_function[old_end] * std::pow(x,y);
                    }
                    else{
                        new_function[i] = 0.0;
                    }
                }
                else{
                    new_function[i] = old_function[old_end];
                }
            }
            else if(eup == 3){
                new_function[i] = 0.0;
            }
            else if(eup == 4){
                if(number_of_old_points > 1){
                    x = (new_points[i] - old_points[old_end]) / (old_points[old_end-old_step]-old_points[old_end]);
                    if(x > -1.0){
                        new_function[i] = (1-std::abs(x))*old_function[old_end];
                    }
                    else{
                        new_function[i] = 0.0;
                    }
                }
                else{
                    new_function[i] = old_function[old_end];
                }
            }
        }
        else{
            //In any other case, we are in the domain of the old points
            hunted_value = hunt(old_points,number_of_old_points,new_points[i],number_of_old_points);
            if(new_points[i] == old_points[0]){
                //std::cout << "A" << std::endl;
                new_function[i] = old_function[0];
            }
            else if(new_points[i] == old_points[number_of_old_points - 1]){
                //std::cout << "B" << std::endl;
                new_function[i] = old_function[number_of_old_points - 1];
            }
            else{
                if((hunted_value < 0) || (hunted_value >= number_of_old_points - 1)){
                    std::cerr << "Error in remaping function" << std::endl;
                    std::cerr << "Hunt function out of range" << std::endl;
                    exit(0);
                }
                y = (new_points[i] - old_points[hunted_value]) / (old_points[hunted_value+1] - old_points[hunted_value]);
                if((old_function[hunted_value] > 0) && (old_function[hunted_value+1] > 0)){
                    //std::cout << "C" << std::endl;
                    new_function[i] = std::exp((1.0 - y)*std::log(old_function[hunted_value])+y*std::log(old_function[hunted_value+1]));
                }
                else{
                    //std::cout << "D" << std::endl;
                    new_function[i] = ((1-y)*old_function[hunted_value])+(y*old_function[hunted_value+1]);
                }
            }
        }
        i = i + new_step;
    }
    return new_function;
}

std::vector<double> common::interpolation_function(std::vector<double> old_function, std::vector<double> old_points, int number_of_old_points, std::vector<double> new_points, int number_of_new_points, std::vector<double> remapped){
    //double find_dust_kappa_interpol = 0.0;
    std::vector<double> new_function(number_of_new_points);
    double specie_last_frequency = old_points.back();
    double specie_first_frequency = old_points[0];
    double eps=-1.0;
    int inumax, inumin,inu=-1,x;
    double margin = 1e-4;
    for (int i = 0; i < number_of_new_points; ++i) {
        if(((new_points[i] - specie_last_frequency)*(new_points[i] - specie_first_frequency)) < 0.0){
            // Yes, the frequencies lies within the boundaries.
            x = hunt(old_points,number_of_old_points,new_points[i],number_of_old_points);
            if(x < 0 || x >= number_of_old_points - 1){
                std::cerr << "Error in interpolation function" << std::endl;
                std::cerr << "Hunt function out of range" << std::endl;
                exit(0);
            }
            eps = (new_points[i] - old_points[x]) / (old_points[x+1] - old_points[x]);
            if(eps < 0.0 || eps > 1.0){
                std::cerr << "Error in interpolation function at eps calculation" << std::endl;
                exit(0);
            }
            new_function[i] = ((1.0 - eps) * old_function[x]) + eps * old_function[x+1];
        }
        else{
            //The frequencies are out of range of the original grid
            if(new_points[2] > new_points[1]){
                inumax = number_of_new_points-1;
                inumin = 0;
            }
            else{
                inumax = 0;
                inumin = number_of_new_points-1;
            }
            if(new_points[i] <= new_points[inumin]){
                if(new_points[i] <= (1.0-margin)*new_points[inumin]){
                    std::cerr << "Error in interpolation function" << std::endl;
                    std::cerr << "Frequency out of range ..." << std::endl;
                    exit(0);
                }
                else{
                    eps = (new_points[i]-((1.0 - margin)*new_points[inumin])) / (margin*new_points[inumin]);
                    inu = inumin;
                }
            }
            if(new_points[i] >= new_points[inumax]){
                if(new_points[i] >= (1.0+margin)*new_points[inumax]){
                    std::cerr << "Error in interpolation function" << std::endl;
                    std::cerr << "Frequency out of range ..." << std::endl;
                    exit(0);
                }
                else{
                    eps = ((1.0 + margin)*new_points[inumax] - new_points[i]) / (margin*new_points[inumax]);
                    inu = inumax;
                }
            }
            if(inu >= 0){
                if(eps < 0.0) eps = 0.0;
                if(eps > 1.0) eps = 1.0;
                new_function[i] = eps * remapped[inu];
            }
            else{
                x = hunt(new_points,number_of_new_points,new_points[i],inu);
                if(x < 0 || x > number_of_new_points - 1){
                    std::cerr << "Error in interpolation function" << std::endl;
                    std::cerr << "Hunt function out of range" << std::endl;
                    exit(0);
                }
                new_function[i] = (1.0 - eps) * remapped[x] + eps * remapped[x + 1];
            }
        }
    }
    return new_function;
}

int common::hunt(const std::vector<double>& xx, int n, double x, int jlo){
    int jm, jhi, inc;
    int ascnd;

    ascnd = (xx[n-1] >= xx[0]);

    if (jlo <= -1 || jlo > n-1){
        jlo = - 1;
        jhi = n;
    }

    else{
        inc = 1;
        if ((x >= xx[jlo]) == ascnd){
            if(jlo == n - 1) {
                return jlo;
            }
            jhi = jlo + 1;
            while ((x >= xx[jhi]) == ascnd){
                jlo = jhi;
                inc += inc;
                jhi = jlo + inc;
                if(jhi > n - 1){
                    jhi = n;
                    break;
                }
            }
        }
        else{
            if(jlo == 0){
                jlo = - 1;
                return jlo;
            }
            jhi = jlo - 1;
            while ((x < xx[jlo]) == ascnd) {
                jhi = jlo;
                inc <<= 1;
                if (inc >= jhi) {
                    jlo = -1;
                    break;
                } else {
                    jlo = jhi - inc;
                }
            }
        }
    }
    while (jhi -jlo != 1){
        jm = (jhi+jlo) >> 1;
        if((x > xx[jm]) == ascnd){
            jlo = jm;
        }
        else{
            jhi = jm;
        }
    }

    if(x == xx[n-1]){
        jlo = n - 2;
    }

    if(x == xx[0]){
        jlo = 0;
    }

    return jlo;
}

double common::black_body_planck_function(double temperature,double frequency){
    //This function calculates the Blackbody thermal radiation with the planck function
    //The planck function is definided as:
    //B_v(T) = [2 * h * v^3 / c^2] / [e^(h*v/K_b*T)-1]
    //Where h is the planck constant with value = 6.62607015e−34 J⋅Hz^−1, but we want it in erg⋅Hz^-1 = 6.6260701e-27 erg⋅Hz^-1
    //      v is the frequency in Hz
    //      c is the speed of light with value = 299792458 m / s^2, but we want it in cm / s^2 = 29979245800 cm / s^2
    //      K_b is the Boltzmann constant with value = 1.380649e−23 J⋅K^−1, but we want it in erg⋅K^-1 = 1.380649e-16 erg⋅K^-1
    //double h = 6.6260701e-27;
    //double c = 29979245800;
    //double K_b = 380649e-16;
    /*
    double h = 6.6262000e-27;
    double c = 2.9979245800000e10;
    double K_b = 1.3807e-16;
    std::cout << (2*h*std::pow(frequency,3))/std::pow(c,2) << std::endl;
    std::cout << (std::exp(h*frequency/K_b*temperature)-1) << std::endl;
    return ((2*h*std::pow(frequency,3))/std::pow(c,2))/(std::exp(h*frequency/K_b*temperature)-1);
     */
    double xx = 4.7989e-11 * frequency / temperature;
    if (xx > 709.78271) {
        return 0;
    }
    else{
        return 1.47455e-47 * std::pow(frequency,3) / (std::exp(xx)-1) + 1.e-290;
    }
    //TODO : Entiendo la funcion original, pero no entiendo el algoritmo copiado, ¿qué tiene que ver?
}

std::vector<std::string> common::tokenize(std::string s, std::string del){
    int start = 0;
    int end = s.find(del);
    std::string word;
    std::vector<std::string> words;
    while (end != -1) {
        word = s.substr(start, end - start);
        words.push_back(word);
        start = end + del.size();
        end = s.find(del, start);
    }
    word = s.substr(start, end - start);
    words.push_back(word);
    return words;
}
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
    //This values are used in the logarithmic interpolation
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
        new_end = number_of_new_points - 1;
    }
        //Descending case for new points
    else{
        new_step = -1;
        new_start = number_of_new_points - 1;
        new_end = 0;
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
            x = hunt(old_points,number_of_old_points,new_points[i],number_of_old_points/2);
            std::cout << "JLO: " << x << std::endl;
            if(new_points[i] == old_points[0]){
                new_function[i] = old_function[0];
            }
            else if(new_points[i] == old_points[number_of_old_points - 1]){
                new_function[i] = old_function[number_of_old_points - 1];
            }
            else{
                if((x < 0) || (x > number_of_old_points - 1)){
                    std::cerr << "Error in remaping function" << std::endl;
                    exit(0);
                }
                y = (new_points[i] - old_points[x]) / (old_points[x+1] - old_points[x]);
                if((old_function[x] > 0) && (old_function[x+1] > 0)){
                    new_function[i] = std::exp((1.0 - y)*std::log(old_function[i])+y*std::log(old_function[i+1]));
                }
                else{
                    new_function[i] = (1-y)*old_function[x]+y*old_function[x+1];
                }
            }
        }
        i = i + new_step;
    }
    return new_function;
}

int common::hunt(std::vector<double> xx, int n, double x, int jlo){
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
        jlo = - 0;
    }

    return jlo;
}
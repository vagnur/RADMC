#include <common.h>

static double common::get_light_speed(){
    return this -> light_speed;
}

static double common::get_light_speed_microns(){
    return this -> light_speed_microns;
}

static std::vector<double> common::remap_function(int number_of_old_points, std::vector<double> old_points, std::vector<double> old_function, int number_of_new_points, std::vector<double> new_points, int elow, int eup){
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
    for (int i = new_start; i <= new_end; i = i + new_step) {
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
                new_fuction[i] = 0.0;
            }
            else if(elow == 4){
                //If elow = 4, we do a smoother version of elow = 3
                if(number_of_new_points > 1){
                    x = (new_points[i])-old_points[old_start]) / (old_points[old_start+old_step]-old_points[old_start]);
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
            //TODO : Acá se llama a la función hunt (hay que implementarla u________________u)
            // x = hunt(<parameters>)
            x = 10;
        }
    }
}
#include <photon.hh>


photon::photon(int number_of_species, int number_of_stars, const std::vector<double>& star_energy, const std::vector<double>& star_position, const std::vector<int>& number_of_grid_points, const std::vector<double>& difference) {
    this -> orientation.resize(3);
    this -> direction.resize(3);
    this -> distance.resize(3);
    this -> alpha_A_specie.resize(number_of_species);
    this -> alpha_S_specie.resize(number_of_species);
    this -> cumulative_alpha.resize(number_of_species + 1);
    this -> energy.resize(number_of_species);
    this -> idenfify_star(number_of_stars);
    double photon_energy = star_energy[this -> star_source];
    this -> ray_position[0] = star_position[0];
    this -> ray_position[1] = star_position[1];
    this -> ray_position[2] = star_position[2];
    this -> grid_position[0] = this -> found_point(ray_position[0], "cartesian", number_of_grid_points[0], difference[0]);
    this -> grid_position[1] = this -> found_point(ray_position[1], "cartesian", number_of_grid_points[1], difference[1]);
    this -> grid_position[2] = this -> found_point(ray_position[2], "cartesian", number_of_grid_points[2], difference[2]);
    this -> get_random_direction();
    //TODO : Acorde a lo estudiado, es mejor hacer un random entre 0 y el número de frequencias (a primer ojo, no tiene sentido muestrear como se hace)
    this -> get_random_frequency_inu();
    this -> get_tau_path();
    this -> on_grid = true;
}

double photon::found_point(double x, std::string type, int number_of_points, double diference){
    if(type == "cartesian"){
        return std::floor(x/diference) + (number_of_points/2);
    }
    if(type == "spherical"){
        return std::floor(x/diference);
    }
}

void photon::walk_full_path(){
    this -> walk_next_event();
    while (this -> on_grid){
        if (this -> is_scattering){
            if(scattering_mode == 1){
                this -> get_random_direction();
            }
            else if (scattering_mode == 2){
                int specie = this -> find_specie_to_scattering();
                this -> getHenveyGreensteinDirection();
            }
        }
        else{
            this -> do_absoprtion_event();
            this -> get_random_direction();
        }
        this -> get_tau_path();
        this -> walk_next_event();
    }
}

void photon::do_absorption_event(){
    int ix = this -> grid_position[0];
    int iy = this -> grid_position[1];
    int iz = this -> grid_position[2];
    this -> divideAbsorvedEnergy(photon,stars,dustDensity, dustOpacity);
    this -> addTemperatureDecoupled(photon, dustDensity, grid, emissivityDb, dustTemperature);
    for (int i=0 ; i<dustDensity->numSpec ; i++){
        photon->tempLocal[i] = dustTemperature->temperatures[i][iz][iy][ix];
    }
    photon->iFrequency = pickRandomFreqDb(emissivityDb, photon, dustDensity->numSpec, freqData->numFrequencies, photon->tempLocal, photon->enerPart);
}

void photon::addTemperatureDecoupled(int number_of_species){
    int ix = this -> grid_position[0];
    int iy = this -> grid_position[1];
    int iz = this -> grid_position[2];
    double cumen;
    for (int iSpec=0 ; iSpec < number_of_species ; ++iSpec){
        //TODO : Obtener cumulEner de la especie ispec y densidad en ispec y cell
        cumen = cumulEner[iz][iy][ix] / (densities[iz][iy][ix] * cellVolumes);
        //TODO : Este vector inicia en 0 para cada especie en el objeto de polvo...
        //TODO : Hay que cachar si puedo cargar en cada protón y luego sumar todo o depende de las temps anteriores
        temperatures[iSpec][iz][iy][ix] = this -> computeDusttempEnergyBd(emissivityDb, cumen, iSpec);
        //dustTemperature->temperatures[iSpec][iz][iy][ix] = computeDusttempEnergyBd(emissivityDb, cumen, iSpec);
    }
}

double photon::computeDusttempEnergyBd(const std::vector<double>& dbTemp, std::vector<std::vector<double>>& dbLogEnerTemp, const std::vector<std::vector<double>>& dbEnerTemp, int number_of_temperatures, double energy, int iSpec){
    //printf("in computeDusttempEnergyBd\n");
    double tempReturn=0.0;
    int itemp = 0;
    double logEner = std::log(energy);
    float eps;
    //TODO : Obtener dblogenergtemp de la especie
    itemp = common::hunt(dbLogEnerTemp[iSpec], number_of_temperatures, (double) logEner, dbLogEnerTemp[iSpec][number_of_temperatures/2]);
    //printf("itemp=%d\n", *itemp);
    if (itemp >= number_of_temperatures-1){
        std::cerr << "ERROR : Too high temperature discovered" << std::endl;
        exit(0);
    }

    if (itemp <= -1){
        //Temperature presumably below lowest temp in dbase
        //TODO : Obtener enertemp de la especie
        eps = energy/dbEnerTemp[iSpec][0];
        if (eps >= 1){
            //printf("exit\n");
            //exit(1);
            std::cerr << "ERROR : Too high temperature discovered" << std::endl;
            exit(0);
        }
        tempReturn = eps * dbTemp[0];
    }else{
        //TODO : Obtener enertemp de la especie
        eps = (energy - dbEnerTemp[iSpec][itemp]) / (dbEnerTemp[iSpec][itemp+1] - dbEnerTemp[iSpec][itemp]);
        if ((eps>1) || (eps<0)){
            std::cerr << "ERROR : Temperature found out of range..." << std::endl;
            printf("exit\n");
            //exit(1);
        }
        //TODO : Obtener dbTemp
        tempReturn = (1.0 - eps) * dbTemp[itemp] + eps * dbTemp[itemp+1];
    }
    return tempReturn;
}

void photon::divideAbsorvedEnergy(int number_of_species, const std::vector<double>& star_energies, const std::vector<std::vector<std::vector<double>>>& density, const std::vector<std::vector<double>>& kappa_A){
    //printf("in divideAbsorvedEnergy\n");
    int ix = this -> grid_position[0];
    int iy = this -> grid_position[1];
    int iz = this -> grid_position[2];
    double alphaA=0;
    if (number_of_species == 1){
        //TODO : Obtener energia de la estrella
        this -> energy[0] = star_energies[this -> star_source];
    }else{
        for (int i=0 ; i < number_of_species ; ++i){
            //TODO : Obtener densidad de la especie i
            this -> energy[i] = density[iz][iy][ix] * kappa_A[i][this -> ray_inu];
            //photon->enerPart[i] = dustDensity->densities[i][iz][iy][ix] * dustOpacity->kappaA[i][photon->iFrequency];
            alphaA += this -> energy[i];
        }
        for (int i=0 ; i < number_of_species ; ++i){
            this -> energy[i] = star_energies[this -> star_source] * this -> energy[i] / alphaA;
        }
    }
}

void photon::getHenveyGreensteinDirection(const std::vector<double>& g){
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    int valuesOrientations[3] = {0,1,1};
    double newx,newy,newz;
    double temp, g2, sign;
    bool equalZero = true;
    //printf("dustOpacity->g[iSpec][photon->iFrequency]=%10.10lg\n",dustOpacity->g[iSpec][photon->iFrequency]);
    //TODO : Obtener g de la especie obtenida
    double g_value = g[this -> ray_inu];

    //get random numbers != 0.5 , and 2*random-1 != 0
    double rnX = 0.5;
    double rnY = 0.5;
    double rnZ = 0.5;

    while (rnX != 0.5){
        rnX=dis(gen);
    }

    newx = 2*rnX - 1.0;
    if (g_value > 0){
        g2 = g_value*g_value;
        temp = (1.0 - g2) / (1.0 + g_value*newx);
        newx = (1.0 + g2 - temp * temp) / (2.0 * g_value);
        newx = std::max(newx,-1.0);
        newx = std::min(newx,1.0);
    }

    float l2 = 2.0;
    while(l2>1.0){
        rnY=dis(gen);//curand_uniform(&photon->state);
        rnZ=dis(gen);//curand_uniform(&photon->state);
        newy =2*rnY-1;
        newz =2*rnZ-1;
        l2 = newy*newy + newz*newz;
        if(l2 < 0.0001){
            l2 = 2.0;
        }
        //l2 = (l2 < 0.0001 ? 2.0 : l2);
    }

    double linv = std::sqrt((1.0 - newx*newx)/l2);
    //float linv = sqrtf((1.0-newx*newx)/l2);
    newy= newy*linv;
    newz= newz*linv;

    double oldx = this -> direction[0];
    double oldy = this -> direction[1];
    double oldz = this -> direction[2];

    //rotateVector
    double l = std::sqrt((oldx*oldx)+(oldx*oldy));
    //float l = sqrtf(oldx*oldx+oldy*oldy);
    double vx = l*newx-oldz*newz;
    double vy=newy;
    double vz=oldz*newx+l*newz;
    newx = vx;
    newz = vz;
    if (l>1e-6){
        double dx=oldx/l;
        double dy=oldy/l;
        newx=dx*vx-dy*vy;
        newy=dy*vx+dx*vy;
    }
    this -> checkUnitVector(newx,newy,newz);
    //get orientations
    //obtain orientations. It is 0 (left,down) or 1 (right, up)
    int ix = std::floor(newx)+1.0;
    int iy = std::floor(newy)+1.0;
    int iz = std::floor(newz)+1.0;

    this -> orientation[0]=valuesOrientations[ix];
    this -> orientation[1]=valuesOrientations[iy];
    this -> orientation[2]=valuesOrientations[iz];
    //printf("floor: %d, %d, %d\n",ix,iy,iz);
    /*if (ix==2 || iy==2 ||iz==2 ){
      printf("actual: dirx=%lf diry=%lf dirz=%lf\nnew: dirx=%lf diry=%lf dirz=%lf\n",photon->direction[0],photon->direction[1],photon->direction[2],newx,newy,newz);
    }

    if (ix==2 || iy==2 ||iz==2 ){
      checkUnitVector(newx,newy,newz);
      //printf("newx=%2.8lg newy=%2.8lg newz=%2.8lg\n",newx,newy,newz);
    }*/
    this -> direction[0]=newx;
    this -> direction[1]=newy;
    this -> direction[2]=newz;
}

int photon::find_specie_to_scattering(int number_of_species){

    //TODO : Esto pareciera ser general a todos los fotones, si es el caso, dejarlo en la clase superior
    this -> cumulative_alpha[0] = 0.0;

    for (int i = 0; i < number_of_species; ++i) {
        this -> cumulative_alpha[i+1] = this -> cumulative_alpha[i] + this -> alpha_S_specie[i];
    }
    for (int i = 0; i < number_of_species; ++i) {
        this -> cumulative_alpha[i] = this -> cumulative_alpha[i] / this -> cumulative_alpha[number_of_species];
    }

    this -> cumulative_alpha[number_of_species] = 1.0;

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    double rn = dis(gen);
    int iSpec = common::hunt(this -> cumulative_alpha, number_of_species + 1, (double) rn, cumulative_alpha[number_of_species/2]);
    return iSpec;

}

void photon::walk_next_event(int number_of_species, const std::vector<double>& star_energies){
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double minor_distance, fraction;
    bool carry_on = true;
    while(carry_on){
        this -> prev_ray_position[0] = this -> ray_position[0];
        this -> prev_ray_position[1] = this -> ray_position[1];
        this -> prev_ray_position[2] = this -> ray_position[2];

        this -> prev_grid_position[0] = this -> grid_position[0];
        this -> prev_grid_position[1] = this -> grid_position[1];
        this -> prev_grid_position[2] = this -> grid_position[2];

        minor_distance = this -> advance_next_position();

        this -> calculate_opacity_coefficients();

        if(this -> tau_path_gone + this -> dtau > this -> tau_path_total){
            fraction = (this -> tau_path_total - tau_path_gone) / this -> dtau;

            this -> ray_position[0] = this -> prev_ray_position[0] + fraction * (this -> ray_position[0] - prev_ray_position[0]);
            this -> ray_position[1] = this -> prev_ray_position[1] + fraction * (this -> ray_position[0] - prev_ray_position[1]);
            this -> ray_position[2] = this -> prev_ray_position[2] + fraction * (this -> ray_position[0] - prev_ray_position[2]);

            this -> grid_position[0] = this -> prev_grid_position[0];
            this -> grid_position[1] = this -> prev_grid_position[1];
            this -> grid_position[2] = this -> prev_grid_position[2];

            double dum = (1.0 - this -> albedo) * (this -> tau_path_total - this -> tau_path_gone) * star_energies[this -> star_source] / this -> alpha_A_total;
            for (int i=0 ; i < number_of_species ; ++i){
                double addTmp = dum * this -> alpha_A_specie[i];
                //TODO : Hay que ver cómo hacer este vector
                cumulative_energy_specie[i][this -> grid_position[2]][this ->grid_position[1]][this->grid_position[0]] = cumulative_energy_specie[i][this -> grid_position[2]][this ->grid_position[1]][this->grid_position[0]] + add_tmp;
            }
            carry_on = false;
        }
        else{
            double dum = (1.0 - this -> albedo) * this -> dtau * star_energies[this -> star_source] / this -> alpha_A_total;
            for (int i=0 ; i < number_of_species ; ++i){
                double addTmp = dum * this -> alpha_A_specie[i];
                cumulative_energy_specie[i][this -> grid_position[2]][this ->grid_position[1]][this->grid_position[0]] = cumulative_energy_specie[i][this -> grid_position[2]][this ->grid_position[1]][this->grid_position[0]] + add_tmp;
            }
            this -> tau_path_gone =  this -> tau_path_gone + this -> dtau;
            carry_on = this -> on_grid;
        }
    }
    double rn = dis(gen);
    this -> is_scattering = rn < this -> albedo;
}

void photon::calculate_opacity_coefficients(double minor_distance, int number_of_species, std::vector<std::vector<std::vector<std::vector<double>>>> densities, std::vector<std::vector<double>> kappa_A, std::vector<std::vector<double>> kappa_S){
    int ix = this -> prev_grid_position[0];
    int iy = this -> prev_grid_position[0];
    int iz = this -> prev_grid_position[0];

    double opacity_coefficient_alpha_A_total = 0;
    double opacity_coefficient_alpha_S_total = 0;

    for (int i = 0; i < number_of_species; ++i) {
        this -> alpha_A_specie[i] = densities[i][iz][iy][ix] * kappa_A[i][this->ray_inu];
        this -> alpha_S_specie[i] = densities[i][iz][iy][ix] * kappa_S[i][this->ray_inu];
        opacity_coefficient_alpha_A_total = opacity_coefficient_alpha_A_total + this -> alpha_A_specie[i];
        opacity_coefficient_alpha_S_total = opacity_coefficient_alpha_S_total + this -> alpha_S_specie[i];
    }
    this -> alpha_A_total = opacity_coefficient_alpha_A_total;
    this -> alpha_S_total = opacity_coefficient_alpha_S_total;
    this -> alpha_total = opacity_coefficient_alpha_A_total + opacity_coefficient_alpha_S_total;
    this -> albedo = opacity_coefficient_alpha_S_total / this -> alpha_total;
    this -> dtau = this -> alpha_total * minor_distance;

}

double photon::advance_next_position(){
    std::vector<int> signs = {-1,1};
    this -> get_cell_walls();
    this -> distance[0] = (this -> cell_walls[0] - this -> ray_position[0]) / this -> direction[0];
    this -> distance[1] = (this -> cell_walls[1] - this -> ray_position[1]) / this -> direction[1];
    this -> distance[2] = (this -> cell_walls[2] - this -> ray_position[2]) / this -> direction[2];

    double min_distance = std::min(std::min(distance[0],distance[1]), distance[2]);
    int count = 0;
    std::vector<int> indexes = {-1,-1,-1};
    for (int i=0 ; i<3 ; ++i){
        if (this -> distance[i] == min_distance){
            indexes[count]=i;
            count++;
        }
    }

    this -> ray_position[0] = this -> ray_position[0] + min_distance * this -> direction[0];
    this -> ray_position[1] = this -> ray_position[1] + min_distance * this -> direction[1];
    this -> ray_position[2] = this -> ray_position[2] + min_distance * this -> direction[2];

    //avoid bug assign cellWall to ray position
    //update grid position with signs
    for (int i=0 ; i<count ; i++){
        this -> ray_position[indexes[i]] = this -> cell_walls[indexes[i]];
        this -> grid_position[indexes[i]] = this -> grid_position[indexes[i]] + signs[this -> orientation[indexes[i]]];
    }

    this -> is_on_grid();
    return min_distance;
}

void photon::is_on_grid(int number_of_points_x, int number_of_points_y, int number_of_points_z){
    bool on_x = (this -> grid_position[0] >= 0) && (this -> grid_position[0] < number_of_points_x);
    bool on_y = (this -> grid_position[1] >= 0) && (this -> grid_position[1] < number_of_points_y);
    bool on_z = (this -> grid_position[2] >= 0) && (this -> grid_position[2] < number_of_points_z);
    this -> on_grid = on_x && on_y && on_z;
}

void photon::get_cell_walls(std::vector<double> grid_cell_walls_x,std::vector<double> grid_cell_walls_y,std::vector<double> grid_cell_walls_z){
    this -> cell_walls[0] = grid_cell_walls_x[this-> grid_position[0] + this -> orientation[0]];
    this -> cell_walls[1] = grid_cell_walls_y[this -> grid_position[1] + this -> orientation[1]];
    this -> cell_walls[2] = grid_cell_walls_z[this -> grid_position[2] + this -> orientation[2]];
}

void photon::identify_star(int number_of_stars, const std::vector<double>& luminosities_cum){
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    int star_lum = dis(gen);
    //TODO : Obtener luminosities cum
    int star = common::hunt(luminosities_cum, number_of_stars+1, star_lum, luminosities_cum[number_of_stars/2]);
    this -> star_source = star;
}

void photon::get_random_direction(){
    std::vector<int> values_orientations = {0,1,1};
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double l2,dirx,diry,dirz,linv;
    int ix,iy,iz;
    l2 = 2.0;
    while (l2 > 1.0){
        dirx = 2.0 * dis(gen) - 1.0;
        diry = 2.0 * dis(gen) - 1,0;
        dirz = 2.0 * dis(gen) - 1.0;
        l2 = dirx*dirx + diry*diry + dirz*dirz;
        if(l2 < 1e-4){
            l2 = 2.0;
        }
    }
    //TODO : Verificar si el l2 = 0 afecta (ver codigo esteban)
    linv = 1.0 / std::sqrt(l2);

    dirx = dirx * linv;
    diry = diry * linv;
    dirz = dirz * linv;

    this -> chek_unity_vector(dirx,diry,dirz);

    ix = std::floor(dirx) + 1.0;
    iy = std::floor(diry) + 1.0;
    iz = std::floor(dirz) + 1.0;

    this -> orientation[0] = values_orientations[ix];
    this -> orientation[1] = values_orientations[iy];
    this -> orientation[2] = values_orientations[iz];

    this -> direction[0] = dirx;
    this -> direction[1] = diry;
    this -> direction[2] = dirz;
}

void photon::chek_unity_vector(double x, double y, double z){
    double module = std::sqrt(x*x + y*y + z*z);
    if (std::fabs(module - 1.0) > 1e-6){
        std::cerr << "ERROR : Error unity vector " << x << " " << y << " " << z << std::endl;
        exit(0);
    }
}

void photon::get_random_frequency_inu(const std::vector<double>& star_cumulative_spectrum, int number_of_frequencies){
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double rn = dis(gen);
    int ray_inu = common::hunt(star_cumulative_spectrum,number_of_frequencies+1,rn,star_cumulative_spectrum[int(number_of_frequencies/2)]);
    this -> ray_inu = ray_inu;
}

void photon::get_tau_path(){
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double rn = dis(gen);
    this -> tau_path_total = -1.0 * std::log(1.0 - rn);
    this -> tau_path_gone = 0.0;
}

int photon::get_star_source(void){
    return this -> star_source;
}

void photon::set_star_source(int star_source){
    this -> star_source = star_source;
}

photon::~photon(void) {
    ;
}

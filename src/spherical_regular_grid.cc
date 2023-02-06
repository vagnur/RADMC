#include "spherical_regular_grid.hh"

spherical_regular_grid::spherical_regular_grid() {
    ;
}

void spherical_regular_grid::initialize_grid(){
    this -> type = "spherical";
    this -> read_grid_file();
    this -> adjust_theta();
    this -> adjust_phi();
    this -> calculate_points_delta();
    this -> calculate_cell_volume();
    this -> amrray_sint1 = sin(this -> y_points[0]);
    this -> amrray_cost1 = cos(this -> y_points[0]);
    this -> amrray_sint2 = sin(this -> y_points[this -> number_of_points_y-1]);
    this -> amrray_cost2 = cos(this -> y_points[this -> number_of_points_y-1]);
    this -> amrray_sinp1 = sin(this -> z_points[0]);
    this -> amrray_cosp1 = cos(this -> z_points[0]);
    this -> amrray_sinp2 = sin(this -> z_points[this -> number_of_points_z-1]);
    this -> amrray_cosp2 = cos(this -> z_points[this -> number_of_points_z-1]);
    this -> amr_cyclic_xyz.resize(3);
    amr_cyclic_xyz[0] = false;
    amr_cyclic_xyz[1] = false;
    amr_cyclic_xyz[2] = false;
    if(this -> z_points[0] == 0.0 && this -> z_points[this->number_of_points_z-1] == common::get_two_pi()){
        this -> amr_cyclic_xyz[2] = true;
    }
    this->compute_trigometrics_for_edges();
}

void spherical_regular_grid::compute_trigometrics_for_edges(){
    double theta1,theta2,phi1,phi2;
    int nn;
    //
    nn = this -> number_of_points_y;
    this-> amrray_finegrid_sintsq1.resize(nn);
    this-> amrray_finegrid_sintsq2.resize(nn);
    this-> amrray_finegrid_costsq1.resize(nn);
    this-> amrray_finegrid_costsq2.resize(nn);
    for (int i = 0; i < nn; ++i) {
        if(this -> y_points[i] < common::get_pi_half()){
            theta1 = this -> y_points[i];
            theta2 = this -> y_points[i+1];
        }
        else{
            theta1 = this -> y_points[i+1];
            theta2 = this -> y_points[i];
        }
        this -> amrray_finegrid_sintsq1[i] = pow(sin(theta1),2);
        this -> amrray_finegrid_sintsq2[i] = pow(sin(theta2),2);
        this -> amrray_finegrid_costsq1[i] = 1.0 - amrray_finegrid_sintsq1[i];
        this -> amrray_finegrid_costsq2[i] = 1.0 - amrray_finegrid_sintsq2[i];
        if(theta2 == common::get_pi_half()){
            this -> amrray_finegrid_sintsq2[i] = 1.0;
            this -> amrray_finegrid_costsq2[i] = 0.0;
        }
        if(theta1 == 0.0){
            this -> amrray_finegrid_sintsq1[i] = 0.0;
            this -> amrray_finegrid_costsq1[i] = 1.0;
        }
        if(theta2 == common::get_pi()){
            this -> amrray_finegrid_sintsq2[i] = 0.0;
            this -> amrray_finegrid_costsq2[i] = -1.0;
        }
    }
    //
    nn = this -> number_of_points_z;
    this-> amrray_finegrid_sinp1.resize(nn);
    this-> amrray_finegrid_sinp2.resize(nn);
    this-> amrray_finegrid_cosp1.resize(nn);
    this-> amrray_finegrid_cosp2.resize(nn);
    for (int i = 0; i < nn; ++i) {
        phi1 = this -> z_points[i];
        phi2 = this -> z_points[i+1];
        this -> amrray_finegrid_sinp1[i] = sin(phi1);
        this -> amrray_finegrid_sinp2[i] = sin(phi2);
        this -> amrray_finegrid_cosp1[i] = sin(phi1);
        this -> amrray_finegrid_cosp2[i] = sin(phi2);
    }
}

void spherical_regular_grid::adjust_theta() {

    if((this -> y_points[0] != 0.0 ) && (std::abs(this -> y_points[0]) < 1.0e-4 * std::abs(this -> y_points[1] - this -> y_points[0]))){
        std::cout << "Adjusting theta(0) to exactly 0" << std::endl;
        this -> y_points[0] = 0.0;
    }

    if((this -> y_points[this -> number_of_points_y-1] != common::get_pi_half()) && (std::abs(this -> y_points[this -> number_of_points_y-1] - common::get_pi_half()) < 1.0e-4*std::abs(this -> y_points[this -> number_of_points_y-1] - this -> y_points[this -> number_of_points_y - 2]))){
        std::cout << "Adjusting theta(ny) to exactly pi/2..." << std::endl;
        this -> y_points[this -> number_of_points_y-1] = common::get_pi_half();
    }

    if((this -> y_points[this -> number_of_points_y-1] != common::get_pi()) && (std::abs(this -> y_points[this -> number_of_points_y-1] - common::get_pi()) < 1.0e-4*std::abs(this -> y_points[this -> number_of_points_y-1] - this -> y_points[this -> number_of_points_y - 2]))){
        std::cout << "Adjusting theta(ny) to exactly pi..." << std::endl;
        this -> y_points[this -> number_of_points_y-1] = common::get_pi();
    }

    for (int i = 1; i < this -> number_of_points_y - 1; ++i) {
        if((this -> y_points[i] != common::get_pi_half()) && (std::abs(this -> y_points[i]-common::get_pi_half()) < 0.5e-4*std::abs(this -> y_points[i + 1] - this -> y_points[i - 1]))){
            std::cout << "Adjusting theta(" << i <<") to exactly pi/2..." << std::endl;
            this -> y_points[i] = common::get_pi_half();
        }
    }

}

void spherical_regular_grid::adjust_phi() {

    if((this -> z_points[0] != 0.0) && (std::abs(this -> z_points[0]) < 1.0e-2 )){
        std::cout << "Adjusting phi(0) to exactly 0..." << std::endl;
        this -> z_points[0] = 0.0;
    }

    if((this -> z_points[this -> number_of_points_z - 1] != common::get_two_pi()) && (std::abs(this -> z_points[this -> number_of_points_z - 1] - common::get_two_pi()) < 1.0e-2)){
        std::cout << "Adjusting phi(nz) to exactly 2*PI..." << std::endl;
        this -> z_points[this -> number_of_points_z - 1] = common::get_two_pi();
    }

}

void spherical_regular_grid::calculate_cell_volume(){
    this -> cell_volume.resize(this -> number_of_points_z - 1, std::vector<std::vector<double>>(this -> number_of_points_y - 1, std::vector<double>(this -> number_of_points_x - 1)));
    double volumen;

    for (int k = 0; k < this -> number_of_points_z - 1; ++k) {
        for (int j = 0; j < this -> number_of_points_y - 1; ++j) {
            for (int i = 0; i < this -> number_of_points_x - 1; ++i) {
                volumen = (1.0/3.0) * (std::pow(this -> x_points[i+1],3) - std::pow(this -> x_points[i],3)) * (std::cos(this -> y_points[j]) - std::cos(this -> y_points[j+1])) * (this -> z_points[k+1] - this -> z_points[k]);
                this -> cell_volume[k][j][i] = volumen;
            }
        }
    }
}

std::vector<int> spherical_regular_grid::found_ray_position_in_grid(const std::vector<double>& ray_position){
    std::vector<int> grid_points(3);
    if(ray_position[0] >= this -> x_points[0] && ray_position[0] <= this -> x_points[this->number_of_points_x-1]) {
        grid_points[0] = std::floor(ray_position[0] / (this->diference_x + this->x_points[0]));
    }
    else{
        grid_points[0] = -1;
    }
    if(ray_position[1] >= this -> y_points[0] && ray_position[1] <= this -> y_points[this->number_of_points_y-1]) {
        grid_points[1] = std::floor(ray_position[1] / (this->diference_y + this->y_points[0]));
    }
    else{
        grid_points[1] = -1;
    }
    if(ray_position[2] >= this -> z_points[0] && ray_position[2] <= this -> z_points[this->number_of_points_z-1]) {
        grid_points[2] = std::floor(ray_position[2] / (this->diference_z + this->z_points[0]));
    }
    else{
        grid_points[2] = -1;
    }
    return grid_points;
}

void spherical_regular_grid::calculate_points_delta(){
    this->diference_x = this->x_points[1] - this->x_points[0];
    this->diference_y = this->y_points[1] - this->y_points[0];
    this->diference_z = this->z_points[1] - this->z_points[0];
}

spherical_regular_grid::~spherical_regular_grid() {
    ;
}

double spherical_regular_grid::get_cell_volume(std::vector<int> grid_position) const{
    double volume = (1.0/3.0) * (std::pow(this -> x_points[grid_position[0]+1],3) - std::pow(this -> x_points[grid_position[0]],3)) * (std::cos(this -> y_points[grid_position[1]]) - std::cos(this -> y_points[grid_position[1]+1])) * (this -> z_points[grid_position[2]+1] - this -> z_points[grid_position[2]]);
    return volume;
}

void spherical_regular_grid::calculate_photon_cell_walls(photon &photon_i){
    ;
}

//TODO : La grilla esférica tiene que tener un radio que no puede ser 0, por lo que si la estrella
//TODO : parte entre el centro y el radio inicial de la grilla, entonces el fotón estará "fuera" de la grilla
//TODO : Digo fuera, por que estará en el centro que no está considerado en la grilla
//TODO : Por así decirlo, el fotón estará más a la "izquierda" del primer radio conocido.
void spherical_regular_grid::calculate_photon_new_position(photon &photon_i){
    //FOTON DENTRO
    //TODO : Acá cambia la lógica de cuándo el fotón está dentro.
    if(photon_i.get_on_grid_condition()){
        this -> move_photon_inside(photon_i);
    }
    //FOTON FUERA
    //TODO : Esta función considera ambos casos, cuando el fotón está más a la izquierda y cuando está más a la derecha
    //TODO : Los procesos están juntos, por lo que es necesario separarlos
    else{
        //TODO : Puede el fotón venir desde adentro e ingresar al "centro" de la esfera y volver a entrar?
        //Si el foton está afuera, tiene una oportunidad para entrar
        //TODO : Me disgusta esta parte, si la muralla de r inicial está muy lejos el fotón tiene pocas chances de entrar con lo que se pierde haciendo nara
        this->move_photon_outside(photon_i);
    }
}

//void spherical_regular_grid::move_photon_inside(photon &photon_i) {
//    ;
//}

void spherical_regular_grid::move_photon_inside(photon &photon_i){
    int amrray_icross;
    int amrray_ix_next, amrray_iy_next, amrray_iz_next;
    double x0dotdir,ray_dsend,cross_ds, det,r02,r0,theta0,phi0,ds_try,sgnz,cossindum,ct2,st2,ds;
    double pa,pb,pc,eps,sdet,dum1,dum2;
    double eps_thres = 1e-4;
    double small = 1e-12;
    bool topquadrant,crossequator;
    std::vector<double> ray_position;
    std::vector<int> ray_grid_position;
    std::vector<double> direction;
    std::vector<double> spherical_coordinates;
    ray_grid_position = photon_i.get_grid_position();
    ray_position = photon_i.get_ray_position();
    direction = photon_i.get_direction();
    x0dotdir = ray_position[0] * direction[0] + ray_position[1] * direction[1] + ray_position[2] * direction[2];

    amrray_icross = 0;
    ray_dsend = 1e99;
    cross_ds = ray_dsend;

    spherical_coordinates = common::convert_cartesian_to_spherical_coordinates(ray_position[0],ray_position[1],ray_position[2]);
    r02 = spherical_coordinates[0]*spherical_coordinates[0];
    r0 = spherical_coordinates[0];
    theta0 = spherical_coordinates[1];
    phi0 = spherical_coordinates[2];
    // Now try out crossings with all 6 cell walls. Since we are starting
    // within_ a cell the true crossing is the one with the smallest ds.
    // We do not need to check if the crossing happens within the limits
    // of the cell wall, as we had to do above.
    // That is why this is much easier than the determination of the
    // crossings from outside the grid where this simple selection criterion
    // could not be used.

    // Try out a crossing with inner radius cell wall.
    // Only if ray is moving inward. It is always the smaller of the two.
    if(x0dotdir < 0.0){
        det = x0dotdir*x0dotdir + this -> x_points[0]*this -> x_points[0] -r02;
        if(det > 0.0){
            ds_try = -x0dotdir - sqrt(det);
            if(ds_try < 0.0) {
                ds_try = 0.0;
            }
            if(ds_try < cross_ds){
                amrray_icross = 1;
                cross_ds = ds_try;
            }
        }
    }

    // Try out crossing with outer radius.
    // Since we are within the cell the solution must be
    // the largest of the two.
    det = x0dotdir*x0dotdir + this -> x_points[1]*this -> x_points[1] - r02;
    ds_try = -x0dotdir + sqrt(det);
    if(ds_try < 0.0) {
        ds_try = 0.0;
    }
    if(ds_try < cross_ds){
        amrray_icross = 2;
        cross_ds = ds_try;
    }

    if(this -> present_dimensions[1]){
        // First check if we are in the top or bottom quadrant
        if(this -> y_points[this->number_of_points_y-2] < common::get_pi_half()){
            topquadrant = true;
            sgnz = 1.0;
            if(ray_position[2] < 0.0){
                //TODO : Mejorar mensaje
                std::cerr << "Error" << std::endl;
                exit(0);
            }
        }
        else{
            topquadrant = false;
            sgnz = -1.0;
            if(ray_position[2] > 0.0){
                std::cerr << "Error in spherical coordinates" << std::endl;
                std::cerr << ray_position[2] << " " << this -> y_points[0] << " " << r0 << " " << theta0 << " " << phi0 << std::endl;
                exit(0);
            }
        }

        //TODO : Ver este vector y generar
        cossindum = this -> amrray_finegrid_costsq2[ray_grid_position[1]]; // ray_position[1]

        if(cossindum == 0.0){
            // This is simple: just the crossing with the z=0-plane
            if(direction[2] != 0.0){
                ds_try = -ray_position[2] / direction[2];
                if((ds_try > 0.0) && (ds_try < cross_ds)){
                    // Yes, we have a valid crossing
                    if(ray_position[0] < 0.0){
                        // We cross from below
                        amrray_icross = 3;
                        cross_ds = ds_try;
                        crossequator = true;
                    }
                    else{
                        // We cross from above
                        amrray_icross = 4;
                        cross_ds = ds_try;
                        crossequator = true;
                    }
                }
            }
        }
        // Now check for crossing with theta=constant cone: the one farthest
        // away from the midplane. Since we are definitely inside this cell
        // this means that we are currently outside this cone.

        st2 = this -> amrray_finegrid_sintsq1[ray_grid_position[1]];
        ct2 = this -> amrray_finegrid_costsq1[ray_grid_position[1]];

        // Be careful if st2=0.d0: this means that this theta=const wall of
        // the cell is in fact the z-axis, which is infinitely thin and
        // should therefore never yield a hit. Since I do not want to be
        // compiler dependent, I do not want this to depend on the precise
        // rounding-off way of formulae, so I do a real check here. Note
        // that I only have to check this for the theta1 wall, because that
        // is by definition the one closest to the z-axis.

        if(st2 > 0.0){
            // Compute some of the coefficients
            pa = ct2 * (direction[0]*direction[0] + direction[1]*direction[1]) - st2 * direction[2]*direction[2];
            pb = 2.0 * ct2 * (ray_position[0] * direction[0] + ray_position[1] * direction[1]) - 2.0 * st2 * ray_position[2] * direction[2];
            pc = ct2 * (ray_position[0]*ray_position[0] + ray_position[1]*ray_position[1]) - st2*ray_position[2]*ray_position[2];
            // Compute eps == 4*pa*pc/pb^2
            eps = 4.0 * pa * pc / (pb*pb+1e-99);
            //Now check out if there is a solution
            det = pb*pb - 4.0*pa*pc;
            if(det > 0.0){
                if(abs(eps) < eps_thres && pb < 0.0){
                    ds_try = (pc / abs(pb)) * (1.0 + 0.25*eps + 0.125*eps*eps + (5.0/64.0)* pow(eps,3));
                    if(ds_try > 0.0 && ds_try < cross_ds){
                        if(topquadrant){
                            amrray_icross = 3;
                        }
                        else{
                            amrray_icross = 4;
                        }
                    }
                    cross_ds = ds_try;
                }
                else{
                    if(pa != 0.0){
                        sdet = sqrt(det);
                        ds_try = -0.5 * (pb + sdet) / pa;
                        if((ds_try > 0.0) && (ds_try < cross_ds)){
                            if(topquadrant){
                                amrray_icross = 3;
                            }
                            else{
                                amrray_icross = 4;
                            }
                            cross_ds = ds_try;
                        }
                    }
                }
            }
        }
        //
        st2 = this -> amrray_finegrid_sintsq2[ray_grid_position[1]];
        ct2 = this -> amrray_finegrid_costsq2[ray_grid_position[1]];

        if(ct2 != 0.0){
            pa = ct2 * (direction[0]*direction[0] + direction[1]*direction[1]) - st2*direction[2]*direction[2];
            pb = 2.0 * ct2 * (ray_position[0]*direction[0] + ray_position[1]*direction[1]) - 2.0 * st2 * ray_position[2]*direction[2];
            pc = ct2 * (ray_position[0]*ray_position[0] + ray_position[1]*ray_position[1]) - st2*ray_position[2]*ray_position[2];
            eps = 4.0 * pa * pc / (pb*pb+1e-99);
            det = pb*pb - 4.0*pa*pc;
            if(det > 0.0){
                if(abs(eps) < eps_thres && pb > 0.0){
                    ds_try = -(pc/pb) * (1.0 + 0.25*eps + 0.125*eps*eps + (5.0/64.0)*pow(eps,3));
                    if(ds_try > 0.0 && ds_try < cross_ds){
                        if(topquadrant) {
                            amrray_icross = 4;
                        }
                        else{
                            amrray_icross = 3;
                        }
                        cross_ds = ds_try;
                    }
                }
                else{
                    if(pa != 0.0){
                        sdet = sqrt(det);
                        ds_try = -0.5 * (pb - sdet) / pa;
                        if(ds_try > 0.0 && ds_try < cross_ds){
                            if(topquadrant) {
                                amrray_icross = 4;
                            }
                            else{
                                amrray_icross = 3;
                            }
                            cross_ds = ds_try;
                        }
                    }
                }
            }
        }
    }

    if(present_dimensions[2]){
        dum1 = ray_position[0] * this -> amrray_finegrid_sinp1[ray_grid_position[2]] - ray_position[1]*this -> amrray_finegrid_cosp1[ray_grid_position[2]];
        dum2 = direction[1] * this -> amrray_finegrid_cosp1[ray_grid_position[2]] - direction[0] * this -> amrray_finegrid_sinp1[ray_grid_position[2]];
        if(dum2 < -small){
            ds_try = dum1 / dum2;
        }
        else{
            ds_try = -1.0;
        }
        if(ds_try >= 0.0 && ds_try < cross_ds){
            amrray_icross = 5;
            cross_ds = ds_try;
        }
        dum1 = ray_position[0] * this -> amrray_finegrid_sinp2[ray_grid_position[2]] - ray_position[1]*this -> amrray_finegrid_cosp2[ray_grid_position[2]];
        dum2 = direction[1] * this -> amrray_finegrid_cosp2[ray_grid_position[2]] - direction[0] * this -> amrray_finegrid_sinp2[ray_grid_position[2]];
        if(dum2 > small){
            ds_try = dum1 / dum2;
        }
        else{
            ds_try = -1.0;
        }
        if(ds_try >= 0.0 && ds_try < cross_ds){
            amrray_icross = 6;
            cross_ds = ds_try;
        }
    }

    ray_position[0] = ray_position[0] + cross_ds * direction[0];
    ray_position[1] = ray_position[1] + cross_ds * direction[1];
    ray_position[2] = ray_position[2] + cross_ds * direction[2];

    if(this -> present_dimensions[1]){
        if(topquadrant){
            if(ray_position[2] < 0){
                ray_position[2] = 0.0;
            }
        }
        else{
            if(ray_position[2] > 0){
                ray_position[2] = 0.0;
            }
        }
        if(crossequator){
            if(topquadrant){
                if(amrray_icross != 4){
                    crossequator = false;
                }
            }
            else{
                if(amrray_icross != 3){
                    crossequator = false;
                }
            }
            if(crossequator){
                ray_position[2] = 0.0;
            }
        }
    }
    int idir,ilr;
    switch (amrray_icross) {
        case(0):
            amrray_ix_next = ray_grid_position[0];
            amrray_iy_next = ray_grid_position[1];
            amrray_iz_next = ray_grid_position[2];
            //TODO : Aqui hace algo pa ver la celda donde queda, verificar como hacerlo pa nuestra lógica
            break;
        case(1):
            idir = 1;
            ilr = 2;
            amrray_ix_next = ray_grid_position[0] - 1;
            amrray_iy_next = ray_grid_position[1];
            amrray_iz_next = ray_grid_position[2];
            //TODO : Aqui hace algo pa ver la celda donde queda, verificar como hacerlo pa nuestra lógica
            if(amrray_ix_next < 0){
                amrray_ix_next = -1;
                amrray_iy_next = -1;
                amrray_iz_next = -1;
            }
            break;
        case(2):
            idir = 1;
            ilr = 1;
            amrray_ix_next = ray_grid_position[0] + 1;
            amrray_iy_next = ray_grid_position[1];
            amrray_iz_next = ray_grid_position[2];
            //TODO : Aqui hace algo pa ver la celda donde queda, verificar como hacerlo pa nuestra lógica
            if(amrray_ix_next > this -> number_of_points_x - 1){
                amrray_ix_next = -1;
                amrray_iy_next = -1;
                amrray_iz_next = -1;
            }
            break;
        case(3):
            idir = 2;
            ilr = 2;
            amrray_ix_next = ray_grid_position[0];
            amrray_iy_next = ray_grid_position[1] - 1;
            amrray_iz_next = ray_grid_position[2];
            //TODO : Aqui hace algo pa ver la celda donde queda, verificar como hacerlo pa nuestra lógica
            if(amrray_iy_next < 0){
                amrray_ix_next = -1;
                amrray_iy_next = -1;
                amrray_iz_next = -1;
            }
            break;
        case(4):
            idir = 2;
            ilr = 1;
            amrray_ix_next = ray_grid_position[0];
            amrray_iy_next = ray_grid_position[1] + 1;
            amrray_iz_next = ray_grid_position[2];
            //TODO : Aqui hace algo pa ver la celda donde queda, verificar como hacerlo pa nuestra lógica
            if(amrray_iy_next > this -> number_of_points_y - 1) {
                amrray_ix_next = -1;
                amrray_iy_next = -1;
                amrray_iz_next = -1;
            }
            break;
        case(5):
            idir = 3;
            ilr = 2;
            amrray_ix_next = ray_grid_position[0];
            amrray_iy_next = ray_grid_position[1];
            amrray_iz_next = ray_grid_position[2] - 1;
            //TODO : Aqui hace algo pa ver la celda donde queda, verificar como hacerlo pa nuestra lógica
            if(amrray_iz_next < 0){
                if(this->amr_cyclic_xyz[2]){
                    amrray_iz_next = this -> number_of_points_z - 1;
                }
                else {
                    amrray_ix_next = -1;
                    amrray_iy_next = -1;
                    amrray_iz_next = -1;
                }
            }
            break;
        case(6):
            idir = 3;
            ilr = 1;
            amrray_ix_next = ray_grid_position[0];
            amrray_iy_next = ray_grid_position[1];
            amrray_iz_next = ray_grid_position[2] + 1;
            //TODO : Aqui hace algo pa ver la celda donde queda, verificar como hacerlo pa nuestra lógica
            if(amrray_iz_next > this -> number_of_points_z - 1) {
                if(this -> amr_cyclic_xyz[2]){
                    amrray_iz_next = 0;
                }
                else {
                    amrray_ix_next = -1;
                    amrray_iy_next = -1;
                    amrray_iz_next = -1;
                }
            }
            break;
    }
    //TODO : Calcular distancia
    ds = sqrt(pow((ray_position[0] - photon_i.get_prev_grid_position()[0]),2) + pow((ray_position[1] - photon_i.get_prev_grid_position()[1]),2) + pow((ray_position[2] - photon_i.get_prev_grid_position()[2]),2));
    photon_i.set_min_distance(ds);
    std::vector<int> grid_new_position = {amrray_ix_next, amrray_iy_next, amrray_iz_next};
    photon_i.set_grid_position(grid_new_position);
    photon_i.set_ray_position(ray_position);
    photon_i.is_on_grid(this->number_of_points_x,this->number_of_points_y,this->number_of_points_z);
}

//TODO : Considero esta función inecesaria, el fotón comienza con una dirección aleatoria.
//TODO : Se podría simplemente iniciar en un valor aleatorio de theta y phi , siempre va a entrar por la izquierda al
//TODO : radio más pequeño y por la derecha al radio más grande.
//TODO : Habría que ver qué pasa si la estrella es una esfera, pero si es puntual pareciera que podemos hacer lo anterior
void spherical_regular_grid::move_photon_outside(photon &photon_i){
    //TODO : La función simplemente es una traducción directa del código original.
    //TODO : Necesita ser rediseñada, es claro que hay procesos repetidos que se pueden generalizar
    //TODO : Es claro que hay un exceso de variables que se pueden generalizar en vectores
    //TODO : Para entender la función buscar intersección de una línea con una esfera.
    //TODO : Notar que la grilla esférica considera, implicitamente, que tiene centro en el origen (0,0,0)
    std::vector<double> ray_position;
    std::vector<double> direction;
    std::vector<double> spherical_coordinates;
    int amrray_icross;
    double ray_dsend = 1e99;
    double x0dotdir, x00dotdir, s00,x00,y00,z00, r00, r002, dum, r2, det, ds_r1,ds_r2,ds_t1,ds_t2,ds_p1,ds_p2, cross_ds, dummy;
    double r02,r0,theta0,phi0;
    double x_r1, y_r1, z_r1, x_r2, y_r2, z_r2, x_t1, y_t1, z_t1, x_t2, y_t2, z_t2, x_p1, y_p1, z_p1, x_p2, y_p2, z_p2;
    double r_r1,theta_r1,phi_r1,r_r2,theta_r2,phi_r2,r_t1,theta_t1,phi_t1,r_t2,theta_t2,phi_t2,r_p1,theta_p1,phi_p1,r_p2,theta_p2,phi_p2;
    // We want to know the edges that the photon can intercept, so we create booleans to know the edges that the photon can intercept
    bool val_r1 = true, val_r2 = true, val_t1 = true, val_t2 = true, val_p1 = true, val_p2 = true;
    //We get the position, direction and spherical coordinates of the photon
    ray_position = photon_i.get_ray_position();
    direction = photon_i.get_direction();
    spherical_coordinates = common::convert_cartesian_to_spherical_coordinates(ray_position[0],ray_position[1],ray_position[2]);
    r02 = spherical_coordinates[0]*spherical_coordinates[0];
    r0 = spherical_coordinates[0];
    theta0 = spherical_coordinates[1];
    phi0 = spherical_coordinates[2];
    //We calculate the dot product of the position and the direction
    x0dotdir = ray_position[0] * direction[0] + ray_position[1] * direction[1] + ray_position[2] * direction[2];
    if(x0dotdir >= 0.0){
        // Ray is outward pointing, so no shift is necessary
        s00 = 0.0;
        x00 = ray_position[0];
        y00 = ray_position[1];
        z00 = ray_position[2];
    }
    else{
        // Ray is inward pointing, so do the shift for ensuring stability of
        // the quadratic equations in case of zooming into extremely small
        // regions close to the center
        s00 = spherical_coordinates[0];
        x00 = ray_position[0] + s00 * direction[0];
        y00 = ray_position[1] + s00 * direction[1];
        z00 = ray_position[2] + s00 * direction[2];
    }
    x00dotdir = x00 * direction[0] + y00 * direction[1] + z00 * direction[2];
    r002 = x00*x00 + y00*y00 + z00*z00;
    r00 = sqrt(x00*x00 + y00*y00 + z00*z00);
    // If r00 is very very much smaller than r0, then the remaining non-zeroness
    // is likely a numerical error, and in fact a perfectly radial inward
    // ray is in fact meant. In that case do not attempt to compute crossings
    // with theta and phi coordinates, because the formulae would become
    // singular.

    // If tiny is too small, so the photon can only intercept with the r edge
    if(r00 < tiny * spherical_coordinates[0]){
        val_t1 = false;
        val_t2 = false;
        val_p1 = false;
        val_p2 = false;
    }
    // Same is true for radially outward moving rays
    dum = pow(ray_position[0] - r0 * direction[0],2) +
          pow(ray_position[1] - r0 * direction[1],2) +
          pow(ray_position[2] - r0 * direction[2],2);
    if(dum < pow(tiny * r0,2)){
        val_t1 = false;
        val_t2 = false;
        val_p1 = false;
        val_p2 = false;
    }

    // If cyclic boundary conditions in phi are chosen, or if phi is not
    // considered, then we do not need to find crossings with the phi
    // grid boundaries
    
    if(this -> amr_cyclic_xyz[2]){
        val_p1 = false;
        val_p2 = false;
    }

    // Now solve the equation for the crossing with inner R=const surface
    // and check if the crossing happens within the theta, phi grid

    // s = -(x00.dirx) +- sqrt( (x00.dirx)^2 + r^2 - r00^2 ) + s00

    // s00 must be added because we do the calculation now from a
    // different point along the ray (not the true starting point).  Without
    // adding the s00 we obtain the s as measured from this different
    // point. By adding s00 we get the actual s, as measured from the true
    // starting point. The reason for this complicated stuff is that the
    // quadratic equations may become inaccurate if the ray points very
    // sharply toward the center and crosses the grid only very close to the
    // center. If the ratio of the radius where this crossing happens and the
    // radius where the ray starts is larger than about 1d-8, then the square
    // of this is 1d-16, meaning that we are in the numerical noise of the
    // double precision. A ratio of 1d-6 is therefore the most extreme one can
    // afford UNLESS we do the shift as described above. The shift is a linear
    // calculation, so with a ratio of 1d-8 we still have 8 digits of accuracy
    // left in the shift. Once the shift is done, the quadratic equations are
    // "renormalized" and such problems with accuracy do no longer take
    // place. So when would we need this? Well, 1 parsec / 1 Rsun = 4.4d7, so
    // if we model a huge cloud around a star we may already get into this
    // regime of dynamic range, and the above shift/renormalization of the
    // equations becomes necessary.

    this -> cross_inner_radius(x_r1,y_r1,z_r1,val_r1,ds_r1,x00dotdir,r002,s00,ray_position,direction,r_r1,theta_r1,phi_r1);

    // Now solve the equation for the crossing with outer R=const surface
    // and check if the crossing happens within the theta, phi grid
    //
    // s = -(x00.dirx) +- sqrt( (x00.dirx)^2 + r^2 - r00^2 ) + s00
    //
    // The same issues with stiffness are accounted for here, although
    // it is not expected to be necessary

    this -> cross_outer_radius(x_r2,y_r2,z_r2,val_r2,ds_r2,x00dotdir,r002,r02,s00,ray_position,direction,r_r2,theta_r2,phi_r2);
    // Now let us try to find the crossings with the theta1=const cone.

    if(val_t1){
        this -> cross_inner_theta_cone(x_t1,y_t1,z_t1,val_t1,ds_t1,ray_position,direction,x00,y00,z00,s00,r_t1,theta_t1,phi_t1);
    }
    // Now let us try to find the crossings with the theta2=const cone.
    if(val_t2){
        this ->cross_outer_theta_cone(x_t2,y_t2,z_t2,val_t2,ds_t2,s00,x00,y00,z00,ray_position,direction,r_t2,theta_t2,phi_t2);
    }
    // Now let us try to find the crossings with the phi1=const cone.
    if(val_p1){
        this ->cross_inner_phi_cone(x_p1,y_p1,z_p1,val_p1,ds_p1,ray_position,direction,r_p1,theta_p1,phi_p1);
    }
    // Now let us try to find the crossings with the phi2=const cone.
    if(val_p2){
        this ->cross_outer_phi_cone(x_p2,y_p2,z_p2,val_p2,ds_p2,ray_position,direction,r_p2,theta_p2,phi_p2);
    }
    // Now determine which, if any, of the crossings is the first one
    // that will be encountered.


    //Determine if/where the ray would cross the grid

    amrray_icross = 0;
    cross_ds = 1e99;
    if(val_r1) {
        amrray_icross = 1;
        cross_ds = ds_r1;
    }
    if(val_r2 && (ds_r2 < cross_ds)){
        amrray_icross   = 2;
        cross_ds = ds_r2;
    }
    if(val_t1 && (ds_t1 < cross_ds)){
        amrray_icross = 3;
        cross_ds = ds_t1;
    }
    if(val_t2 && (ds_t2 < cross_ds)){
        amrray_icross = 4;
        cross_ds = ds_t2;
    }
    if(val_p1 && (ds_p1 < cross_ds)){
        amrray_icross = 5;
        cross_ds = ds_p1;
    }
    if(val_p2 && (ds_p2 < cross_ds)) {
        amrray_icross = 6;
        cross_ds = ds_p2;
    }

    // If amrray_icross==0 then we do not enter the grid at all, so we have arrived

    if(amrray_icross == 0){
        //We do not enter the grid at all
        //TODO : Verificar retorno, pero bajo la lógica actual hay que evitar que el fotón trabaje
        //TODO : Dejar fuera de la grilla por el radio externo?
    }

    // We could enter the grid, but we have to check if the cross_ds is smaller than ray_dsend

    if(cross_ds > ray_dsend){
        // Nope, we arrived at the end point
        ray_position[0] = ray_position[0] + ray_dsend * ray_position[0];
        ray_position[1] = ray_position[1] + ray_dsend * ray_position[1];
        ray_position[2] = ray_position[2] + ray_dsend * ray_position[2];
        //TODO : Verificar retorno, pero bajo la lógica actual hay que evitar que el fotón trabaje
        //TODO : Dejar fuera de la grilla por el radio externo?
    }

    if(amrray_icross == 1) {
        dummy = r_r1 / sqrt(x_r1 * x_r1 + y_r1 * y_r1 + z_r1 * z_r1);
        this ->entering_from_radius(photon_i,dummy,0,ray_position,x_r1,y_r1,z_r1,r_r1,theta_r1,phi_r1);
        amrray_icross = -1;
    }
    else if(amrray_icross == 2){
        dummy = 1.0;
        this ->entering_from_radius(photon_i,dummy,this->number_of_points_x-1,ray_position,x_r2,y_r2,z_r2,r_r2,theta_r2,phi_r2);
        amrray_icross = -2;
    }
    else if(amrray_icross == 3){
        this ->entering_from_theta(photon_i,0,ray_position,x_t1,y_t1,z_t1,r_t1,phi_t1);
        amrray_icross = -3;
    }
    if(amrray_icross == 4){
        this ->entering_from_theta(photon_i,this->number_of_points_y-1,ray_position,x_t2,y_t2,z_t2,r_t2,phi_t2);
        amrray_icross = -4;
    }
    if(amrray_icross == 5){
        this ->entering_from_phi(photon_i,0,ray_position,x_p1,y_p1,z_p1,r_p1,theta_p1);
        amrray_icross = -5;
    }
    if(amrray_icross == 6){
        this ->entering_from_phi(photon_i,this -> number_of_points_z-1,ray_position,x_p2,y_p2,z_p2,r_p2,theta_p2);
        amrray_icross = -6;
    }
    photon_i.is_on_grid(this -> number_of_points_x, this -> number_of_points_y, this -> number_of_points_z);
    //TODO : Generalizar los métodos, probar y completar para la pos del fotón
    //TODO : Verificar camino de fotones que entran y de fotones que salen
}

void spherical_regular_grid::cross_inner_radius(double &x_r1, double &y_r1, double &z_r1,bool &val_r1, double &ds_r1,double x00dotdir, double r002, double s00, const std::vector<double> &ray_position, const std::vector<double> &direction, double &r_r1, double &theta_r1, double &phi_r1){
    double r2,det;
    r2 = this -> x_points[0] * this -> x_points[0];
    det = x00dotdir*x00dotdir + r2 - r002;
    if(det <= 0.0 ){
        ds_r1 = -1e99;
    }
    else{
        ds_r1 = -x00dotdir + sqrt(det) + s00;
    }
    if(ds_r1 <= 0.0){
        val_r1 = false;
    }
    else{
        x_r1 = ray_position[0] + ds_r1 * direction[0];
        y_r1 = ray_position[1] + ds_r1 * direction[1];
        z_r1 = ray_position[2] + ds_r1 * direction[2];
        r_r1 = this -> x_points[0];
        theta_r1 = acos(z_r1/r_r1);
        if(theta_r1 <= this -> y_points[0] - this -> tiny || theta_r1 >= this -> y_points[this->number_of_points_y-1] + tiny){
            val_r1 = false;
        }

        //Now compute phi
        phi_r1 = common::convert_cartesian_to_spherical_coordinates(x_r1,y_r1,z_r1)[2];

        if(not(this -> amr_cyclic_xyz[2])) {
            if (phi_r1 < this->z_points[0] - this->tiny) {
                val_r1 = false;
            }
            if (phi_r1 > this->z_points[this->number_of_points_z-1] + this->tiny) {
                val_r1 = false;
            }
        }
    }
}

void spherical_regular_grid::cross_outer_radius(double &x_r2, double &y_r2, double &z_r2,bool &val_r2, double &ds_r2,double x00dotdir, double r002, double r02, double s00, const std::vector<double> &ray_position, const std::vector<double> &direction,double &r_r2,double &theta_r2,double &phi_r2){
    double r2,det;
    r2 = this -> x_points[this -> number_of_points_x-1]*this -> x_points[this -> number_of_points_x-1];
    det = x00dotdir*x00dotdir + r2 - r002;
    ds_r2 = -1e99;
    if(det > 0.0){
        if(x00dotdir < 0.0 && r2 <= r02){
            ds_r2 = -x00dotdir - sqrt(det) + s00;
        }
    }
    if(ds_r2 < 0.0){
        val_r2 = false;
    }
    else {
        x_r2 = ray_position[0] + ds_r2 * direction[0];
        y_r2 = ray_position[1] + ds_r2 * direction[1];
        z_r2 = ray_position[2] + ds_r2 * direction[2];

        r_r2 = this->x_points[this->number_of_points_x-1];

        theta_r2 = acos(z_r2 / r_r2);

        phi_r2 = common::convert_cartesian_to_spherical_coordinates(x_r2, y_r2, z_r2)[2];

        if(not(this -> amr_cyclic_xyz[2])) {
            if (phi_r2 < this->z_points[0] + this->tiny) {
                val_r2 = false;
            }
            if (phi_r2 > this->z_points[this->number_of_points_z-1] + this->tiny) {
                val_r2 = false;
            }
        }
    }
}

void spherical_regular_grid::cross_inner_theta_cone(double &x_t1, double &y_t1, double &z_t1,bool &val_t1, double &ds_t1,const std::vector<double> &ray_position, const std::vector<double> &direction, double x00, double y00, double z00, double s00,double &r_t1,double &theta_t1,double &phi_t1){
    double st2,ct2,pa,pb,pc,det,sdet,pabc;
    st2 = this -> amrray_sint1 * this -> amrray_sint1;
    ct2 = this -> amrray_cost1 * this -> amrray_cost1;
    pa = ct2 * (direction[0]*direction[0] + direction[1]*direction[1]) - st2*direction[2]*direction[2];
    pb = 2*ct2*(direction[0]*x00 + direction[1]*y00) - 2*st2*direction[2]*z00;
    pc = ct2*(x00*x00+y00*y00)-st2*z00*z00;
    det = pb*pb + 4*pa*pc;
    if(det <= 0.0){
        ds_t1 = -1e99;
        val_t1 = false;
    }
    else{
        sdet = sqrt(det);
        pabc = 4*pa*pc/(pb*pb);
        if(direction[0]*ray_position[0] + direction[1]*ray_position[1] <= 0.0){
            if(fabs(pabc) > 1e-4 || pb < 0.0){
                ds_t1 = -0.5 * (pb - sdet ) / pa + s00;
            }
            else{
                ds_t1 = -0.5 * pb * (0.5*pabc + (1.0/8.0)*pabc*pabc + (1.0/16.0)*pow(pabc,3) + (5.0/128.0)*pow(pabc,4)) / pa + s00;
            }
            z_t1 = ray_position[2] + ds_t1 * direction[2];
            if(ds_t1 <= 0.0 || z_t1*this -> theta1_sgn < 0.0){
                val_t1 = false;
            }
        }
        else{
            if(fabs(pabc) > 1e-4 || pb > 0.0){
                ds_t1 = -0.5 * (pb + sdet) / pa + s00;
            }
            else{
                ds_t1 = -0.5 * pb * (0.5*pabc + (1.0/8.0)*pabc*pabc + (1.0/16.0)*pow(pabc,3) + (5.0/128.0)*pow(pabc,4)) / pa + s00;
            }
            z_t1 = ray_position[2] + ds_t1 * direction[2];
            if(ds_t1 <= 0.0 || z_t1*this -> theta1_sgn < 0.0){
                if(fabs(pabc) > 1e-4 || pb < 0.0){
                    ds_t1 = -0.5 * (pb + sdet) / pa + s00;
                }
                else{
                    ds_t1 = -0.5 * pb * (0.5*pabc + (1.0/8.0)*pabc*pabc + (1.0/16.0)*pow(pabc,3) + (5.0/128.0)*pow(pabc,4)) / pa + s00;
                }
                z_t1 = ray_position[2] + ds_t1 * direction[2];
                if(ds_t1 <= 0.0 || z_t1*this -> theta1_sgn < 0.0) {
                    val_t1 = false;
                }
            }
        }
    }
    double oneminus = 0.99999999999999;
    double oneplus = 1.00000000000001;
    if(val_t1){
        x_t1 = ray_position[0] + ds_t1 * direction[0];
        y_t1 = ray_position[1] + ds_t1 * direction[1];
        z_t1 = ray_position[2] + ds_t1 * direction[2];

        theta_t1 = this -> y_points[0];

        r_t1 = sqrt(x_t1*x_t1 + y_t1*y_t1 + z_t1*z_t1);

        if(r_t1 < this -> x_points[0]*oneminus || r_t1 > this -> x_points[this -> number_of_points_x-1]*oneplus){
            val_t1 = false;
        }

        phi_t1 = common::convert_cartesian_to_spherical_coordinates(x_t1,y_t1,z_t1)[2];

        if(not(this -> amr_cyclic_xyz[2])) {
            if (phi_t1 < this->z_points[0] - tiny) {
                val_t1 = false;
            }
            if (phi_t1 > this->z_points[this->number_of_points_z-1] + tiny) {
                val_t1 = false;
            }
        }
    }
}

void spherical_regular_grid::cross_outer_theta_cone(double &x_t2, double &y_t2, double &z_t2,bool &val_t2, double &ds_t2,double s00, double x00, double y00, double z00, const std::vector<double> &ray_position, const std::vector<double> &direction,double &r_t2,double &theta_t2,double &phi_t2){
    double st2,ct2,pa,pb,pc,pabc,det,sdet;
    st2 = this -> amrray_sint2*this -> amrray_sint2;
    ct2 = this -> amrray_cost2*this -> amrray_cost2;
    pa = ct2 * (direction[0]*direction[0] + direction[1]*direction[1]) - st2 * direction[2]*direction[2];
    pb = 2*ct2 * (direction[0]*x00 + direction[1]*y00) - 2*st2*direction[2]*z00;
    pc = ct2 * (x00*x00 + y00*y00) - st2*z00*z00;
    det = pb * pb - 4 * pa*pc;
    if(det <= 0.0){
        ds_t2 = -1e99;
        val_t2 = false;
    }
    else{
        sdet = sqrt(det);
        pabc = 4*pa*pc/(pb*pb);
        if(this -> y_points[this -> number_of_points_y-1] > common::get_pi_half() && ray_position[0]*direction[0] + ray_position[1]*direction[1] <= 0.0){
            if(fabs(pabc) > 1e-4 || pb < 0.0){
                ds_t2 = -0.5 * (pb - sdet) / pa + s00;
            }
            else{
                ds_t2 = -0.5 * pb * (0.5*pabc + (1.0/8.0)*pabc*pabc + (1.0/16.0)*pow(pabc,3) + (5.0/128.0)*pow(pabc,4)) / 4.0 + s00;
            }
            z_t2 = ray_position[2] + ds_t2 * direction[2];
            if((ds_t2 <= 0.0) || (z_t2 * theta2_sgn < 0.0)){
                val_t2 = false;
            }
        }
        else{
            if(fabs(pabc) > 1e-4 || pb > 0.0){
                ds_t2 = -0.5 * (pb + sdet ) / pa + s00;
            }
            else{
                ds_t2 = -0.5 * pb * (0.5*pabc + (1.0/8.0)*pabc*pabc + (1.0/16.0)*pow(pabc,3) + (5.0/128.0)*pow(pabc,4)) / 4.0 + s00;
            }
            z_t2 = ray_position[2] + ds_t2 * direction[2];
            if((ds_t2 <= 0.0) || (z_t2 * theta2_sgn < 0.0)){
                val_t2 = false;
            }
        }
    }
    double oneminus = 0.99999999999999;
    double oneplus = 1.00000000000001;
    if(val_t2){
        x_t2 = ray_position[0] + ds_t2 * direction[0];
        y_t2 = ray_position[1] + ds_t2 * direction[1];
        z_t2 = ray_position[2] + ds_t2 * direction[2];

        theta_t2 = this -> y_points[this -> number_of_points_y-1];

        r_t2 = sqrt(x_t2*x_t2 + y_t2*y_t2 + z_t2*z_t2);

        if(r_t2 < this -> x_points[0]*oneminus || r_t2 > this -> x_points[this->number_of_points_x-1]*oneplus){
            val_t2 = false;
        }
        phi_t2 = common::convert_cartesian_to_spherical_coordinates(x_t2,y_t2,z_t2)[2];

        if(not(this -> amr_cyclic_xyz[2])) {
            if (phi_t2 < this->z_points[0] - tiny) {
                val_t2 = false;
            }
            if (phi_t2 > this->z_points[this->number_of_points_z-1] + tiny) {
                val_t2 = false;
            }
        }
    }
}

void spherical_regular_grid::cross_inner_phi_cone(double &x_p1, double &y_p1, double &z_p1,bool &val_p1, double &ds_p1,const std::vector<double> &ray_position, const std::vector<double> &direction,double &r_p1,double &theta_p1,double &phi_p1){
    double dum1,dum2,err;
    dum1 = ray_position[0]*this -> amrray_sinp1 - ray_position[1]*this -> amrray_cosp1;
    dum2 = direction[0]*this -> amrray_cosp1 - direction[1]*this -> amrray_sinp1;
    if(dum2 == 0.0){
        val_p1 = false;
        ds_p1 = -1e99;
    }
    else{
        ds_p1 = dum1/dum2;
    }
    if(ds_p1 <= 0.0){
        val_p1 = false;
    }
    x_p1 = ray_position[0] + ds_p1 * direction[0];
    y_p1 = ray_position[1] + ds_p1 * direction[1];
    z_p1 = ray_position[2] + ds_p1 * direction[2];

    r_p1 = sqrt(x_p1*x_p1 + y_p1*y_p1 + z_p1*z_p1);
    theta_p1 = acos(z_p1/r_p1);
    phi_p1 = common::convert_cartesian_to_spherical_coordinates(x_p1,y_p1,z_p1)[2];

    if((phi_p1 < 0.0) || (phi_p1 > common::get_two_pi())){
        exit(0);
    }

    err = fabs(phi_p1 - this -> z_points[0]);

    if(err < 1e-4) {
        phi_p1 = this->z_points[0];
    }
    else{
        err = fabs(phi_p1 - common::get_two_pi() - this -> z_points[0]);
        if(err < 1e-4){
            phi_p1 = this -> z_points[0];
        }
        else{
            val_p1 = false;
        }
    }

    double oneminus = 0.99999999999999;
    double oneplus = 1.00000000000001;

    if((theta_p1 <= this -> y_points[0] - this -> tiny) || (theta_p1 >= this -> y_points[this -> number_of_points_y-1] + this -> tiny)){
        val_p1 = false;
    }
}

void spherical_regular_grid::cross_outer_phi_cone(double &x_p2, double &y_p2, double &z_p2,bool &val_p2, double &ds_p2,const std::vector<double> &ray_position, const std::vector<double> &direction,double &r_p2,double &theta_p2,double &phi_p2){
    double dum1,dum2,err;
    dum1 = ray_position[0]*this -> amrray_sinp2 - ray_position[1]*this -> amrray_cosp2;
    dum2 = direction[0]*this -> amrray_cosp2 - direction[1]*this -> amrray_sinp2;
    if(dum2 == 0.0){
        val_p2 = false;
        ds_p2 = -1e99;
    }
    else{
        ds_p2 = dum1 / dum2;
    }
    if(ds_p2 <= 0.0){
        val_p2 = false;
    }
    x_p2 = ray_position[0] + ds_p2 * direction[0];
    y_p2 = ray_position[1] + ds_p2 * direction[1];
    z_p2 = ray_position[2] + ds_p2 * direction[2];

    r_p2 = sqrt(x_p2 * x_p2 + y_p2 * y_p2 + z_p2 * z_p2);
    theta_p2 = acos(z_p2 / r_p2);
    phi_p2 = common::convert_cartesian_to_spherical_coordinates(x_p2, y_p2, z_p2)[2];

    if((phi_p2 < 0.0) || (phi_p2 > common::get_two_pi())){
        exit(0);
    }

    err = fabs(phi_p2 - this -> z_points[0]);

    if(err < 1e-4) {
        phi_p2 = this->z_points[0];
    }
    else{
        err = fabs(phi_p2 - common::get_two_pi() - this -> z_points[0]);
        if(err < 1e-4){
            phi_p2 = this -> z_points[0];
        }
        else{
            val_p2 = false;
        }
    }

    double oneminus = 0.99999999999999;
    double oneplus = 1.00000000000001;

    if((theta_p2 <= this -> y_points[0] - this -> tiny) || (theta_p2 >= this -> y_points[this -> number_of_points_y-1] + this -> tiny)){
        val_p2 = false;
    }
}

void spherical_regular_grid::entering_from_radius(photon &photon_i, double dummy, int ixi, std::vector<double> &ray_position, double x_r, double y_r, double z_r, double r_r, double theta_r, double phi_r){
    int ix,iy,iz;
    ray_position[0] = x_r;
    ray_position[1] = y_r;
    ray_position[2] = z_r;
    // This only when entering from inner radius
    // The calculation of the crossing with the inner radius may suffer
    // from linear precision loss if the starting point of the ray
    // is several factors of 10 further out. While the s00 trick (see
    // above) prevents quadratic loss of precision, which would easily
    // lead to catastrophic problems, it cannot avoid linear loss of
    // precision. While linear loss of precision is not so damaging,
    // it might in very special circumstances cause problems if you have
    // very finely space gridding in radius, for instance. So let us
    // correct for this. Note that we already checked beforehand that
    // the error is not too large.

    ray_position[0] = dummy * ray_position[0];
    ray_position[1] = dummy * ray_position[1];
    ray_position[2] = dummy * ray_position[2];

    iy = common::hunt(this->y_points, this->number_of_points_y-1, theta_r, this->number_of_points_y-1);
    iz = common::hunt(this->z_points, this->number_of_points_z-1, phi_r, this->number_of_points_z-1);
    ix = ixi;
    if (iy > this->number_of_points_y - 1) {
        iy = number_of_points_y - 1;
    }
    if (iy < 0) {
        iy = 0;
    }
    if (iz > this->number_of_points_z - 1) {
        if (this->amr_cyclic_xyz[2]) {
            iz = 0;
        } else {
            iz = this->number_of_points_z - 1;
        }
    }
    if (iz < 0) {
        if(this -> amr_cyclic_xyz[2]){
            iz = this -> number_of_points_z - 1;
        }
        else{
            iz = 0;
        }
    }
    std::vector<int> grid_position = {ix,iy,iz};
    photon_i.set_grid_position(grid_position);
    photon_i.set_ray_position(ray_position);
}

void spherical_regular_grid::entering_from_theta(photon &photon_i, int iyi, std::vector<double> &ray_position,double x_t,double y_t,double z_t,double r_t,double phi_t){
    int ix,iy,iz;
    ray_position[0] = x_t;
    ray_position[1] = y_t;
    ray_position[2] = z_t;

    ix = common::hunt(this -> x_points,this -> number_of_points_x-1, r_t, this -> number_of_points_x-1);
    iz = common::hunt(this -> z_points,this -> number_of_points_z-1, phi_t, this -> number_of_points_z-1);
    iy = iyi;
    if(ix < 0){
        ix = 0;
    }
    if(ix > this -> number_of_points_x - 1){
        ix = this -> number_of_points_x - 1;
    }
    if (iz > this->number_of_points_z - 1) {
        if (this->amr_cyclic_xyz[2]) {
            iz = 0;
        } else {
            iz = this->number_of_points_z - 1;
        }
    }
    if (iz < 0) {
        if(this -> amr_cyclic_xyz[2]){
            iz = this -> number_of_points_z - 1;
        }
        else{
            iz = 0;
        }
    }
    std::vector<int> grid_position = {ix,iy,iz};
    photon_i.set_grid_position(grid_position);
    photon_i.set_ray_position(ray_position);
}

void spherical_regular_grid::entering_from_phi(photon &photon_i, int izi, std::vector<double> &ray_position, double x_p, double y_p, double z_p, double r_p, double theta_p){
    int ix,iy,iz;
    ray_position[0] = x_p;
    ray_position[1] = y_p;
    ray_position[2] = z_p;

    ix = common::hunt(this -> x_points, this -> number_of_points_x-1,r_p,this -> number_of_points_x-1);
    iy = common::hunt(this -> y_points, this -> number_of_points_y-1,theta_p,this -> number_of_points_y-1);
    iz = izi;
    if(ix < 0){
        ix = 0;
    }
    if(ix > this -> number_of_points_x - 1){
        ix = this -> number_of_points_x - 1;
    }
    if(iy < 0){
        iy = 0;
    }
    if(iy > this -> number_of_points_y - 1){
        iy = this -> number_of_points_y - 1;
    }
    std::vector<int> grid_position = {ix,iy,iz};
    photon_i.set_grid_position(grid_position);
    photon_i.set_ray_position(ray_position);
}
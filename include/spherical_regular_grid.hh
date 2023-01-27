#ifndef RADMC_SPHERICAL_REGULAR_GRID_HH
#define RADMC_SPHERICAL_REGULAR_GRID_HH

#include <regular_grid.hh>
#include <common.hh>

class spherical_regular_grid: public regular_grid {

private:

    double tiny = 1e-14;
    double amrray_sint1;
    double amrray_sint2;
    double amrray_cost1;
    double amrray_cost2;
    double amrray_sinp1;
    double amrray_sinp2;
    double amrray_cosp1;
    double amrray_cosp2;
    double theta1_sgn = 1;
    double theta2_sgn = -1;
    std::vector<std::vector<std::vector<double>>> cell_volume;
    std::vector<bool> amr_cyclic_xyz;

public:

    spherical_regular_grid();

    void initialize_grid();

    double get_cell_volume(std::vector<int> grid_position) const;

    void adjust_theta();

    void adjust_phi();

    void calculate_cell_volume();

    std::vector<int> found_ray_position_in_grid(const std::vector<double>& ray_position);

    void calculate_points_delta();

    void calculate_photon_cell_walls(photon &photon_i);

    void calculate_photon_new_position(photon &photon_i);

    void move_photon_inside(photon &photon_i);

    //Metodos para mover el fotón afuera
    void move_photon_outside(photon &photon_i);

    void
    entering_from_radius(photon &photon_i, double dummy, int ixi, std::vector<double> &ray_position, double x_r, double y_r, double z_r,
                         double r_r, double theta_r, double phi_r);

    void entering_from_theta(photon &photon_i, int iyi, std::vector<double> &ray_position, double x_t, double y_t, double z_t, double r_t,
                             double phi_t);

    void entering_from_phi(photon &photon_i, int izi, std::vector<double> &ray_position, double x_p, double y_p, double z_p, double r_p,
                           double theta_p);

    void
    cross_inner_radius(double &x_r1, double &y_r1, double &z_r1, bool &val_r1, double &ds_r1, double x00dotdir,
                       double r002,
                       double s00, const std::vector<double> &ray_position, const std::vector<double> &direction,
                       double &r_r1, double &theta_r1, double &phi_r1);

    void
    cross_outer_radius(double &x_r2, double &y_r2, double &z_r2, bool &val_r2, double &ds_r2, double x00dotdir,
                       double r002,
                       double r02, double s00, const std::vector<double> &ray_position,
                       const std::vector<double> &direction, double &r_r2, double &theta_r2, double &phi_r2);

    void cross_inner_theta_cone(double &x_t1, double &y_t1, double &z_t1, bool &val_t1, double &ds_t1,
                                const std::vector<double> &ray_position, const std::vector<double> &direction,
                                double x00,
                                double y00, double z00, double s00, double &r_t1, double &theta_t1, double &phi_t1);

    void
    cross_outer_theta_cone(double &x_t2, double &y_t2, double &z_t2, bool &val_t2, double &ds_t2, double s00,
                           double x00,
                           double y00, double z00, const std::vector<double> &ray_position,
                           const std::vector<double> &direction, double &r_t2, double &theta_t2, double &phi_t2);

    void cross_inner_phi_cone(double &x_p1, double &y_p1, double &z_p1, bool &val_p1, double &ds_p1,
                              const std::vector<double> &ray_position, const std::vector<double> &direction,
                              double &r_p1,
                              double &theta_p1, double &phi_p1);

    void cross_outer_phi_cone(double &x_p2, double &y_p2, double &z_p2, bool &val_p2, double &ds_p2,
                              const std::vector<double> &ray_position, const std::vector<double> &direction,
                              double &r_p2,
                              double &theta_p2, double &phi_p2);

    ~spherical_regular_grid();
};


#endif //RADMC_SPHERICAL_REGULAR_GRID_HH

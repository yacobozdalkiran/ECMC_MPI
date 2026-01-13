//
// Created by ozdalkiran-l on 1/8/26.
//

#include "observables.h"

//Returns the value of the mean plaquette of the gauge field
double observables::mean_plaquette(const GaugeField &field, const Geometry &geo) {
    double sum = 0.0;
    size_t counter = 0;
    SU3 staple;
    for (size_t site = 0; site < geo.V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            compute_staple(field, geo, site, mu, staple);
            double p = (field.view_link_const(site, mu) * staple).trace().real() / 18.0; //N_plaquettes*N_colors
            sum += p;
            counter += 1;
        }
    }
    double mean = sum / static_cast<double>(counter);
    return mean;
}

//Returns the value of the mean plaquette of the gauge field ignoring the frozen sites plaquettes
double observables::mean_plaquette(const GaugeField &field, const mpi::GeometryFrozen &geo) {
    double sum = 0.0;
    size_t counter = 0;
    SU3 staple;
    for (int t = 1; t<geo.L-1; t++){
        for (int z = 1; z<geo.L-1; z++) {
            for (int y = 1; y<geo.L-1; y++) {
                for (int x = 1; x<geo.L-1; x++) {
                    size_t site = geo.index(x,y,z,t);
                    for (int mu = 0; mu < 4; mu++) {
                        compute_staple(field, geo, site, mu, staple);
                        double p = (field.view_link_const(site, mu) * staple).trace().real() / 18.0; //N_plaquettes*N_colors
                        sum += p;
                        counter += 1;
                    }
                }
            }
        }
    }
    double mean = sum / static_cast<double>(counter);
    return mean;
}

//Returns the value of the Wilson action of the gauge field
double observables::wilson_action(const GaugeField &field, const Geometry &geo) {
    double action = 0.0;
    SU3 U0, U1, U2, U3;
    SU3 id = SU3::Identity();
    SU3 plaquette;
    for (size_t site = 0; site < geo.V; site++) {
        for (int nu =0; nu<4; nu++) {
            for (int mu = 0; mu<nu; mu++) {
                U0 = field.view_link_const(site, mu);
                U1 = field.view_link_const(geo.get_neigh(site,mu,0), nu);
                U2 = field.view_link_const(geo.get_neigh(site,nu,0), mu).adjoint();
                U3 = field.view_link_const(site, nu).adjoint();
                plaquette = U0*U1*U2*U3;
                action += (id - plaquette).trace().real();
            }
        }
    }
    action *= (1.0/3.0); //1/Ncolors
    //action /= (lat.V*4.0);
    return action;
}

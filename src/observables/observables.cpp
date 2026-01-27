//
// Created by ozdalkiran-l on 1/16/26.
//

#include "observables.h"


//Returns the value of the mean plaquette of the gauge field
double observables::mean_plaquette(const GaugeField &field, const Geometry &geo) {
    double sum = 0.0;
    size_t counter = 0;
    SU3 U1, U2, U3, U4;
    for (size_t site = 0; site < geo.V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu+1; nu<4; nu++) {
                U1 = field.view_link_const(site, mu);
                U2 = field.view_link_const(geo.get_neigh(site,mu,0), nu);
                U3 = field.view_link_const(geo.get_neigh(site,nu,0), mu).adjoint();
                U4 = field.view_link_const(site, nu).adjoint();
                sum += (U1*U2*U3*U4).trace().real()/3.0;
                counter++;
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
    action /= 3.0; //Ncolors
    //action /= (lat.V*4.0);
    return action;
}

//Applies a random gauge transform to the gauge field
void observables::gauge_transform(GaugeField &field, const Geometry &geo, std::mt19937_64 rng) {
    GaugeField gauge_field(geo);
    gauge_field.hot_start(rng);
    //We will use only the mu=0 links of the gauge_transf field as the transforming gauge field
    for (int t = 0; t<geo.T; t++) {
        for (int z=0; z<geo.L; z++) {
            for (int y=0; y<geo.L; y++) {
                for (int x=0; x<geo.L; x++) {
                    size_t site = geo.index(x,y,z,t);
                    for (int mu=0; mu<4; mu++) {
                        size_t npmu = geo.get_neigh(site, mu, 0);
                        SU3 Uold = field.view_link_const(site, mu);
                        field.view_link(site, mu) = gauge_field.view_link_const(site, 0)*Uold*gauge_field.view_link_const(npmu, 0).adjoint();
                    }
                }
            }
        }
    }
}

//Maximal drift from unitarity in field
double observables::max_drift_det(const GaugeField &field, const Geometry &geo) {
    double res = 0.0;
    for (size_t site = 0; site<geo.V; site++) {
        for (int mu=0; mu<4; mu++) {
            double drift = std::abs(field.view_link_const(site, mu).determinant()-1.0);
            if (drift>res) res=drift;
        }
    }
    return res;
}

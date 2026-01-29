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

//Computes G_{\mu\nu}(site) clover
SU3 observables::clover_site(const GaugeField &field, const Geometry &geo, size_t site, int mu, int nu) {
    if (mu==nu) std::cerr << "mu = nu => G(site, mu, nu) = 0\n";
    size_t x = site;
    size_t xpmu = geo.get_neigh(x,mu,0); //x+mu
    size_t xpnu = geo.get_neigh(x,nu,0); //x+nu
    size_t xmmu = geo.get_neigh(x,mu,1); //x-mu
    size_t xmnu = geo.get_neigh(x,nu,1); //x-nu
    size_t xmmupnu = geo.get_neigh(xmmu,nu,0); //x-mu+nu
    size_t xmmumnu = geo.get_neigh(xmmu,nu,1); //x-mu-nu
    size_t xpmumnu = geo.get_neigh(xpmu,nu,1); //x+mu-nu
    SU3 clover = SU3::Zero();
    clover += field.view_link_const(x, mu) * field.view_link_const(xpmu, nu) * field.view_link_const(xpnu, mu).adjoint() * field.view_link_const(x, nu).adjoint();
    clover += field.view_link_const(x, nu) * field.view_link_const(xmmupnu, mu).adjoint() * field.view_link_const(xmmu,nu).adjoint() * field.view_link_const(xmmu,mu);
    clover += field.view_link_const(xmmu, mu).adjoint() * field.view_link_const(xmmumnu, nu).adjoint() * field.view_link_const(xmmumnu, mu) * field.view_link_const(xmnu, nu);
    clover += field.view_link_const(xmnu, nu).adjoint() * field.view_link_const(xmnu, mu) * field.view_link_const(xpmumnu, nu) * field.view_link_const(x, mu).adjoint();
    SU3 F = (clover - clover.adjoint()).eval();
    F *= 0.125;
    return F;
}

//Computes the local clover topological charge and energy at site
std::pair<double,double> observables::local_q_e_clover(const GaugeField &field, const Geometry &geo, size_t site) {
    SU3 F01 = clover_site(field, geo, site, 0, 1);
    SU3 F02 = clover_site(field, geo, site, 0, 2);
    SU3 F03 = clover_site(field, geo, site, 0, 3);
    SU3 F12 = clover_site(field, geo, site, 1, 2);
    SU3 F13 = clover_site(field, geo, site, 1, 3);
    SU3 F23 = clover_site(field, geo, site, 2, 3);

    // 2. On utilise la forme explicite de epsilon_{mu nu rho sigma}
    // Q ~ Tr(F01*F23 - F02*F13 + F03*F12)
    double q = (F01 * F23).trace().real() -
               (F02 * F13).trace().real() +
               (F03 * F12).trace().real();

    //On calcule aussi l'énergie locale
    double e_local = 0.5*(F01.squaredNorm() + F02.squaredNorm() + F03.squaredNorm() +
                  F12.squaredNorm() + F13.squaredNorm() + F23.squaredNorm());
    // 3. Le facteur global est 1/(4*pi^2) car le 1/32 a été absorbé
    // par les combinaisons et le facteur 1/8 de F.
    return {q * (1.0 / (4.0 * M_PI * M_PI)),e_local};
}

//Return the clover topological charge and energy density
std::pair<double,double> observables::topo_q_e_clover(const GaugeField &field, const Geometry &geo) {
    double q = 0.0;
    double e = 0.0;
    for (size_t site = 0; site < geo.V; site++) {
        auto local = local_q_e_clover(field, geo, site);
        q += local.first;
        e += local.second;
    }
    return {q, e/static_cast<double>(geo.V)};
}


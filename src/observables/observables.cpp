//
// Created by ozdalkiran-l on 1/8/26.
//

#include "observables.h"
#include "../mpi/HalosExchange.h"

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

//Returns the value of the mean plaquette of the gauge field ignoring the frozen sites plaquettes
double observables::mean_plaquette(const GaugeField &field, const mpi::GeometryFrozen &geo) {
    double sum = 0.0;
    size_t counter = 0;
    SU3 U1, U2, U3, U4;
    for (int t = 1; t<geo.L-1; t++){
        for (int z = 1; z<geo.L-1; z++) {
            for (int y = 1; y<geo.L-1; y++) {
                for (int x = 1; x<geo.L-1; x++) {
                    size_t site = geo.index(x,y,z,t);
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

//Returns the sum of retr/3 of plaquettes of the bulk (i.e. coords btw 1 and L-2)
double mpi::observables::sum_plaquette_bulk(const GaugeField &field, const mpi::GeometryFrozen &geo) {
    double sum = 0.0;
    SU3 U1, U2, U3, U4;
    for (int t = 1; t<geo.L-1; t++){
        for (int z = 1; z<geo.L-1; z++) {
            for (int y = 1; y<geo.L-1; y++) {
                for (int x = 1; x<geo.L-1; x++) {
                    size_t site = geo.index(x,y,z,t);
                    for (int mu = 0; mu < 4; mu++) {
                        for (int nu = mu+1; nu<4; nu++) {
                            U1 = field.view_link_const(site, mu);
                            U2 = field.view_link_const(geo.get_neigh(site,mu,0), nu);
                            U3 = field.view_link_const(geo.get_neigh(site,nu,0), mu).adjoint();
                            U4 = field.view_link_const(site, nu).adjoint();
                            sum += (U1*U2*U3*U4).trace().real()/3.0;
                        }
                    }
                }
            }
        }
    }
    return sum;
}

//Returns the sum of retr/3 of plaquettes on the boundaries of the field using the halos
double mpi::observables::sum_plaquette_boundaries(const GaugeField &field, const mpi::GeometryFrozen &geo,
    const HaloObs &halo_obs) {
    double sum_boundaries = 0.0;
    int L = geo.L;
    SU3 U1, U2, U3, U4;
    for (int t=0; t<L; t++) {
        for (int z=0; z<L; z++) {
            for (int y=0; y<L; y++) {
                for (int x=0; x<L; x++) {
                    //We treat only the boundary sites
                    if (x==0 || y==0 || z==0 || t==0 || x==L-1 || y==L-1 || z==L-1 || t==L-1) {
                        for (int mu=0; mu<4; mu++) {
                            for (int nu=mu+1; nu<4; nu++) {
                                U1 = halo_obs.get_link_with_halo_obs(field, geo, x,y,z,t,mu);
                                U2 = halo_obs.get_link_with_halo_obs(field, geo, x+(mu==0), y+(mu==1), z+(mu==2), t+(mu==3), nu);
                                U3 = halo_obs.get_link_with_halo_obs(field, geo, x+(nu==0), y+(nu==1), z+(nu==2), t+(nu==3), mu).adjoint();
                                U4 = halo_obs.get_link_with_halo_obs(field, geo, x, y, z, t, nu).adjoint();
                                sum_boundaries += (U1*U2*U3*U4).trace().real()/3.0;
                            }
                        }
                    }
                }
            }
        }
    }
    return sum_boundaries;
}

//Returns the sum of retr of plaquettes of the local field (includes halo filling and exchange)
double mpi::observables::mean_plaquette_local(const GaugeField &field, const mpi::GeometryFrozen &geo,
                                              HaloObs &halo_obs, mpi::MpiTopology &topo) {
    double local_mean_plaquette = 0.0;
    MPI_Request reqs[16];
    mpi::observables::fill_halo_obs_send(field, geo, halo_obs);
    mpi::observables::exchange_halos_obs(halo_obs, topo, reqs);
    //During the exchanges, we compute the plaquette in the bulk
    local_mean_plaquette += mpi::observables::sum_plaquette_bulk(field, geo);
    //Now we wait for all the nodes to end the exchanges
    MPI_Waitall(16, reqs, MPI_STATUSES_IGNORE);
    //Finally we compute the mean plaquette of the boundary sites
    local_mean_plaquette += sum_plaquette_boundaries(field, geo, halo_obs);
    return local_mean_plaquette;
}

//Returns the mean plaquette of the global lattice
double mpi::observables::mean_plaquette_global(const GaugeField &field, const mpi::GeometryFrozen &geo,
    HaloObs &halo_obs, mpi::MpiTopology &topo) {
    double local_mean_plaquette = mean_plaquette_local(field, geo, halo_obs, topo);
    double global_mean_plaquette = 0.0;
    MPI_Allreduce(&local_mean_plaquette, &global_mean_plaquette, 1, MPI_DOUBLE, MPI_SUM, topo.cart_comm);
    global_mean_plaquette /= 6.0*geo.V*topo.size;
    return global_mean_plaquette;
}

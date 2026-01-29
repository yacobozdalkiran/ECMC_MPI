//
// Created by ozdalkiran-l on 1/29/26.
//

#include "gradient_flow.h"

#include "../su3/utils.h"

GradientFlow::GradientFlow(double epsilon_, const GaugeField &field, const Geometry &geo) :
epsilon(epsilon_), field_c(field), force_0(field), force_1(field),geo_p(&geo){}

//Updates the force field according to staples and links of field_c, into force_0 or force_1 depending on i
void GradientFlow::compute_force(int i) {
    SU3 staple{};
    SU3 tmp{};
    for (size_t site=0; site<geo_p->V; site++) {
        for (int mu=0; mu<4; mu++) {
            field_c.compute_staple(*geo_p, site, mu, staple);
            tmp = staple * field_c.view_link_const(site, mu).adjoint();
            proj_lie_su3(tmp);
            if (i == 0) force_0.view_link(site, mu) = tmp;
            if (i == 1) force_1.view_link(site, mu) = tmp;
        }
    }
}

//First RK3 step
void GradientFlow::compute_w1() {
    SU3 tmp{};
    for (size_t site=0; site<geo_p->V; site++) {
        for (int mu=0; mu<4; mu++) {
            tmp = exp_analytic(force_0.view_link_const(site, mu), 0.25*epsilon)*field_c.view_link_const(site, mu);
            field_c.view_link(site, mu) = tmp;
        }
    }
}

//Second RK3 step
void GradientFlow::compute_w2() {
    SU3 tmp{};
    SU3 Z{};
    for (size_t site=0; site<geo_p->V; site++) {
        for (int mu=0; mu<4; mu++) {
            Z = (8.0/9.0)*epsilon*force_1.view_link_const(site, mu) - (17.0/36.0)*epsilon*force_0.view_link_const(site, mu);
            tmp = exp_analytic(Z, 1.0)*field_c.view_link_const(site, mu);
            field_c.view_link(site, mu) = tmp;
        }
    }
}

//Third RK3 step
void GradientFlow::compute_w3() {
    SU3 tmp{};
    SU3 Z{};
    for (size_t site=0; site<geo_p->V; site++) {
        for (int mu=0; mu<4; mu++) {
            Z = 0.75*epsilon*force_0.view_link_const(site, mu) - (15.0/12.0)*epsilon*force_1.view_link_const(site, mu);
            tmp = exp_analytic(Z, 1.0)*field_c.view_link_const(site, mu);
            field_c.view_link(site, mu) = tmp;
        }
    }
}

//Performs a full RK3 step, the field updated of epsilon is in field_c
void GradientFlow::rk3_step() {
    compute_force(0);
    compute_w1();
    compute_force(1);
    compute_w2();
    compute_force(0);
    compute_w3();
}

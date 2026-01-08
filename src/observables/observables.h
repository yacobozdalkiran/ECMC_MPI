//
// Created by ozdalkiran-l on 1/8/26.
//

#ifndef INC_4D_MPI_OBSERVABLES_H
#define INC_4D_MPI_OBSERVABLES_H

#include "../gauge/GaugeField.h"

namespace observables {
    double mean_plaquette(const GaugeField &field, const Geometry &geo);
    double wilson_action(const GaugeField &field, const Geometry &geo);
}

#endif //INC_4D_MPI_OBSERVABLES_H
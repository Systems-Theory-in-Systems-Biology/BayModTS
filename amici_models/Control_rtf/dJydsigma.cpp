#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_Control_rtf {

void dJydsigma_Control_rtf(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_y_obs - 1.0*std::pow(-my_obs + y_obs, 2)/std::pow(sigma_y_obs, 3);
            break;
    }
}

} // namespace model_Control_rtf
} // namespace amici

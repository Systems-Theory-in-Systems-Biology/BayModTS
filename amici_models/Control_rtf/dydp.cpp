#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_Control_rtf {

void dydp_Control_rtf(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 0:
            dydp[0] = 1 - std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t1);
            break;
        case 1:
            dydp[0] = (1 - std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t11))*std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t2);
            break;
        case 2:
            dydp[0] = Asus*(std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)*std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t1)/std::pow(t1, 2);
            break;
        case 3:
            dydp[0] = Atrans*(std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)*std::exp((1.0/t2 + 1.0/t11)*(std::log(std::pow(10, Tshift) + 1) - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange)))/M_LN10)/std::pow(t11, 2);
            break;
        case 4:
            dydp[0] = -Atrans*(1 - std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t11))*(std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)*std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t2)/std::pow(t2, 2);
            break;
    }
}

} // namespace model_Control_rtf
} // namespace amici

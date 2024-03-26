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

void w_Control_rtf(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl){
    amici_y = Asus*(1 - std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t1)) + Atrans*(1 - std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t11))*std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t2) + p0;  // w[0]
}

} // namespace model_Control_rtf
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Control_rtf {

static constexpr std::array<sunindextype, 7> dwdp_colptrs_Control_rtf_ = {
    0, 1, 2, 3, 4, 5, 5
};

void dwdp_colptrs_Control_rtf(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Control_rtf_));
}
} // namespace model_Control_rtf
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Control_rtf {

static constexpr std::array<sunindextype, 5> dwdp_rowvals_Control_rtf_ = {
    0, 0, 0, 0, 0
};

void dwdp_rowvals_Control_rtf(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Control_rtf_));
}
} // namespace model_Control_rtf
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "w.h"
#include "dwdp.h"

namespace amici {
namespace model_Control_rtf {

void dwdp_Control_rtf(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    damici_y_dAsus = 1 - std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t1);  // dwdp[0]
    damici_y_dAtrans = (1 - std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t11))*std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t2);  // dwdp[1]
    damici_y_dt1 = Asus*(std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)*std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t1)/std::pow(t1, 2);  // dwdp[2]
    damici_y_dt11 = Atrans*(std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)*std::exp((1.0/t2 + 1.0/t11)*(std::log(std::pow(10, Tshift) + 1) - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange)))/M_LN10)/std::pow(t11, 2);  // dwdp[3]
    damici_y_dt2 = -Atrans*(1 - std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t11))*(std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)*std::exp((std::log(std::pow(10, Tshift) + 1)/M_LN10 - std::log(std::pow(10, Tshift) + std::pow(10000000000, t/Trange))/M_LN10)/t2)/std::pow(t2, 2);  // dwdp[4]
}

} // namespace model_Control_rtf
} // namespace amici

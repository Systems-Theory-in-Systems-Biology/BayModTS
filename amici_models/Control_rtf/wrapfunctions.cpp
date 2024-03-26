#include "wrapfunctions.h"
#include "Control_rtf.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Control_rtf::Model_Control_rtf()
    );
}

} // namespace generic_model

} // namespace amici

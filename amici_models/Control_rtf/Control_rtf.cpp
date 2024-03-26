#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Control_rtf {

// clang-format off

std::array<const char*, 6> parameterNames = {
    "Amplitude sustained response", // p[0]
"Amplitude transient response", // p[1]
"Sustained response time constant", // p[2]
"Transient response time constant", // p[3]
"Transient decay time constant", // p[4]
"noiseParameter1_y_obs", // p[5]
};

std::array<const char*, 3> fixedParameterNames = {
    "Offset", // k[0]
"Response time", // k[1]
"Time span of the data", // k[2]
};

std::array<const char*, 0> stateNames = {
    
};

std::array<const char*, 1> observableNames = {
    "RTD output value", // y[0]
};

std::array<const ObservableScaling, 1> observableScalings = {
    ObservableScaling::lin, // y[0]
};

std::array<const char*, 2> expressionNames = {
    "y", // w[0]
"flux_r0", // w[1]
};

std::array<const char*, 6> parameterIds = {
    "Asus", // p[0]
"Atrans", // p[1]
"t1", // p[2]
"t11", // p[3]
"t2", // p[4]
"noiseParameter1_y_obs", // p[5]
};

std::array<const char*, 3> fixedParameterIds = {
    "p0", // k[0]
"Tshift", // k[1]
"Trange", // k[2]
};

std::array<const char*, 0> stateIds = {
    
};

std::array<const char*, 1> observableIds = {
    "y_obs", // y[0]
};

std::array<const char*, 2> expressionIds = {
    "amici_y", // w[0]
"flux_r0", // w[1]
};

std::array<int, 0> stateIdxsSolver = {
    
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Control_rtf

} // namespace amici

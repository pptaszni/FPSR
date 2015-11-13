#ifndef ENDOGENOUS
#define ENDOGENOUS
#include <vector>
#include <string>
#include "EquationCommon.hpp"

class EndogenousMethod
{
public:
    EndogenousMethod() {}
    ~EndogenousMethod() {}
    void Start();
    void SolveSampleEquation();
    void SolveSampleMatrixEquation();
    void SaveResults(const std::vector<state_type> out_states,
        const std::vector<time_type> out_time,
        std::string filename);
    void SaveResults(const std::vector<matrix_state_type> out_states,
        const std::vector<time_type> out_time,
        std::string filename);
};

#endif // ENDOGENOUS

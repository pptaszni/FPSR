#ifndef ENDOGENOUS
#define ENDOGENOUS
#include <vector>
#include "EquationCommon.hpp"

class EndogenousMethod
{
public:
    EndogenousMethod() {}
    ~EndogenousMethod() {}
    void Start();
    void SolveSampleEquation();
    void SaveResults(const std::vector<state_type> out_states,
        const std::vector<time_type> out_time);
};

#endif // ENDOGENOUS

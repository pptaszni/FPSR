#ifndef EQUATIONCOMMON
#define EQUATIONCOMMON
#include <vector>

typedef std::vector<double> state_type;
typedef double time_type;

struct SolutionObserver
{
    std::vector<state_type>& out_states_;
    std::vector<double>& out_time_;

    SolutionObserver(std::vector<state_type> &states, std::vector<time_type> &time);

    void operator()(const state_type &x, time_type t);
};

#endif // EQUATIONCOMMON

// Minimal script used to calculate V & R for a range of dT given C0 and an Alloy.

#include <optional>
#include <fstream>
#include <tuple>
#include <string>
#include <cmath>
#include "alloy.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
#include "optimiser.h"

struct Result{bool hasDiverged{}; double f1{}; double f2{}; double V{}; double R{};};

Result VRSolver(double dT, double C0, const alloy::Alloy& A)
{
    double f1{}, f2{}, dV{}, dR{};
    double V{approx::getTipVelocity(dT, C0, A)};
    double R{approx::getTipRadius(dT, C0, A)};
    diff::Jacobian J{};

    for (int i{0}; i<1'000; ++i)
    {
        std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
        J = diff::calculateGrads<models::LGK>(V, R, dT, C0, A);
        std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
        if (std::isnan(dV) || std::isnan(dR)) // solver diverges
            return Result{true};
        V += dV;
        R += dR;
    }
    
    return Result{false, f1, f2, V, R};
}

int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outf{dataPath + "/result.csv"};
    outf << "diverged,dT,C0,V,R,f1,f2\n";
    
    double C0{5.0};
    for(double dT{0}; dT<=100; dT+=0.1)
    {
        Result result{VRSolver(dT, C0, alloy::SnAg)};
        if (result.hasDiverged)
            outf << "True\n";
        else
            outf << "False" << ',' << dT << ',' << C0 << ',' << result.V << ',' << result.R << ',' << result.f1 << ','
                << result.f2 << '\n';
    }
    return 0;
}

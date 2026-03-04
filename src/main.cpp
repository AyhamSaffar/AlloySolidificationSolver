#include <iostream>
#include <tuple>
#include "alloy.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
//TODO #include "logger.h"


int main()
{
    double V{1e-5}, R{1e-5}, dT{0.5}, C0{5.0};
    alloy::Alloy A{alloy::SnAg};

    std::tuple<double, double> LGKResiduals{models::LGK(V, R, dT, C0, A)};
    std::cout << "LGK f1 =" << std::get<0>(LGKResiduals) << ", LGK f2 =" << std::get<1>(LGKResiduals) << '\n';

    diff::Jacobian J{diff::calculateGrads<models::LGK>(V, R, dT, C0, A)};
    std::cout << J << '\n';

    return 0;
}

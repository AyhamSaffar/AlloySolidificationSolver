#include <catch2/catch_test_macros.hpp>
#include <tuple>
#include "alloy.h"
#include "differentials.h"
#include "models.h"

TEST_CASE("diffs given reasonable values", "[advanced]")
{
    REQUIRE(diff::calculateGrads<models::LGK>(0.0280516,-1.64116e-07,0.1,5, alloy::SnAg).df1dV == 0.0);
    REQUIRE(1+1==2);
}


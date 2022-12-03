    scalar CoNum = -great;
/*
    forAll(fluidRegions, regionI)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                fluidRegions[regionI],
                runTime,
                rhoFluid[regionI],
                phiFluid[regionI]
            ),
            CoNum
        );
    }
*/

	CoNum = max
        (
            compressibleCourantNo
            (
                mesh,
                runTime,
                rho,
                phi
            ),
            CoNum
        );

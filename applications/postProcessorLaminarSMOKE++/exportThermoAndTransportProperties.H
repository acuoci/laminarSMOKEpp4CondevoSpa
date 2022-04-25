if (export_thermo_and_transport_properties == true)
{
	Cp.write();
	Cpv.write();
	lambda.write();
	mu.write();
	alpha.write();
	MW.write();

	forAll(Y, i)
	{
		volScalarField& Dmixi = Dmix[i];
		Dmixi.write();
	}
}

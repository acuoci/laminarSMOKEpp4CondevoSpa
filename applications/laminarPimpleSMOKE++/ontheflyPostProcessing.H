// Radiative heat loss
{
	if (runTime.outputTime() == true)
	{
		radiation->Qloss(T,Qrad);
	}
}

// OpenSMOKE++ reacting mixture
{
	mixture.onthefly_postprocessing(mesh, rho);
}

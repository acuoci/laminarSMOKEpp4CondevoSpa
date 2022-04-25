if (write_mole_fractions == true)
{
	Info<< "Write mole fractions on file..." << endl;

	for (unsigned int i=0;i<mixture.thermodynamicsMap().NumberOfSpecies();i++)
		mixture.X()[i].write();
}


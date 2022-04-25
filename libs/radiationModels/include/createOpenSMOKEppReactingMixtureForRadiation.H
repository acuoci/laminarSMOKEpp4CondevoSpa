// Create OpenSMOKEppReactingMixtureForRadiation object

std::vector<double> W(mixture.ns());
std::vector<std::string> species_names(mixture.ns());

for (unsigned int i=0;i<species_names.size();i++)
{
	species_names[i] = mixture.species()[i];
	W[i] = mixture.mw_species(i).value();
}

autoPtr<Foam::OpenSMOKEppReactingMixtureForRadiation> myOpenSMOKEppReactingMixtureForRadiation(Foam::OpenSMOKEppReactingMixtureForRadiation::New(T,p, Cpv, Y, W, species_names));


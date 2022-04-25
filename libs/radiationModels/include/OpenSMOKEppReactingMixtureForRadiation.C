#include "OpenSMOKEppReactingMixtureForRadiation.H"

namespace Foam
{

	Foam::autoPtr<OpenSMOKEppReactingMixtureForRadiation> OpenSMOKEppReactingMixtureForRadiation::New ( const volScalarField& T, const volScalarField& p, const volScalarField& Cpv, const PtrList<volScalarField>& Y, const std::vector<double>& W, const std::vector<std::string>& species_names )
	{
	    IOobject radIO
	    (
		"OpenSMOKEppReactingMixtureForRadiation",
		T.time().constant(),
		T.mesh(),
		IOobject::NO_READ,
		IOobject::NO_WRITE,
		true
	    );

	    return autoPtr<OpenSMOKEppReactingMixtureForRadiation>(new OpenSMOKEppReactingMixtureForRadiation(T,p,Cpv,Y,W,species_names));
	}

	Foam::IOobject OpenSMOKEppReactingMixtureForRadiation::createIOobject ( const fvMesh& mesh ) const
	{
	    IOobject io
	    (
		"OpenSMOKEppReactingMixtureForRadiation",
		mesh.time().constant(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    );

	    // Check if field exists and can be read
	    if (io.typeHeaderOk<volScalarField>(true))
	    {
		io.readOpt() = IOobject::NO_READ;
		return io;
	    }
	    else
	    {
		io.readOpt() = IOobject::NO_READ;
		return io;
	    }
	}


	OpenSMOKEppReactingMixtureForRadiation::OpenSMOKEppReactingMixtureForRadiation(const volScalarField& T, const volScalarField& p, const volScalarField& Cpv, const PtrList<volScalarField>& Y, const std::vector<double>& W, const std::vector<std::string>& species_names) :
		IOdictionary
		    (
			IOobject
			(
			    "OpenSMOKEppReactingMixtureForRadiation",
			    T.time().constant(),
			    T.mesh(),
			    IOobject::NO_READ,
			    IOobject::NO_WRITE
			)
		    ),
		T_(T),
		p_(p),
		Cpv_(Cpv),
		Y_(Y),
		W_(W),
		species_names_(species_names)
	{
		Info << "Created object OpenSMOKEppReactingMixtureForRadiation" << endl;
	}

	unsigned int OpenSMOKEppReactingMixtureForRadiation::species(const std::string name) const
	{
		for (unsigned int i=0;i<species_names_.size();i++)
		{
			if (name == species_names_[i])
				return i;
		}
	
		Info << "Species " << name << " not available!" << endl;
		abort(FatalError);

		return 0; 
	}
}

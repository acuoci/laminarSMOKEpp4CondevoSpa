/*-----------------------------------------------------------------------*\
|                                                                         |
|   ╭╮╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╭━━━┳━╮╭━┳━━━┳╮╭━┳━━━╮                               |
|   ┃┃╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱┃╭━╮┃┃╰╯┃┃╭━╮┃┃┃╭┫╭━━╯                               |
|   ┃┃╭━━┳╮╭┳┳━╮╭━━┳━┫╰━━┫╭╮╭╮┃┃╱┃┃╰╯╯┃╰━━┳╮╱╭╮                           |
|   ┃┃┃╭╮┃╰╯┣┫╭╮┫╭╮┃╭┻━━╮┃┃┃┃┃┃┃╱┃┃╭╮┃┃╭━┳╯╰┳╯╰╮                          |
|   ┃╰┫╭╮┃┃┃┃┃┃┃┃╭╮┃┃┃╰━╯┃┃┃┃┃┃╰━╯┃┃┃╰┫╰━┻╮╭┻╮╭╯                          |
|   ╰━┻╯╰┻┻┻┻┻╯╰┻╯╰┻╯╰━━━┻╯╰╯╰┻━━━┻╯╰━┻━━━┻╯╱╰╯                           |
|                                                                         |
|   Authors: Alberto Cuoci                                                |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of laminarSMOKE++ solver.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2022 Alberto Cuoci                                       |
|   laminarSMOKE++ is free software: you can redistribute it and/or       |
|   modify it under the terms of the GNU General Public License           |
|   as published by the Free Software Foundation, either version 3 of     |
|   the License, or (at your option) any later version.                   |
|                                                                         |
|   laminarSMOKE++ is distributed in the hope that it will be useful,     |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with laminarSMOKE++.                                            |
|   If not, see <http://www.gnu.org/licenses/>.                           |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include <queue>

// Default name of main dictionary
const Foam::word Foam::TiO2::dictName("TiO2");

// Constructor
Foam::TiO2::TiO2(const Foam::fvMesh& mesh, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap)
:
	IOdictionary
	(
		IOobject
		(
			dictName,
			mesh.time().constant(),
			mesh,
			IOobject::MUST_READ_IF_MODIFIED,
			IOobject::NO_WRITE
		)
	),

	thermodynamicsMap_(thermodynamicsMap), 

	kineticsMap_(kineticsMap), 

	Ys_
	(
	    IOobject
	    (
		"TiO2:Ys",
		mesh.time().timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	),

	Norm_N_
	(
	    IOobject
	    (
		"TiO2:Norm_N",
		mesh.time().timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	),

	S_
	(
	    IOobject
	    (
		"TiO2:S",
		mesh.time().timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	),

	N_
	(
	    IOobject
	    (
		"TiO2:N",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("TiO2:N", dimensionSet(0., -3., 0., 0., 0.), 0.)
	),

	gamma_
	(
	    IOobject
	    (
		"TiO2:gamma",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("TiO2:gamma", dimensionSet(1., -1., -1., 0., 0.), 0.)
	),

	fv_
	(
	    IOobject
	    (
		"TiO2:fv",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("TiO2:fv", dimensionSet(0., 0., 0., 0., 0.), 0.)
	),

	dp_
	(
	    IOobject
	    (
		"TiO2:dp",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("TiO2:dp", dimensionSet(0., 1., 0., 0., 0.), 0.)
	),

	dc_
	(
	    IOobject
	    (
		"TiO2:dc",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("TiO2:dc", dimensionSet(0., 1., 0., 0., 0.), 0.)
	),

	da_
	(
	    IOobject
	    (
		"TiO2:da",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("TiO2:da", dimensionSet(0., 1., 0., 0., 0.), 0.)
	),

	ss_
	(
	    IOobject
	    (
		"TiO2:ss",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("TiO2:ss", dimensionSet(0., 2., 0., 0., 0.), 0.)
	),

	np_
	(
	    IOobject
	    (
		"TiO2:np",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("TiO2:np", dimensionSet(0., 0., 0., 0., 0.), 0.)
	),

	taus_
	(
	    IOobject
	    (
		"TiO2:taus",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("TiO2:taus", dimensionSet(0., 0., 1., 0., 0.), 0.)
	)
{
	// Kinetic mechanism
	NS_ = thermodynamicsMap_.NumberOfSpecies();
	NR_ = kineticsMap_.NumberOfReactions();

	// User-defined options
	word name_TiO2(this->lookup("TiO2"));
	word name_TiOH4(this->lookup("TiOH4"));
	word name_H2O(this->lookup("H2O"));
	index_tio2_ = thermodynamicsMap_.IndexOfSpecies(name_TiO2)-1;
	index_tioh4_ = thermodynamicsMap_.IndexOfSpecies(name_TiOH4)-1;
	index_h2o_ = thermodynamicsMap_.IndexOfSpecies(name_H2O)-1;	

	// Schmidt number
	schmidt_number_ = readScalar( this->lookup("schmidtNumber"));

	// Diffusion type
	{	
		Foam::string option = this->lookup("diffusionType");
		Info << option << endl;
		if (option == "fixedSchmidtNumber")	diffusion_type_ = DiffusionType::FIXED_SCHMIDT_NUMBER;
		else if (option == "collisionDiameter")	diffusion_type_ = DiffusionType::COLLISION_DIAMETER;
		else FatalError << "Wrong diffusionType option. Allowed options: fixedSchmidtNumber | collisionDiameter" << ::Foam::exit(FatalError); 
	}

	// Thermopheretic effect
	is_thermophoretic_effect_ = Switch(this->lookup(word("thermophoreticEffect")));

	// Activation of processes
	is_nucleation_ = Switch(this->lookup(word("nucleation")));
	is_coagulation_ = Switch(this->lookup(word("coagulation")));
	is_sintering_ = Switch(this->lookup(word("sintering")));
	
	// Feedback on the gaseous phase
	is_gaseous_feedback_ = Switch(this->lookup(word("gaseousFeedback")));

	// On-the-fly analysis
	ontheflyAnalysis_ = false;
	frequency_ = 10;
	frequencyCounter_ = 10;	
	ontheflyAnalysis_ = Switch(this->lookup(word("onTheFlyAnalysis")));
	frequency_  = readLabel(this->lookup("frequency"));
	frequencyCounter_ = frequency_;


	// Setup
	model_tio2_ = new ModelTiO2();

	// Additional user-defined parameters
	{
		// Activation of subprocesses
		model_tio2_->SetNucleation(is_nucleation_);
		model_tio2_->SetCoagulation(is_coagulation_);
		model_tio2_->SetSintering(is_sintering_);

		// Number of monomeric units
		{
			const label np_monomer = readLabel(this->lookup("npMonomer"));
			model_tio2_->SetNumberOfParticles(np_monomer);
		}

		// Number nucleation kinetic parameters
		{
			const scalar epsilon = readScalar( this->lookup("epsilon") );
			model_tio2_->SetNucleationCollisionEnhancementFactor(epsilon);
			const scalar na = readScalar( this->lookup("na") );
			model_tio2_->SetNucleationTemperatureExponent(na);
			const scalar Ta = readScalar( this->lookup("Ta") );
			model_tio2_->SetNucleationActivationTemperature(Ta);
		}

		// Number sintering kinetic parameters
		{
			const scalar As = readScalar( this->lookup("As") );
			model_tio2_->SetSinteringFrequencyFactor(As);
			const scalar ns = readScalar( this->lookup("ns") );
			model_tio2_->SetSinteringTemperatureExponent(ns);
			const scalar Ts = readScalar( this->lookup("Ts") );
			model_tio2_->SetSinteringActivationTemperature(Ts);
		}

		// Numerical parameters
		{
			const scalar threshold_N = readScalar( this->lookup("thresholdN") );
			model_tio2_->SetNThreshold(threshold_N);
			const scalar threshold_fv = readScalar( this->lookup("thresholdFv") );
			model_tio2_->SetFvThreshold(threshold_fv);
		}
	}


	// Soot integrals
	if (ontheflyAnalysis_ == true)
	{
	}

	
	// Summary on the screen
	{
		std::cout << "------------------------------------------------------------------------------------------" << std::endl; 
		std::cout << "                               TiO2 Solver Summary                                        " << std::endl;
		std::cout << "------------------------------------------------------------------------------------------" << std::endl; 

		Info << "TiO2 species:    " << name_TiO2 << " (" << index_tio2_ << ")" << endl;
		Info << "Ti(OH)4 species: " << name_TiOH4 << " (" << index_tioh4_ << ")" << endl;
		Info << "H2O species:     " << name_H2O << " (" << index_h2o_ << ")" << endl;
		Info << "Thermophoresis:  " << is_thermophoretic_effect_ << endl;
		Info << "Diffusion type:  " << diffusion_type_ << endl;
		Info << "Feeback on gas:  " << is_gaseous_feedback_ << endl;
	
		model_tio2_->SummaryOnScreen();

		std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	}
}

Foam::TiO2::~TiO2()
{}

const Foam::volScalarField& Foam::TiO2::Ys() const
{
	return Ys_;
}

Foam::volScalarField& Foam::TiO2::Ys()
{
	return Ys_;
}

const Foam::volScalarField& Foam::TiO2::S() const
{
	return S_;
}

Foam::volScalarField& Foam::TiO2::S()
{
	return S_;
}

const Foam::volScalarField& Foam::TiO2::Norm_N() const
{
	return Norm_N_;
}

Foam::volScalarField& Foam::TiO2::Norm_N()
{
	return Norm_N_;
}

void Foam::TiO2::Setup(const Foam::fvMesh& mesh)
{
	
}

void Foam::TiO2::SolveEquations(const Foam::fvMesh& mesh, const volScalarField& p, const volScalarField& T, const surfaceScalarField& phi, const volScalarField& rho, const volScalarField& mu)
{
	// Equations
	{
		fvScalarMatrix YsEqn
		(
			fvm::ddt(rho,Ys_)
		      + fvm::div(phi, Ys_) 
		      - fvm::laplacian(gamma_, Ys_) 
		);

		// Add thermophoretic effect
		if (is_thermophoretic_effect_ == true)
			YsEqn -= fvc::laplacian(0.55*mu/T*Ys_, T, "laplacian(teta,TiO2)");

		// Solve
		YsEqn.relax();
		YsEqn.solve("Ys");
		
		// Correct
		Ys_.max(0.);
		Ys_.min(1.);
	}

	{
		fvScalarMatrix Norm_NEqn
		(
			fvm::ddt(rho,Norm_N_)
		      + fvm::div(phi, Norm_N_) 
		      - fvm::laplacian(gamma_, Norm_N_) 
		);

		// Add thermophoretic effect
		if (is_thermophoretic_effect_ == true)
			Norm_NEqn -= fvc::laplacian(0.55*mu/T*Norm_N_, T, "laplacian(teta,TiO2)");

		// Solve
		Norm_NEqn.relax();
		Norm_NEqn.solve("Norm_N");

		// Correct
		Norm_N_.max(0.);
		Norm_N_.min(1e0);
	}

	{
		fvScalarMatrix SEqn
		(
			fvm::ddt(rho,S_)
		      + fvm::div(phi, S_) 
		      - fvm::laplacian(gamma_, S_) 
		);

		// Add thermophoretic effect
		if (is_thermophoretic_effect_ == true)
			SEqn -= fvc::laplacian(0.55*mu/T*S_, T, "laplacian(teta,TiO2)");

		// Solve
		SEqn.relax();
		SEqn.solve("S");
		
		// Correct
		S_.max(0.);
		S_.min(1.e9);
	}
}

void Foam::TiO2::UpdateProperties(const volScalarField& p, const volScalarField& T, const volScalarField& rho, const volScalarField& mu)
{
	if (diffusion_type_ == FIXED_SCHMIDT_NUMBER)
	{
		gamma_ = mu/schmidt_number_;
	}
	else if (diffusion_type_ == COLLISION_DIAMETER)
	{
		const scalarField& TCells = T.internalField();
		const scalarField& rhoCells = rho.internalField();
		const scalarField& muCells = mu.internalField();

		const scalarField& YsCells = Ys_.internalField();
		const scalarField& Norm_NCells = Norm_N_.internalField();
		const scalarField& SCells = S_.internalField();
		
		scalarField& gammaCells  = gamma_.ref();

		forAll(TCells, celli)
		{
			model_tio2_->DiffusionCoefficient(	TCells[celli], rhoCells[celli], muCells[celli], schmidt_number_,
								YsCells[celli], Norm_NCells[celli]*6.02214086e26, SCells[celli],
								gammaCells[celli] );	
		}

		forAll(T.boundaryField(), patchi)
		{
		
			const fvPatchScalarField& pT  = T.boundaryField()[patchi];
			const fvPatchScalarField& prho  = rho.boundaryField()[patchi];
			const fvPatchScalarField& pmu  = mu.boundaryField()[patchi];

			const fvPatchScalarField& pYs  = Ys_.boundaryField()[patchi];
			const fvPatchScalarField& pNorm_N = Norm_N_.boundaryField()[patchi];
			const fvPatchScalarField& pS = S_.boundaryField()[patchi];
				
			fvPatchScalarField& pgamma  = gamma_.boundaryFieldRef()[patchi];

			forAll(pT, facei)
			{
				model_tio2_->DiffusionCoefficient(	pT[facei], prho[facei], pmu[facei], schmidt_number_,
									pYs[facei], pNorm_N[facei]*6.02214086e26, pS[facei],
									pgamma[facei] );
			}
		}
	}
}

void Foam::TiO2::OnTheFlyAnalysis(const Foam::fvMesh& mesh, const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	if (ontheflyAnalysis_ == true && frequencyCounter_ == frequency_)
	{
		frequencyCounter_ = 0;
	}

	frequencyCounter_++;
}

void Foam::TiO2::FieldsReconstruction(const Foam::fvMesh& mesh, const volScalarField& rho, const volScalarField& T)
{
	if (mesh.time().outputTime() == true)
	{
		const scalarField& TCells = T.internalField();
		const scalarField& rhoCells = rho.internalField();

		const scalarField& YsCells = Ys_.internalField();
		const scalarField& Norm_NCells = Norm_N_.internalField();
		const scalarField& SCells = S_.internalField();

		scalarField& NCells  = N_.ref();

		scalarField& fvCells = fv_.ref(); 
		scalarField& dpCells = dp_.ref();
		scalarField& dcCells = dc_.ref(); 
		scalarField& daCells = da_.ref();
		scalarField& ssCells = ss_.ref(); 
		scalarField& npCells = np_.ref();
		scalarField& tausCells = taus_.ref();

		forAll(TCells, celli)
		{
			// Calculating soot particle number density [#/m3]
			NCells[celli] = Norm_NCells[celli]*6.02214086e26;

			// Calculate properties
			model_tio2_->Properties(	TCells[celli], rhoCells[celli], 
							YsCells[celli], NCells[celli], SCells[celli],
							fvCells[celli], dpCells[celli], dcCells[celli], daCells[celli], ssCells[celli], npCells[celli], tausCells[celli] );
		}

		forAll(T.boundaryField(), patchi)
		{
			const fvPatchScalarField& pT  = T.boundaryField()[patchi];
			const fvPatchScalarField& prho  = rho.boundaryField()[patchi];

			const fvPatchScalarField& pYs  = Ys_.boundaryField()[patchi];
			const fvPatchScalarField& pNorm_N = Norm_N_.boundaryField()[patchi];
			const fvPatchScalarField& pS = S_.boundaryField()[patchi];
			
			fvPatchScalarField& pN  = N_.boundaryFieldRef()[patchi];

			fvPatchScalarField& pfv = fv_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pdp = dp_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pdc = dc_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pda = da_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pss = ss_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pnp = np_.boundaryFieldRef()[patchi];
			fvPatchScalarField& ptaus = taus_.boundaryFieldRef()[patchi];

			forAll(pT, facei)
			{

				// Calculating soot particle number density [#/m3]
				pN[facei] = pNorm_N[facei]*6.02214086e26;

				// Calculate properties
				model_tio2_->Properties(	pT[facei], prho[facei], 
								pYs[facei], pN[facei], pS[facei],
								pfv[facei], pdp[facei], pdc[facei], pda[facei], pss[facei], pnp[facei], ptaus[facei] );
			}
		}
	}
}

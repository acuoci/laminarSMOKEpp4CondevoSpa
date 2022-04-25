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
|   Copyright(C) 2020 Alberto Cuoci                                       |
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


// Default name of main dictionary
const Foam::word Foam::HMOMModel::dictName("HMOMModel");

// Constructor
Foam::HMOMModel::HMOMModel(const Foam::fvMesh& mesh, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap, OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap)
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

	transportMap_(transportMap),


	Norm_N0_
	(
	    IOobject
	    (
		"hmom:Norm_N0",
		mesh.time().timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	),
	
	Norm_M00_
	(
	    IOobject
	    (
		"hmom:Norm_M00",
		mesh.time().timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	),

	Norm_M10_
	(
	    IOobject
	    (
		"hmom:Norm_M10",
		mesh.time().timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	),

	Norm_M01_
	(
	    IOobject
	    (
		"hmom:Norm_M01",
		mesh.time().timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	),

	fv_
	(
	    IOobject
	    (
		"hmom:fv",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("hmom:fv", dimensionSet(0., 0., 0., 0., 0.), 0.)
	),

	N_
	(
	    IOobject
	    (
		"hmom:N",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("hmom:N", dimensionSet(0., -3., 0., 0., 0.), 0.)
	),

	dp_
	(
	    IOobject
	    (
		"hmom:dp",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("hmom:dp", dimensionSet(0., 1., 0., 0., 0.), 0.)
	),

	dc_
	(
	    IOobject
	    (
		"hmom:dc",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("hmom:dc", dimensionSet(0., 1., 0., 0., 0.), 0.)
	),

	np_
	(
	    IOobject
	    (
		"hmom:np",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("hmom:np", dimensionSet(0., 0., 0., 0., 0.), 0.)
	),

	gamma_
	(
	    IOobject
	    (
		"hmom:gamma",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("hmom:gamma", dimensionSet(1., -1., -1., 0., 0.), 0.)
	),

	soot_fv_
	(
	    IOobject
	    (
		"soot:fv",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("soot:fv", dimensionSet(0., 0., 0., 0., 0.), 0.)
	)

{
	// Default values
	Switch hmom_analysis = false;
	Switch hmom_analysis_onthefly = true;
	Switch hmom_test = false;
	hmom_pah_gas_consumption_ = true;

	int hmom_pah_nc = 16;
	int hmom_number_moments = 4;
	int hmom_collision_diameter_model = 2;
	int hmom_fractal_diameter_model = 1;

	int hmom_nucleation = 1;
	int hmom_surface_growth = 1;
	int hmom_oxidation = 1;
	int hmom_condensation = 1;
	int hmom_coagulation = 1;
	int hmom_coagulation_continous = 1;
	double hmom_schmidt_number = 50.;


	double hmom_sticking_coefficient = 0.002;
	double hmom_soot_density = 1800.;
	int hmom_thermophoretic_model = 1;

	hmom_minimum_normalized_value_ = 1.e-16;
	hmom_minimum_fv_for_integration_ = 1.e-16;

	// Surface density correction
	double hmom_surface_density = 1.7e19;
	Switch hmom_surface_density_correction_coefficient = false;	
	double hmom_surface_density_coeff_a1 = 12.65;		// [-]
	double hmom_surface_density_coeff_a2 = 0.00563; 	// [1/K]
	double hmom_surface_density_coeff_b1 = 1.38; 		// [-]
	double hmom_surface_density_coeff_b2 = 0.00069;		// [1/K]

	// Kinetics (default values)
	double hmom_A1f = 6.72e1;
	double hmom_n1f = 3.33;
	double hmom_E1f = 6.09;
	double hmom_A1b = 6.44e-1;
	double hmom_n1b = 3.79;
	double hmom_E1b = 27.96;
	double hmom_A2f = 1.00e8;
	double hmom_n2f = 1.80;
	double hmom_E2f = 68.42;
	double hmom_A2b = 8.68e4;
	double hmom_n2b = 2.36;
	double hmom_E2b = 25.46;
	double hmom_A3f = 1.13e16;
	double hmom_n3f = -0.06;
	double hmom_E3f = 476.05;
	double hmom_A3b = 4.17e13;
	double hmom_n3b = 0.15;
	double hmom_E3b = 0.00;
	double hmom_A4 = 2.52e9;
	double hmom_n4 = 1.10;
	double hmom_E4 = 17.13;
	double hmom_A5 = 2.20e12;
	double hmom_n5 = 0.00;
	double hmom_E5 = 31.38;
	double hmom_eff6 = 0.13;

	// Post-processing
	ontheflyAnalysis_ = false;
	frequency_ = 10;
	frequencyCounter_ = 10;	


	// Recover data 
	hmom_index_H_ = thermodynamicsMap_.IndexOfSpecies("H")-1;
	hmom_index_OH_ = thermodynamicsMap_.IndexOfSpecies("OH")-1;
	hmom_index_H2O_ = thermodynamicsMap_.IndexOfSpecies("H2O")-1;
	hmom_index_H2_ = thermodynamicsMap_.IndexOfSpecies("H2")-1;
	hmom_index_O2_ = thermodynamicsMap_.IndexOfSpecies("O2")-1;
	hmom_index_C2H2_ = thermodynamicsMap_.IndexOfSpecies("C2H2")-1;

	List<word>  listSpecies(this->lookup("pah_species"));
	hmom_pah_species_indices_.resize(listSpecies.size());
	std::vector<std::string> hmom_pah_species_names(listSpecies.size());
	for (int i=0;i<listSpecies.size();i++)
	{
		hmom_pah_species_indices_[i] = thermodynamicsMap_.IndexOfSpecies(listSpecies[i])-1;
		hmom_pah_species_names[i] = listSpecies[i];
	}

	hmom_pah_nc = readLabel(this->lookup("pah_nc"));
	hmom_collision_diameter_model = readLabel(this->lookup("collision_diameter_model"));
	hmom_fractal_diameter_model = readLabel(this->lookup("fractal_diameter_model"));

	hmom_nucleation = readLabel(this->lookup("nucleation"));
	hmom_surface_growth = readLabel(this->lookup("surface_growth"));
	hmom_oxidation = readLabel(this->lookup("oxidation"));
	hmom_condensation = readLabel(this->lookup("condensation"));
	hmom_coagulation = readLabel(this->lookup("coagulation"));
	hmom_coagulation_continous = readLabel(this->lookup("coagulation_continous"));

	hmom_pah_gas_consumption_ = Switch(this->lookup(word("pah_gas_consumption")));

	hmom_schmidt_number = this->lookupOrDefault<scalar>(word("schmidt_number"),  scalar(50.));

	hmom_sticking_coefficient = this->lookupOrDefault<scalar>(word("sticking_coefficient"),  scalar(0.002));
	hmom_soot_density = this->lookupOrDefault<scalar>(word("soot_density"),  scalar(1800.));
	hmom_thermophoretic_model = this->lookupOrDefault<label>(word("thermophoretic_model"),  label(1));

	hmom_minimum_normalized_value_   = this->lookupOrDefault<scalar>(word("mininum_normalized_value"),  scalar(1.e-16));
	hmom_minimum_fv_for_integration_ = this->lookupOrDefault<scalar>(word("mininum_fv_for_integration"),  scalar(1.e-16));

	hmom_A1f   = this->lookupOrDefault<scalar>(word("A1f"),  scalar(hmom_A1f));
	hmom_A1b   = this->lookupOrDefault<scalar>(word("A1b"),  scalar(hmom_A1b));
	hmom_A2f   = this->lookupOrDefault<scalar>(word("A2f"),  scalar(hmom_A2f));
	hmom_A2b   = this->lookupOrDefault<scalar>(word("A2b"),  scalar(hmom_A2b));
	hmom_A3f   = this->lookupOrDefault<scalar>(word("A3f"),  scalar(hmom_A3f));
	hmom_A3b   = this->lookupOrDefault<scalar>(word("A3b"),  scalar(hmom_A3b));	
	hmom_A4    = this->lookupOrDefault<scalar>(word("A4"),   scalar(hmom_A4));
	hmom_A5    = this->lookupOrDefault<scalar>(word("A5"),   scalar(hmom_A5));	

	hmom_E1f   = this->lookupOrDefault<scalar>(word("E1f"),  scalar(hmom_E1f));
	hmom_E1b   = this->lookupOrDefault<scalar>(word("E1b"),  scalar(hmom_E1b));
	hmom_E2f   = this->lookupOrDefault<scalar>(word("E2f"),  scalar(hmom_E2f));
	hmom_E2b   = this->lookupOrDefault<scalar>(word("E2b"),  scalar(hmom_E2b));
	hmom_E3f   = this->lookupOrDefault<scalar>(word("E3f"),  scalar(hmom_E3f));
	hmom_E3b   = this->lookupOrDefault<scalar>(word("E3b"),  scalar(hmom_E3b));	
	hmom_E4    = this->lookupOrDefault<scalar>(word("E4"),   scalar(hmom_E4));
	hmom_E5    = this->lookupOrDefault<scalar>(word("E5"),   scalar(hmom_E5));	

	hmom_n1f   = this->lookupOrDefault<scalar>(word("n1f"),  scalar(hmom_n1f));
	hmom_n1b   = this->lookupOrDefault<scalar>(word("n1b"),  scalar(hmom_n1b));
	hmom_n2f   = this->lookupOrDefault<scalar>(word("n2f"),  scalar(hmom_n2f));
	hmom_n2b   = this->lookupOrDefault<scalar>(word("n2b"),  scalar(hmom_n2b));
	hmom_n3f   = this->lookupOrDefault<scalar>(word("n3f"),  scalar(hmom_n3f));
	hmom_n3b   = this->lookupOrDefault<scalar>(word("n3b"),  scalar(hmom_n3b));	
	hmom_n4    = this->lookupOrDefault<scalar>(word("n4"),   scalar(hmom_n4));
	hmom_n5    = this->lookupOrDefault<scalar>(word("n5"),   scalar(hmom_n5));	

	hmom_eff6  = this->lookupOrDefault<scalar>(word("Eff6"),   scalar(hmom_eff6));

	hmom_surface_density = this->lookupOrDefault<scalar>(word("surface_density"),  scalar(hmom_surface_density));
	hmom_surface_density_correction_coefficient = Switch(this->lookup(word("surface_density_correction")));
	hmom_surface_density_coeff_a1 = this->lookupOrDefault<scalar>(word("surface_density_a1"),  scalar(hmom_surface_density_coeff_a1));
	hmom_surface_density_coeff_a2 = this->lookupOrDefault<scalar>(word("surface_density_a2"),  scalar(hmom_surface_density_coeff_a2));
	hmom_surface_density_coeff_b1 = this->lookupOrDefault<scalar>(word("surface_density_b1"),  scalar(hmom_surface_density_coeff_b1));
	hmom_surface_density_coeff_b2 = this->lookupOrDefault<scalar>(word("surface_density_b2"),  scalar(hmom_surface_density_coeff_b2));

	Info << "HMOM: surface density (#/m2): " << hmom_surface_density << endl;
	Info << "HMOM: sticking coefficient:   " << hmom_sticking_coefficient << endl;
	Info << "HMOM: soot density (kg/m3):   " << hmom_soot_density  << endl;
	Info << "HMOM: thermophoresis:         " << hmom_thermophoretic_model  << endl;

	Info << "HMOM: minimum normalized moments: " << hmom_minimum_normalized_value_ << endl;
	Info << "HMOM: minimum fv for integration: " << hmom_minimum_fv_for_integration_ << endl;

	// On-the-fly analysis
	ontheflyAnalysis_    = Switch(this->lookup(word("onTheFlyAnalysis")));
	frequency_  = readLabel(this->lookup("frequency"));
	frequencyCounter_ = frequency_;

	// Soot integrals
	if (ontheflyAnalysis_ == true)
	{
		if ( !boost::filesystem::exists("soot.integrals") )
		{
			fSootIntegrals_.reset(new std::ofstream("soot.integrals", std::ios::out));
			fSootIntegrals_().setf(std::ios::scientific);

			fSootIntegrals_() << std::setw(20) << "time";
			fSootIntegrals_() << std::setw(20) << "soot_mass[kg]";
			fSootIntegrals_() << std::setw(20) << "soot_volume[m3]";
			fSootIntegrals_() << std::setw(20) << "soot_particles[#]";		
			fSootIntegrals_() << std::setw(20) << "pah_mass[kg]";
			fSootIntegrals_() << std::endl;
		}
		else
		{
			fSootIntegrals_.reset(new std::ofstream("soot.integrals", std::ios::app));
			fSootIntegrals_().setf(std::ios::scientific);
		}
	}

	// Setup
	hmomAnalyzer_ = new OpenSMOKE::HMOM();

	hmomAnalyzer_->SetSurfaceDensity(hmom_surface_density);
	hmomAnalyzer_->SetSurfaceDensityCorrectionCoefficient(hmom_surface_density_correction_coefficient);
	hmomAnalyzer_->SetSurfaceDensityCorrectionCoefficientA1(hmom_surface_density_coeff_a1);
	hmomAnalyzer_->SetSurfaceDensityCorrectionCoefficientA2(hmom_surface_density_coeff_a2);
	hmomAnalyzer_->SetSurfaceDensityCorrectionCoefficientB1(hmom_surface_density_coeff_b1);
	hmomAnalyzer_->SetSurfaceDensityCorrectionCoefficientB2(hmom_surface_density_coeff_b2);

	hmomAnalyzer_->SetA1f(hmom_A1f);
	hmomAnalyzer_->SetA1b(hmom_A1b);
	hmomAnalyzer_->SetA2f(hmom_A2f);
	hmomAnalyzer_->SetA2b(hmom_A2b);
	hmomAnalyzer_->SetA3f(hmom_A3f);
	hmomAnalyzer_->SetA3b(hmom_A3b);
	hmomAnalyzer_->SetA4(hmom_A4);
	hmomAnalyzer_->SetA5(hmom_A5);

	hmomAnalyzer_->SetE1f(hmom_E1f);
	hmomAnalyzer_->SetE1b(hmom_E1b);
	hmomAnalyzer_->SetE2f(hmom_E2f);
	hmomAnalyzer_->SetE2b(hmom_E2b);
	hmomAnalyzer_->SetE3f(hmom_E3f);
	hmomAnalyzer_->SetE3b(hmom_E3b);
	hmomAnalyzer_->SetE4(hmom_E4);
	hmomAnalyzer_->SetE5(hmom_E5);

	hmomAnalyzer_->Setn1f(hmom_n1f);
	hmomAnalyzer_->Setn1b(hmom_n1b);
	hmomAnalyzer_->Setn2f(hmom_n2f);
	hmomAnalyzer_->Setn2b(hmom_n2b);
	hmomAnalyzer_->Setn3f(hmom_n3f);
	hmomAnalyzer_->Setn3b(hmom_n3b);
	hmomAnalyzer_->Setn4(hmom_n4);
	hmomAnalyzer_->Setn5(hmom_n5);

	hmomAnalyzer_->SetEfficiency6(hmom_eff6);

	hmomAnalyzer_->SetNucleation(hmom_nucleation);
	hmomAnalyzer_->SetSurfaceGrowth(hmom_surface_growth);
	hmomAnalyzer_->SetOxidation(hmom_oxidation);
	hmomAnalyzer_->SetCondensation(hmom_condensation);
	hmomAnalyzer_->SetCoagulation(hmom_coagulation);
	hmomAnalyzer_->SetCoagulationContinous(hmom_coagulation_continous);
	hmomAnalyzer_->SetThermophoreticModel(hmom_thermophoretic_model);
	hmomAnalyzer_->SetPAH(hmom_pah_species_names);

	hmomAnalyzer_->SetNumberCarbonPAH(hmom_pah_nc);
	hmomAnalyzer_->SetFractalDiameterModel(hmom_fractal_diameter_model);
	hmomAnalyzer_->SetCollisionDiameterModel(hmom_collision_diameter_model);
	hmomAnalyzer_->SetSchmidtNumber(hmom_schmidt_number);
	hmomAnalyzer_->SetStickingCoefficient(hmom_sticking_coefficient);
	hmomAnalyzer_->SetSootDensity(hmom_soot_density);

	hmomAnalyzer_->Summary();
}

Foam::HMOMModel::~HMOMModel()
{}

const Foam::volScalarField& Foam::HMOMModel::Norm_N0() const
{
	return Norm_N0_;
}

Foam::volScalarField& Foam::HMOMModel::Norm_N0()
{
	return Norm_N0_;
}

const Foam::volScalarField& Foam::HMOMModel::Norm_M00() const
{
	return Norm_M00_;
}

Foam::volScalarField& Foam::HMOMModel::Norm_M00()
{
	return Norm_M00_;
}

const Foam::volScalarField& Foam::HMOMModel::Norm_M10() const
{
	return Norm_M10_;
}

Foam::volScalarField& Foam::HMOMModel::Norm_M10()
{
	return Norm_M10_;
}

const Foam::volScalarField& Foam::HMOMModel::Norm_M01() const
{
	return Norm_M01_;
}

Foam::volScalarField& Foam::HMOMModel::Norm_M01()
{
	return Norm_M01_;
}

const Foam::volScalarField& Foam::HMOMModel::fv() const
{
	return fv_;
}

Foam::volScalarField& Foam::HMOMModel::fv()
{
	return fv_;
}

void Foam::HMOMModel::Setup(const Foam::fvMesh& mesh)
{
	soot_fv_ = fv_;
}

void Foam::HMOMModel::UpdateProperties(const volScalarField& p, const volScalarField& T)
{

	const scalarField& TCells = T.internalField();

	scalarField& Norm_N0Cells  = Norm_N0_.ref();
	scalarField& Norm_M00Cells = Norm_M00_.ref();
	scalarField& Norm_M01Cells = Norm_M01_.ref();
	scalarField& Norm_M10Cells = Norm_M10_.ref();

	scalarField& fvCells = fv_; 
	scalarField& NCells = N_; 

	forAll(TCells, celli)
	{
		hmomAnalyzer_->SetNormalizedMoments(Norm_M00Cells[celli], Norm_M10Cells[celli], Norm_M01Cells[celli], Norm_N0Cells[celli]);
		fvCells[celli] = hmomAnalyzer_->SootVolumeFraction();
		NCells[celli] = hmomAnalyzer_->SootParticleNumberDensity();
	}

	forAll(T.boundaryField(), patchi)
	{
		const fvPatchScalarField& pT  = T.boundaryField()[patchi];
	
		fvPatchScalarField& pNorm_N0  = Norm_N0_.boundaryFieldRef()[patchi];
		fvPatchScalarField& pNorm_M00 = Norm_M00_.boundaryFieldRef()[patchi];
		fvPatchScalarField& pNorm_M01 = Norm_M01_.boundaryFieldRef()[patchi];
		fvPatchScalarField& pNorm_M10 = Norm_M10_.boundaryFieldRef()[patchi];

		fvPatchScalarField& pfv = fv_.boundaryFieldRef()[patchi];
		fvPatchScalarField& pN = N_.boundaryFieldRef()[patchi];

		forAll(pT, facei)
		{
			hmomAnalyzer_->SetNormalizedMoments(pNorm_M00[facei], pNorm_M10[facei], pNorm_M01[facei], pNorm_N0[facei]);
			pfv[facei] = hmomAnalyzer_->SootVolumeFraction();
			pN[facei] = hmomAnalyzer_->SootParticleNumberDensity();
		}
	}

	// This is needed by the radiation model
	soot_fv_ = fv_;
}

void Foam::HMOMModel::OnTheFlyAnalysis(const Foam::fvMesh& mesh, const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	if (ontheflyAnalysis_ == true && frequencyCounter_ == frequency_)
	{
		const double soot_mass = hmomAnalyzer_->SootDensity()*fvc::domainIntegrate(fv_).value();
		const double soot_volume = fvc::domainIntegrate(fv_).value();
		const double soot_particles = fvc::domainIntegrate(N_).value();
		
		double pah_mass = 0.;
		for (unsigned int j = 0; j < hmomAnalyzer_->pah_species().size(); j++)
			pah_mass += fvc::domainIntegrate(rho*Y[hmom_pah_species_indices_[j]]).value();
		
		// Only if the processor is the master, write to the file
		if (Pstream::master())
		{
	    		fSootIntegrals_() << std::setw(20) << mesh.time().value();
			fSootIntegrals_() << std::setw(20) << soot_mass;
			fSootIntegrals_() << std::setw(20) << soot_volume;
			fSootIntegrals_() << std::setw(20) << soot_particles;
			fSootIntegrals_() << std::setw(20) << pah_mass;
			fSootIntegrals_() << std::endl;
		}
		
		Info << "Soot mass [kg]:     " << soot_mass << endl;
		Info << "Soot volume [m3]:   " << soot_volume << endl;
		Info << "Soot particles [#]: " << soot_particles << endl;		
		Info << "PAH mass [kg]:      " << pah_mass << endl;

		frequencyCounter_ = 0;
	}

	frequencyCounter_++;
}

void Foam::HMOMModel::FieldsReconstruction(const Foam::fvMesh& mesh, const volScalarField& p, const volScalarField& T)
{
	if (mesh.time().outputTime() == true)
	{
		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField();

		scalarField& Norm_N0Cells  = Norm_N0_.ref();
		scalarField& Norm_M00Cells = Norm_M00_.ref();
		scalarField& Norm_M01Cells = Norm_M01_.ref();
		scalarField& Norm_M10Cells = Norm_M10_.ref();

		scalarField& fvCells = fv_.ref(); 
		scalarField& NCells  = N_.ref(); 
		scalarField& dpCells = dp_.ref();
		scalarField& dcCells = dc_.ref(); 
		scalarField& npCells = np_.ref(); 

		forAll(TCells, celli)
		{
			hmomAnalyzer_->SetNormalizedMoments(Norm_M00Cells[celli], Norm_M10Cells[celli], Norm_M01Cells[celli], Norm_N0Cells[celli]);

			// Calculating soot volume fraction [-]
			fvCells[celli] = hmomAnalyzer_->SootVolumeFraction();

			// Calculating soot particle number density [#/m3]
			NCells[celli] = hmomAnalyzer_->SootParticleNumberDensity();

			// Calculating particle diameters [nm]
			dpCells[celli] = hmomAnalyzer_->SootParticleDiameter() * 1.e9;

			// Calculating collision diameters [nm]
			dcCells[celli] = hmomAnalyzer_->SootCollisionParticleDiameter() * 1.e9;

			// Calculating number of primary particles [-]
			npCells[celli] = hmomAnalyzer_->SootNumberOfPrimaryParticles();
		}

		forAll(T.boundaryField(), patchi)
		{
			const fvPatchScalarField& pT  = T.boundaryField()[patchi];
			const fvPatchScalarField& pp  = p.boundaryField()[patchi];

			fvPatchScalarField& pNorm_N0  = Norm_N0_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pNorm_M00 = Norm_M00_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pNorm_M01 = Norm_M01_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pNorm_M10 = Norm_M10_.boundaryFieldRef()[patchi];

			fvPatchScalarField& pfv = fv_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pN  = N_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pdp = dp_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pdc = dc_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pnp = np_.boundaryFieldRef()[patchi];
		
			forAll(pT, facei)
			{
				hmomAnalyzer_->SetNormalizedMoments(pNorm_M00[facei], pNorm_M10[facei], pNorm_M01[facei], pNorm_N0[facei]);

				// Calculating soot volume fraction [-]
				pfv[facei] = hmomAnalyzer_->SootVolumeFraction();

				// Calculating soot particle number density [#/m3]
				pN[facei] = hmomAnalyzer_->SootParticleNumberDensity();

				// Calculating particle diameters [nm]
				pdp[facei] = hmomAnalyzer_->SootParticleDiameter() * 1.e9;

				// Calculating collision diameters [nm]
				pdc[facei] = hmomAnalyzer_->SootCollisionParticleDiameter() * 1.e9;
	
				// Calculating number of primary particles [-]
				pnp[facei] = hmomAnalyzer_->SootNumberOfPrimaryParticles();
			}
		}

		// This is needed by the radiation model
		soot_fv_ = fv_;
	}
}

void Foam::HMOMModel::SolveEquations(const Foam::fvMesh& mesh, const volScalarField& p, const volScalarField& T, const surfaceScalarField& phi, const volScalarField& rho, const volScalarField& mu, const Switch thermophoretic_effect)
{
	gamma_ = mu/hmomAnalyzer_->schmidt_number();

	{
		fvScalarMatrix Norm_N0Eqn
		(
			fvm::ddt(rho,Norm_N0_)
		      + fvm::div(phi, Norm_N0_) 
		      - fvm::laplacian(gamma_, Norm_N0_) 
		);

		// Add thermophoretic effect
		if (thermophoretic_effect == true)
			Norm_N0Eqn -= fvc::laplacian(0.55*mu/T*Norm_N0_, T, "laplacian(teta,hmom)");

		// Solve
		Norm_N0Eqn.relax();
		Norm_N0Eqn.solve("Norm_N0");
		
		// Correct
		Norm_N0_.max(hmom_minimum_normalized_value_);
		Norm_N0_.min(1e0);
	}

	{
		fvScalarMatrix Norm_M00Eqn
		(
			fvm::ddt(rho,Norm_M00_)
		      + fvm::div(phi, Norm_M00_) 
		      - fvm::laplacian(gamma_, Norm_M00_) 
		);

		// Add thermophoretic effect
		if (thermophoretic_effect == true)
			Norm_M00Eqn -= fvc::laplacian(0.55*mu/T*Norm_M00_, T, "laplacian(teta,hmom)");

		// Solve
		Norm_M00Eqn.relax();
		Norm_M00Eqn.solve("Norm_M00");

		// Correct
		Norm_M00_.max(hmom_minimum_normalized_value_);
		Norm_M00_.min(1e0);
	}

	{
		fvScalarMatrix Norm_M10Eqn
		(
			fvm::ddt(rho,Norm_M10_)
		      + fvm::div(phi, Norm_M10_) 
		      - fvm::laplacian(gamma_, Norm_M10_) 
		);

		// Add thermophoretic effect
		if (thermophoretic_effect == true)
			Norm_M10Eqn -= fvc::laplacian(0.55*mu/T*Norm_M10_, T, "laplacian(teta,hmom)");

		// Solve
		Norm_M10Eqn.relax();
		Norm_M10Eqn.solve("Norm_M10");

		// Correct
		Norm_M10_.max(hmom_minimum_normalized_value_);
		Norm_M10_.min(1e0);
	}

	{
		fvScalarMatrix Norm_M01Eqn
		(
			fvm::ddt(rho,Norm_M01_)
		      + fvm::div(phi, Norm_M01_) 
		      - fvm::laplacian(gamma_, Norm_M01_)
		);

		// Add thermophoretic effect
		if (thermophoretic_effect == true)
			Norm_M01Eqn -= fvc::laplacian(0.55*mu/T*Norm_M01_, T, "laplacian(teta,hmom)");

		// Solve
		Norm_M01Eqn.relax();
		Norm_M01Eqn.solve("Norm_M01");

		// Correct
		Norm_M01_.max(hmom_minimum_normalized_value_);
		Norm_M01_.min(1e0);
	}
}

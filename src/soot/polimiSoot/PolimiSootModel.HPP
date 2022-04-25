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
const Foam::word Foam::PolimiSootModel::dictName("PolimiSootModel");

// Constructor
Foam::PolimiSootModel::PolimiSootModel(const Foam::fvMesh& mesh, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap, OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap)
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
	transportMap_(transportMap)
{
	// Reconstruction of soot properties
	ontheflyAnalysis_ = true;
	frequency_ = 10;
	frequencyCounter_ = 10;
	ontheflyAnalysis_    = Switch(this->lookup(word("onTheFlyAnalysis")));
	frequency_  = readLabel(this->lookup("frequency"));
	frequencyCounter_ = frequency_;

	// Thermophoretic effect
	thermophoreticEffect_ = true;
	label thermophoreticEffectMinimumBin = 1;
	thermophoreticEffect_ = Switch(this->lookup(word("thermophoreticEffect")));
	thermophoreticEffectMinimumBin = readLabel(this->lookup("thermophoreticEffectMinimumBin"));

	// Physical diffusion
	physicalDiffusion_ = true;
	label physicalDiffusionBinToStart = 5;
	label physicalDiffusionBinToCut = 10;
	physicalDiffusion_ = Switch(this->lookup(word("physicalDiffusivity")));
	physicalDiffusionBinToStart = readLabel(this->lookup("physicalDiffusivityBinToStart"));
	physicalDiffusionBinToCut = readLabel(this->lookup("physicalDiffusivityBinToCut"));	

	// Minimum section for spherical particles
	label minimumSectionSphericalParticles = readLabel(this->lookup("minimumSectionSphericalParticles"));

	// Minimum section for aggregates
	label minimumSectionAggregates = readLabel(this->lookup("minimumSectionAggregates"));

	// Soot density
	label binIndexZero    = this->lookupOrDefault("binIndexZero", 10);
	label binIndexFinal   = this->lookupOrDefault("binIndexFinal", 20);
	double binDensityZero  = this->lookupOrDefault("binDensityZero", 1500.);
	double binDensityFinal = this->lookupOrDefault("binDensityFinal", 1700.);

	// Fractal diameter
	double fractalDiameter  = this->lookupOrDefault("fractalDiameter", 1.80);

	sootAnalyzer_ = new OpenSMOKE::PolimiSoot_Analyzer(&thermodynamicsMap_);
	sootAnalyzer_->Initialize();
	sootAnalyzer_->SetLabel("BIN");
	sootAnalyzer_->SetFractalDiameter(fractalDiameter);
	sootAnalyzer_->SetMinimumSectionSphericalParticles(minimumSectionSphericalParticles);
	sootAnalyzer_->SetMinimumSectionAggregates(minimumSectionAggregates);
	sootAnalyzer_->SetDensity(binIndexZero, binIndexFinal, binDensityZero, binDensityFinal);
	sootAnalyzer_->SetPhysicalDiffusionBinToStart(physicalDiffusionBinToStart);
	sootAnalyzer_->SetPhysicalDiffusionBinToCut(physicalDiffusionBinToCut);
	sootAnalyzer_->SetThermophoreticEffectMinimumBin(thermophoreticEffectMinimumBin);
	sootAnalyzer_->Setup();

	// Physical diffusivity
	if (physicalDiffusion_ == true)
		for (unsigned int i=0;i<physicalDiffusionCorrectionIndex().size();i++)
		{
			const unsigned int index = physicalDiffusionCorrectionIndex()[i];
			Info << "Physical diffusion coefficient: " << thermodynamicsMap_.NamesOfSpecies()[index] << " : " << physicalDiffusionCorrection()[i] << endl;
		}

	// Thermophoretic effect
	if (thermophoretic_effect_list().size() == 0)
		Info << "WARNING: No BIN species were found in the kinetic mechanism! Are you sure you want to account for the thermophoretic effect?" << endl;
	else
		Info << "Thermophoretic effect: Number of BIN species found equal to " << thermophoretic_effect_list().size() << endl;
}

Foam::PolimiSootModel::~PolimiSootModel()
{}


void Foam::PolimiSootModel::Setup(const Foam::fvMesh& mesh)
{
	// Soot fields 
	{
		sootFields_.resize(6);

		sootFields_.set
		(
			0,
			new volScalarField
			(
				IOobject
				(
					"soot:fv",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("soot:fv", dimensionSet(0, 0, 0, 0, 0), 0.)
			)
		);

		sootFields_.set
		(
			1,
			new volScalarField
			(
				IOobject
				(
					"soot:rho",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("soot:rho", dimensionSet(1, -3, 0, 0, 0), 0.)
			)
		);

		sootFields_.set
		(
			2,
			new volScalarField
			(
				IOobject
				(
					"soot:N",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("soot:N", dimensionSet(0, -3, 0, 0, 0), 0.)
			)
		);

		sootFields_.set
		(
			3,
			new volScalarField
			(
				IOobject
				(
					"soot:X",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("soot:X", dimensionSet(0, 0, 0, 0, 0), 0.)
			)
		);

		sootFields_.set
		(
			4,
			new volScalarField
			(
				IOobject
				(
					"soot:Y",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("soot:Y", dimensionSet(0, 0, 0, 0, 0), 0.)
			)
		);

		sootFields_.set
		(
			5,
			new volScalarField
			(
				IOobject
				(
					"soot:R",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("soot:R", dimensionSet(1, -3, -1, 0, 0), 0.)
			)
		);
	}

	// PAHs fields 
	{
		pahFields_.resize(3);

		pahFields_.set
		(
			0,
			new volScalarField
			(
				IOobject
				(
					"pah:Y12",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("pah:Y12", dimensionSet(0, 0, 0, 0, 0), 0.)
			)
		);

		pahFields_.set
		(
			1,
			new volScalarField
			(
				IOobject
				(
					"pah:Y34",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("pah:Y34", dimensionSet(0, 0, 0, 0, 0), 0.)
			)
		);

		pahFields_.set
		(
			2,
			new volScalarField
			(
				IOobject
				(
					"pah:YLP",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("pah:YLP", dimensionSet(0, 0, 0, 0, 0), 0.)
			)
		);
	}

	// Prepare PolimiSootAnalyzer
	if (ontheflyAnalysis_ == true)
	{
		// Soot integrals
		if ( !boost::filesystem::exists("soot.integrals") )
		{
			fSootIntegrals_.reset(new std::ofstream("soot.integrals", std::ios::out));
			fSootIntegrals_().setf(std::ios::scientific);

			fSootIntegrals_() << std::setw(20) << "time";
			fSootIntegrals_() << std::setw(20) << "soot_mass[kg]";
			fSootIntegrals_() << std::setw(20) << "soot_volume[m3]";
			fSootIntegrals_() << std::setw(20) << "soot_particles[#]";		
			fSootIntegrals_() << std::setw(20) << "prec_mass[kg]";
			fSootIntegrals_() << std::setw(20) << "prec_volume[m3]";
			fSootIntegrals_() << std::setw(20) << "prec_particles[#]";
			fSootIntegrals_() << std::setw(20) << "dimer_mass[kg]";		
			fSootIntegrals_() << std::setw(20) << "pah12_340nm[kg]";
			fSootIntegrals_() << std::setw(20) << "pah34_400nm[kg]";
			fSootIntegrals_() << std::setw(20) << "pahlarge_500nm[kg]";
			fSootIntegrals_() << std::endl;
		}
		else
		{
			fSootIntegrals_.reset(new std::ofstream("soot.integrals", std::ios::app));
			fSootIntegrals_().setf(std::ios::scientific);
		}
	}
}

void Foam::PolimiSootModel::UpdateProperties(const volScalarField& p, const volScalarField& T, const PtrList<volScalarField>& Y)
{
	label ns = Y.size();

	Eigen::VectorXd y(ns);
	Eigen::VectorXd x(ns);

	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();

	scalarField& soot_fv = sootFields_[0].ref();
	
	forAll(TCells, celli)
	{	
		double mw;

		// Extract the mean mass fractions
		for(unsigned int i=0;i<ns;i++)
			y(i) = Y[i].internalField()[celli];
		const double sum = y.sum();
		y /= sum;
	
		// Molecular weight and mole fractions
		thermodynamicsMap_.MoleFractions_From_MassFractions(x.data(), mw, y.data());

		// Density
		const double rhoGas = pCells[celli]*mw/PhysicalConstants::R_J_kmol/TCells[celli];

		// Soot analysis
		sootAnalyzer_->Analysis(TCells[celli], pCells[celli], rhoGas, y, x);
		soot_fv[celli] = sootAnalyzer_->fv_large();
	}

	forAll(T.boundaryField(), patchi)
	{
		const fvPatchScalarField& pT  = T.boundaryField()[patchi];
		const fvPatchScalarField& pp  = p.boundaryField()[patchi];

		fvPatchScalarField& psoot_fv = sootFields_[0].boundaryFieldRef()[patchi];
	
		forAll(pT, facei)
		{
			double mw;

			// Extract the mass fractions
			for(unsigned int i=0;i<ns;i++)
				y(i) = Y[i].boundaryField()[patchi][facei];
			const double sum = y.sum();
			y /= sum;

			// Molecular weight
			thermodynamicsMap_.MoleFractions_From_MassFractions(x.data(), mw, y.data());

			// Density
			const double rhoGas = pp[facei]*mw/PhysicalConstants::R_J_kmol/pT[facei];

			// Soot analysis
			sootAnalyzer_->Analysis(pT[facei], pp[facei], rhoGas, y, x);	
			psoot_fv[facei] = sootAnalyzer_->fv_large();
		}
	}
}

void Foam::PolimiSootModel::OnTheFlyAnalysis(const Foam::fvMesh& mesh, const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	if (ontheflyAnalysis_ == true && frequencyCounter_ == frequency_)
	{
		const double large_bin_mass = largeBinMass(rho, Y);
		const double small_bin_mass = smallBinMass(rho, Y);
		const double large_bin_volume = largeBinVolume(rho, Y);
		const double small_bin_volume = smallBinVolume(rho, Y);
		const double large_bin_particles = largeBinNumberParticles(rho, Y);
		const double small_bin_particles = smallBinNumberParticles(rho, Y);
		const double soot_dimer_mass = sootDimerMass(rho, Y);
		const double pah_1_2_rings_mass = pah_1_2_rings_Mass(rho, Y);
		const double pah_3_4_rings_mass = pah_3_4_rings_Mass(rho, Y);
		const double pah_more_than_4_rings_mass = pah_more_than_4_rings_Mass(rho, Y);

		// Only if the processor is the master, write to the file
		if (Pstream::master())
		{
	    		fSootIntegrals_() << std::setw(20) << mesh.time().value();
			fSootIntegrals_() << std::setw(20) << large_bin_mass;
			fSootIntegrals_() << std::setw(20) << large_bin_volume;
			fSootIntegrals_() << std::setw(20) << large_bin_particles;		
			fSootIntegrals_() << std::setw(20) << small_bin_mass;
			fSootIntegrals_() << std::setw(20) << small_bin_volume;
			fSootIntegrals_() << std::setw(20) << small_bin_particles;
			fSootIntegrals_() << std::setw(20) << soot_dimer_mass;		
			fSootIntegrals_() << std::setw(20) << pah_1_2_rings_mass;
			fSootIntegrals_() << std::setw(20) << pah_3_4_rings_mass;
			fSootIntegrals_() << std::setw(20) << pah_more_than_4_rings_mass;
			fSootIntegrals_() << std::endl;
		}
		
		Info << endl;
		Info << "Polimi Soot Model     " << endl;
		Info << "Soot mass:            " << large_bin_mass << endl;
		Info << "Soot volume:          " << large_bin_volume << endl;
		Info << "Soot particles:       " << large_bin_particles << endl;		
		Info << "Precursors mass:      " << small_bin_mass << endl;
		Info << "Precursors volume:    " << small_bin_volume << endl;
		Info << "Precursors particles: " << small_bin_particles << endl;
		Info << "Dimer mass:           " << soot_dimer_mass << endl;		
		Info << "PAH(340) mass:        " << pah_1_2_rings_mass << endl;
		Info << "PAH(400) mass:        " << pah_3_4_rings_mass << endl;
		Info << "PAH(500) mass:        " << pah_more_than_4_rings_mass << endl;
		Info << endl;

		frequencyCounter_ = 0;
	}

	frequencyCounter_++;
}

void Foam::PolimiSootModel::FieldsReconstruction(const Foam::fvMesh& mesh, const volScalarField& p, const volScalarField& T, const PtrList<volScalarField>& Y)
{
	if (mesh.time().outputTime() == true)
	{
		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

		label ns = Y.size();

		Eigen::VectorXd y(ns);
		Eigen::VectorXd x(ns);
		Eigen::VectorXd c(thermodynamicsMap_.NumberOfSpecies());
		Eigen::VectorXd R(thermodynamicsMap_.NumberOfSpecies());

		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField();

		scalarField& soot_fv = sootFields_[0].ref();
		scalarField& soot_rho = sootFields_[1].ref();
		scalarField& soot_N = sootFields_[2].ref();
		scalarField& soot_x = sootFields_[3].ref();
		scalarField& soot_omega = sootFields_[4].ref();
		scalarField& R_soot = sootFields_[5].ref();

		scalarField& pah_Y12 = pahFields_[0].ref();
		scalarField& pah_Y34 = pahFields_[1].ref();
		scalarField& pah_YLP = pahFields_[2].ref();

		
		forAll(TCells, celli)
		{	
			double mw;

			// Extract the mean mass fractions
			for(unsigned int i=0;i<ns;i++)
				y(i) = Y[i].internalField()[celli];
			const double sum = y.sum();
			y /= sum;
		
			// Molecular weight and mole fractions
			thermodynamicsMap_.MoleFractions_From_MassFractions(x.data(), mw, y.data());

			// Density
			const double rhoGas = pCells[celli]*mw/PhysicalConstants::R_J_kmol/TCells[celli];

			// Soot analysis
			sootAnalyzer_->Analysis(TCells[celli], pCells[celli], rhoGas, y, x);

			soot_fv[celli] = sootAnalyzer_->fv_large();
			soot_rho[celli] = sootAnalyzer_->rho_large();
			soot_N[celli] = sootAnalyzer_->N_large();
			soot_omega[celli] = sootAnalyzer_->omega_large();
			soot_x[celli] = sootAnalyzer_->x_large();

			pah_Y12[celli] = sootAnalyzer_->omega_pah_1_2_rings();
			pah_Y34[celli] = sootAnalyzer_->omega_pah_3_4_rings();
			pah_YLP[celli] = sootAnalyzer_->omega_pah_large();

			// Concentrations
			const double cTot = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
			c = cTot*x;

			// Kinetics
			kineticsMap_.SetTemperature(TCells[celli]);
			kineticsMap_.SetPressure(pCells[celli]);
			kineticsMap_.KineticConstants();
			kineticsMap_.ReactionRates(c.data());
			kineticsMap_.FormationRates(R.data());
			
			for(unsigned int i=0;i<ns;i++)
				R(i) *= thermodynamicsMap_.MW(i);
			R_soot[celli] = largeBinMassFormationRates(R);
		}

		forAll(T.boundaryField(), patchi)
		{
			const fvPatchScalarField& pT  = T.boundaryField()[patchi];
			const fvPatchScalarField& pp  = p.boundaryField()[patchi];

			fvPatchScalarField& psoot_fv = sootFields_[0].boundaryFieldRef()[patchi];
			fvPatchScalarField& psoot_rho = sootFields_[1].boundaryFieldRef()[patchi];
			fvPatchScalarField& psoot_N = sootFields_[2].boundaryFieldRef()[patchi];
			fvPatchScalarField& psoot_x = sootFields_[3].boundaryFieldRef()[patchi];
			fvPatchScalarField& psoot_omega = sootFields_[4].boundaryFieldRef()[patchi];
			fvPatchScalarField& pR_soot = sootFields_[5].boundaryFieldRef()[patchi];

			fvPatchScalarField& ppah_Y12 = pahFields_[0].boundaryFieldRef()[patchi];
			fvPatchScalarField& ppah_Y34 = pahFields_[1].boundaryFieldRef()[patchi];
			fvPatchScalarField& ppah_YLP = pahFields_[2].boundaryFieldRef()[patchi];
		
			forAll(pT, facei)
			{
				double mw;

				// Extract the mass fractions
				for(unsigned int i=0;i<ns;i++)
					y(i) = Y[i].boundaryField()[patchi][facei];
				const double sum = y.sum();
				y /= sum;

				// Molecular weight
				thermodynamicsMap_.MoleFractions_From_MassFractions(x.data(), mw, y.data());

				// Density
				const double rhoGas = pp[facei]*mw/PhysicalConstants::R_J_kmol/pT[facei];

				// Soot analysis
				sootAnalyzer_->Analysis(pT[facei], pp[facei], rhoGas, y, x);	
				psoot_fv[facei] = sootAnalyzer_->fv_large();
				psoot_rho[facei] = sootAnalyzer_->rho_large();
				psoot_N[facei] = sootAnalyzer_->N_large();
				psoot_omega[facei] = sootAnalyzer_->omega_large();
				psoot_x[facei] = sootAnalyzer_->x_large();

				ppah_Y12[facei] = sootAnalyzer_->omega_pah_1_2_rings();
				ppah_Y34[facei] = sootAnalyzer_->omega_pah_3_4_rings();
				ppah_YLP[facei] = sootAnalyzer_->omega_pah_large();

				// Concentrations
				const double cTot = pp[facei]/PhysicalConstants::R_J_kmol/pT[facei];
				c = cTot*x;

				// Kinetics
				kineticsMap_.SetTemperature(pT[facei]);
				kineticsMap_.SetPressure(pp[facei]);
				kineticsMap_.KineticConstants();
				kineticsMap_.ReactionRates(c.data());
				kineticsMap_.FormationRates(R.data());

				for(unsigned int i=0;i<ns;i++)
				R(i) *= thermodynamicsMap_.MW(i);

				pR_soot[facei] = largeBinMassFormationRates(R);
			}
		}
	}
}

// Returns the integral mass of soot particles (large BINS) [kg] 
double Foam::PolimiSootModel::largeBinMass(const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->bin_indices_large_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->bin_indices_large_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]).value();
	}
	return sum;
}

// Returns the integral volume of soot particles (large BINs) [m3] 
double Foam::PolimiSootModel::largeBinVolume(const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->bin_indices_large_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->bin_indices_large_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]/sootAnalyzer_->bin_density_large()[i]).value();
	}
	return sum;
}

// Returns the integral number of soot particles (large BINs) [#] 
double Foam::PolimiSootModel::largeBinNumberParticles(const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->bin_indices_large_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->bin_indices_large_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]/sootAnalyzer_->bin_density_large()[i]/sootAnalyzer_->bin_V_large()[i]).value();
	}
	return sum;
}

// Returns the integral mass of soot particles (small BINS) [kg] 
double Foam::PolimiSootModel::smallBinMass(const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->bin_indices_small_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->bin_indices_small_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]).value();
	}
	return sum;
}

// Returns the integral volume of soot particles (small BINs) [m3] 
double Foam::PolimiSootModel::smallBinVolume(const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->bin_indices_small_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->bin_indices_small_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]/sootAnalyzer_->bin_density_small()[i]).value();
	}
	return sum;
}

// Returns the integral number of soot particles (small BINs) [#] 
double Foam::PolimiSootModel::smallBinNumberParticles(const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->bin_indices_small_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->bin_indices_small_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]/sootAnalyzer_->bin_density_large()[i]/sootAnalyzer_->bin_V_small()[i]).value();
	}
	return sum;
}

// Returns the integral mass of dimer particles [kg] 
double Foam::PolimiSootModel::sootDimerMass(const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->soot_dimer_indices_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->soot_dimer_indices_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]).value();
	}
	return sum;
}

// Returns the integral mass of PAHs with 1 or 2 aromatic rings (340 nm) [kg] 
double Foam::PolimiSootModel::pah_1_2_rings_Mass(const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->pah_1_2_rings_indices_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->pah_1_2_rings_indices_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]).value();
	}
	return sum;
}

// Returns the integral mass of PAHs with 3 or 4 aromatic rings (400 nm) [kg] 
double Foam::PolimiSootModel::pah_3_4_rings_Mass(const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->pah_3_4_rings_indices_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->pah_3_4_rings_indices_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]).value();
	}
	return sum;
}

// Returns the integral mass of PAHs with more than 4 rings (500 nm) [kg] 
double Foam::PolimiSootModel::pah_more_than_4_rings_Mass(const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	return smallBinMass(rho, Y);
}

// Returns the total formation rate of soot particles (large BINs) [kg/m3/s]
double Foam::PolimiSootModel::largeBinMassFormationRates(const Eigen::VectorXd& R_times_W)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->bin_indices_large_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->bin_indices_large_global()[i];
		sum += R_times_W(j);
	}
	return sum;
}

// Returns the total formation rate of soot particles (small BINs) [kg/m3/s]
double Foam::PolimiSootModel::smallBinMassFormationRates(const Eigen::VectorXd& R_times_W)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->bin_indices_small_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->bin_indices_small_global()[i];
		sum += R_times_W(j);
	}
	return sum;
}

// Returns the total formation rate of PAHs with 1 or 2 aromatic rings (340 nm) [kg/m3/s] 
double Foam::PolimiSootModel::pah12RingsFormationRates(const Eigen::VectorXd& R_times_W)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->pah_1_2_rings_indices_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->pah_1_2_rings_indices_global()[i];
		sum += R_times_W(j);
	}
	return sum;
}

// Returns the total formation rate of PAHs with 3 or 4 aromatic rings (400 nm) [kg/m3/s] 
double Foam::PolimiSootModel::pah34RingsFormationRates(const Eigen::VectorXd& R_times_W)
{
	double sum = 0.;
	for(int i=0;i<sootAnalyzer_->pah_3_4_rings_indices_global().size();i++)
	{
		const unsigned int j = sootAnalyzer_->pah_3_4_rings_indices_global()[i];
		sum += R_times_W(j);
	}
	return sum;
}

// Returns the total formation rate of PAHs with more that aromatic rings (500 nm) [kg/m3/s] 
double Foam::PolimiSootModel::pahMoreThan4RingsFormationRates(const Eigen::VectorXd& R_times_W)
{
	largeBinMassFormationRates(R_times_W);
}



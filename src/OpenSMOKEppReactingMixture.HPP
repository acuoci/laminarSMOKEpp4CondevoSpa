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

// OpenSMOKE++ main class
#include "OpenSMOKEppReactingMixture.H"

// OpenFOAM boundary conditions
#include "zeroGradientFvPatchFields.H"
#include "fixedJumpFvPatchFields.H"
#include "fixedJumpAMIFvPatchFields.H"

// OpenSMOKE++ boundary conditions
#include "patches/gradientEnthalpy/gradientEnthalpyFvPatchScalarField.H"
#include "patches/fixedEnthalpy/fixedEnthalpyFvPatchScalarField.H"
#include "patches/mixedEnthalpy/mixedEnthalpyFvPatchScalarField.H"
#include "patches/enthalpyJump/enthalpyJumpFvPatchScalarField.H"
#include "patches/enthalpyJumpAMI/enthalpyJumpAMIFvPatchScalarField.H"

// Customized boundary conditions for multiregion solvers (only for OF>=8)
#if MULTIREGIONSOLVER == 1
#if WM_PROJECT_VERSION >= 8
#include "patches/temperatureCoupledBaseOpenSMOKEpp/temperatureCoupledBaseOpenSMOKEpp.H"
#include "patches/externalWallHeatFluxTemperatureOpenSMOKEpp/externalWallHeatFluxTemperatureOpenSMOKEppFvPatchScalarField.H"
#include "patches/turbulentTemperatureCoupledBaffleMixedOpenSMOKEpp/turbulentTemperatureCoupledBaffleMixedOpenSMOKEppFvPatchScalarField.H"
#endif
#endif

// Default name of main dictionary
const Foam::word Foam::OpenSMOKEppReactingMixture::dictName("SimulationOptions");

// Constructor
Foam::OpenSMOKEppReactingMixture::OpenSMOKEppReactingMixture(const fvMesh& mesh)
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

	T_
	(
		IOobject
		(
			"T",
			mesh.time().timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	),

	p_
	(
		IOobject
		(
			"p",
			mesh.time().timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	),

	psi_
	(
		IOobject
		(
			"mix:psi",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("mix:psi", dimensionSet(0, -2, 2, 0, 0), 0.)
	),

	Cp_
	(
		IOobject
		(
			"mix:Cp",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("mix:Cp", dimensionSet(0, 2, -2, -1, 0), 0.)
	),

	Cpv_
	(
		IOobject
		(
			"mix:Cpv",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("mix:Cpv", dimensionSet(0, 2, -2, -1, 0), 0.)
	),

	Hs_
	(
		IOobject
		(
			"mix:Hs",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("mix:Hs", dimensionSet(0, 2, -2, 0, 0), 0.),
		this->hBoundaryTypes(),
		this->hBoundaryBaseTypes()
	),

	mu_
	(
		IOobject
		(
			"mix:mu",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("mix:mu", dimensionSet(1, -1, -1, 0, 0), 0.)
	),

	lambda_
	(
		IOobject
		(
			"mix:lambda",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("mix:lambda", dimensionSet(1, 1, -3, -1, 0), 0.)
	),

	alpha_
	(
		IOobject
		(
			"mix:alpha",
			mesh.time().timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("mix:alpha", dimensionSet(0, 2, -1, 0, 0), 0.)
	),

	MW_
	(
		IOobject
		(
			"mix:MW",
			mesh.time().timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("mix:MW", dimensionSet(1, 0, 0, 0, -1), 0.)
	),

	cTot_
	(
		IOobject
		(
			"mix:cTot",
			mesh.time().timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("mix:cTot", dimensionSet(0, -3, 0, 0, 1), 0.)
	),

	Jc_
	(
	    IOobject
	    (
		"mix:Jc",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedVector("mix:Jc", dimensionSet(1, -2, -1, 0, 0), vector(0.0, 0.0, 0.0) )
	),

	phic_
	(
		IOobject
		(
			"phic",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("phic", dimensionSet(1, 0, -1, 0, 0), 0.)
	),

	Cps_Js_
	(
		IOobject
		(
			"mix:Cps_Js",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("mix:Cps_Js", dimensionSet(1, 0, -3, -1, 0), vector( 0.,0.,0.) )
	),

	hs_Js_
	(
		IOobject
		(
			"mix:hs_Js",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("mix:hs_Js", dimensionSet(1, 0, -3, 0, 0), vector( 0.,0.,0.) )
	),

	HsConduction_
	(
		IOobject
		(
			"mix:HsConduction",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("mix:HsConduction", dimensionSet(1, 0, -3, 0, 0), vector(0., 0., 0.))
	),

	Q_
	(
		IOobject
		(
			"mix:Q",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("mix:Q", dimensionSet(1, -1, -3, 0, 0), 0.0)
	),

	Qrad_
	(
		IOobject
		(
			"mix:Qrad",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("mix:Qrad", dimensionSet(1, -1, -3, 0, 0), 0.0)
	),

	cpuChemistry_
	(
		IOobject
		(
			"mix:cpuChem",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("mix:cpuChem", dimensionSet(0, 0, 1, 0, 0), 0.0)
	),

	soret_effect_(lookupOrDefault<Switch>("soretEffect", true)),

	thermophoretic_effect_(lookupOrDefault<Switch>("thermophoreticEffect", true)),

	Tmin_(lookupOrDefault<scalar>("Tmin", 270.)),

	Tmax_(lookupOrDefault<scalar>("Tmax", 3500.)),

	diffusion_via_mole_fraction_gradient_(lookupOrDefault<Switch>("moleFractionGradientBasedDiffusion", true)), 

	diffusion_fluxes_correction_(lookupOrDefault<Switch>("diffusionFluxCorrection", true)),

	closure_on_inert_(lookupOrDefault<Switch>("closureOnInert", true)),

	diffusion_in_energy_equation_(lookupOrDefault<Switch>("diffusionInEnergyEquation", true)),

	conduction_via_temperature_gradient_(lookupOrDefault<Switch>("temperatureGradientBasedConductiveFlux", true)),
	
	viscous_heating_source_term_(lookupOrDefault<Switch>("viscousHeatingSourceTerm", false)),

	solveForTauEquation_(lookupOrDefault<Switch>("solveForTauEquation", false)),

	solveForCsiEquation_(lookupOrDefault<Switch>("solveForCsiEquation", false)),

	solveForMomentumEquation_(lookupOrDefault<Switch>("solveForMomentumEquation", true)),

	solveForEnergyEquation_(lookupOrDefault<Switch>("solveForEnergyEquation", true)),
	
	solveForSpeciesEquations_(lookupOrDefault<Switch>("solveForSpeciesEquations", true)),

	homogeneousReactions_(lookupOrDefault<Switch>("homogeneousReactions", true)),

	direct_integration_minimum_temperature_for_chemistry_(lookupOrDefault<scalar>("minTemperatureForChemistry", 270.)),

	hmom_(lookupOrDefault<Switch>("HMOM", false)),

	polimiSoot_(lookupOrDefault<Switch>("PolimiSoot", false)),

	drg_(lookupOrDefault<Switch>("DRG", false)),

	clustering_(lookupOrDefault<Switch>("Clustering", false)),

	repairing_(lookupOrDefault<Switch>("RepairSolution", false)),

	tio2_(lookupOrDefault<Switch>("TiO2", false))

{
	// Energy equation formulation
	{
		word option = lookupOrDefault<word>("energyEquation", "temperature");
		if (option == "temperature")	energyEquation_ = ENERGY_EQ_TEMPERATURE;
		else if (option == "enthalpy")	energyEquation_ = ENERGY_EQ_ENTHALPY;
		else FatalError << "Wrong energyEquation option. Allowed options: temperature (default) | enthalpy" << ::Foam::exit(FatalError); 
	}

	// Diffusivity model
	{
		word option = lookupOrDefault<word>("diffusivityModel", "mix-averaged");
		if (option == "mix-averaged")		diffusivityModel_ = DIFFUSIVITY_MODEL_MIX_AVERAGED;
		else if (option == "lewis-numbers")	diffusivityModel_ = DIFFUSIVITY_MODEL_LEWIS_NUMBERS;
		else FatalError << "Wrong diffusivityModel option. Allowed options: mix-averaged (default) | lewis-numbers" << ::Foam::exit(FatalError); 
	}

	// Mesh data
	{
		ncells_ = mesh.nCells();
	}


	// Import kinetic mechanism
	{
		const dictionary& kineticsDictionary = this->subDict("Kinetics");

		//- Kinetics folder
		Foam::string kinetics_folder = kineticsDictionary.lookup("folder");

		//- Reading homogeneous data
		Info<< " * Loading homogeneous kinetics...\n" << endl;
		{	
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

			// Import XML file
			boost::filesystem::path path_kinetics = kinetics_folder;
			boost::property_tree::ptree ptree;
    			boost::property_tree::read_xml( (path_kinetics / "kinetics.xml").c_str(), ptree );

			// Import maps
			thermodynamicsMap_ = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree); 
			transportMap_ = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree); 
			kineticsMap_ = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMap_, ptree); 

			// Bundling of diffusion coefficients (to save CPU time)
			if (diffusivityModel_ == DIFFUSIVITY_MODEL_MIX_AVERAGED)
			{
				const double diffusionCoefficientsBundling = lookupOrDefault<scalar>("diffusionCoefficientsBundling", 0.);
				if (diffusionCoefficientsBundling > 0.)
					transportMap_->ImportSpeciesBundlingFromXMLFile(ptree, diffusionCoefficientsBundling);
			}

			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
		}

		// Number of species and reactions
		ns_ = thermodynamicsMap_->NumberOfSpecies();
		nr_ = kineticsMap_->NumberOfReactions();

		// Names of species
		for (int i=0;i<ns_;++i)
			species_.append(thermodynamicsMap_->NamesOfSpecies()[i]);

		//- Inert species
		word inertSpecies(kineticsDictionary.lookup("inertSpecies"));
		if (!species().found(inertSpecies))
			FatalError << "Inert species " << inertSpecies << " not found among the available species " << species() << ::Foam::exit(FatalError);
		inertIndex_ = species_[inertSpecies];

		//- Polimi Soot Model
		if (polimiSoot_ == true)
		{
			polimiSootModel_.reset( new PolimiSootModel(mesh, *thermodynamicsMap_, *kineticsMap_, *transportMap_) );
			polimiSootModel_().Setup(mesh);
		}

		//- HMOM Model
		if (hmom_ == true)
		{
			hmomModel_.reset( new HMOMModel(mesh, *thermodynamicsMap_, *kineticsMap_, *transportMap_) );
			hmomModel_().Setup(mesh);
		}

		//- DRG Model
		if (drg_ == true)
		{
			drgModel_.reset( new DRG(mesh, *thermodynamicsMap_, *kineticsMap_) );
			drgModel_().Setup(mesh);
		}

		//- Clustering Model
		if (clustering_ == true)
		{
			clusteringModel_.reset( new Clustering(mesh) );
		//	clusteringModel_().Setup(mesh);
		}

		//- Repairing Model
		if (repairing_ == true)
		{
			repairingModel_.reset( new Repairing(mesh) );
			repairingModel_().Setup(thermodynamicsMap_->NamesOfSpecies());
		}

		//- TiO2 Model
		if (tio2_ == true)
		{
			tio2Model_.reset( new TiO2(mesh, *thermodynamicsMap_, *kineticsMap_) );
			tio2Model_().Setup(mesh);
		}

		//- Chemistry Linear Model
		#if STEADYSTATE == 1
		{
			chemistryLinearModel_.reset( new ChemistryLinearModel(mesh, *thermodynamicsMap_, *kineticsMap_) );
			if (chemistryLinearModel_().jacobianUpdate() != 1)
			{
				Jstored_.resize(mesh.nCells());
				for(int i=0;i<mesh.nCells();i++)
					Jstored_[i].resize(ns_+1);
			}
		}
		#endif

		// Allocate memory
		Info << " * Memory allocation (species mass/mole fractions, diffusion coefficients, specific heats, enthalpies, mass fluxes, ...)" << endl;
		Y_.resize(ns_);
		X_.resize(ns_);
		Dmix_.resize(ns_);
		CpSpecies_.resize(ns_);
		hsSpecies_.resize(ns_);
		J_.resize(ns_);

		// Save the formation enthalpies of species
		Info << " * Saving formation enthalpies of species..." << endl;
		h0_species_.resize(ns_);
		for(int i=0;i<ns_;++i)
		{
			// Set temperature and pressure for thermodynamic and transport maps
			thermodynamicsMap_->SetPressure(101325.);
			thermodynamicsMap_->SetTemperature(298.15);

			// Calculate the formation enthalpy [J/kg]
			thermodynamicsMap_->hMolar_Species(h0_species_.data());		// formation enthalpies of species [J/kmol]
			for(unsigned int i=0;i<ns_;++i)
		        	h0_species_(i) /= thermodynamicsMap_->MW(i);		// formation enthalpies of species [J/kg]
		}

		// Allocate meory for additional fields for correcting the diffusion fluxes
		if (diffusion_via_mole_fraction_gradient_ == true)
			mwDiffusionCorrection_.resize(ns_);

		// Allocate memory for Soret diffusion coefficients
		if (soret_effect_ == true)
		{
			soret_index_.resize(ns_, -1);
			soret_list_.resize(ns_, false);
			Dsoret_.resize(transportMap_->iThermalDiffusionRatios().size());
		}

		// Active species
		active_.resize(ns_, true);

		// Lewis numbers
		{
			lewis_numbers_.resize(ns_);
			lewis_numbers_.setConstant(1.);

			if (diffusivityModel_ == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
			{
				const dictionary& lewisNumberDictionary = this->subDict("LewisNumbers");

				// Default policy
				word default_policy = lewisNumberDictionary.lookup("default");
			
				// Read user-defined Lewis numbers
				for(int i=0;i<ns_;i++)
					lewis_numbers_(i) = lewisNumberDictionary.lookupOrDefault<scalar>(thermodynamicsMap_->NamesOfSpecies()[i], 0.);
			
				// Assign default values
				for(int i=0;i<ns_;i++)
					if (lewis_numbers_(i) == 0.)
					{
						if (default_policy == "mw-based")
							lewis_numbers_(i) = std::min( 0.15 + 0.0287*thermodynamicsMap_->MW(i), 50. );	// linear correlation with molecular weight (max 50)
						else
							lewis_numbers_(i) = std::atof(default_policy.c_str());				// constant default value
					}

				// Print Lewis numbers on the screen
				Info << " * Diffusivities of species will be calculated using the following Lewis numbers: " << endl;
				for(int i=0;i<ns_;i++)
					Info << "   " << thermodynamicsMap_->NamesOfSpecies()[i] << "\t" << lewis_numbers_(i) << endl;
			}
		}
		

		// Molecular weights of species as dimensioned scalars [kg/kmol]
		forAll(species_, i)
			mw_species_.append( dimensionedScalar("MW", dimensionSet(0,0,0,0,0,0,0), scalar(thermodynamicsMap_->MW(i))) );
		forAll(species_, i)
			mw_species_[i].dimensions().reset(dimensionSet(1,0,0,0,-1,0,0));

		// In case of closure on inert, set the inert as unactive species
		if (closure_on_inert_ == true)
			active_[inertIndex_] = false;

	}

	tmp<volScalarField> tYdefault;
	forAll(species_, i)
	{
		IOobject header
		(
			species_[i],
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		// Check if field exists and can be read
		if (header.typeHeaderOk<volScalarField>(true))
		{
			DebugInfo << "Reading " << species_[i] << endl;

			Y_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						species_[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);
		}
		else
		{
			// Read Ydefault if not already read
			if (!tYdefault.valid())
			{
				word YdefaultName("Ydefault");

				IOobject timeIO ( YdefaultName, mesh.time().timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE );
				IOobject constantIO ( YdefaultName,mesh.time().constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE );
				IOobject time0IO ( YdefaultName,Time::timeName(0),mesh,IOobject::MUST_READ,IOobject::NO_WRITE );

				if (timeIO.typeHeaderOk<volScalarField>(true))
				{
					tYdefault = new volScalarField(timeIO, mesh);
				}
				else if (constantIO.typeHeaderOk<volScalarField>(true))
				{
					tYdefault = new volScalarField(constantIO, mesh);
				}
				else
				{
					tYdefault = new volScalarField(time0IO, mesh);
				}
			}

			Y_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						species_[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					tYdefault()
				)
			);
		}
	}

	Info << " * Allocating memory for mixture fractions of species ..." << endl;
	{
		IOobject::writeOption  write_option = IOobject::NO_WRITE;
		if ( lookupOrDefault<Switch>("outputMoleFractions", false) == true)
			write_option = IOobject::AUTO_WRITE;

		forAll(species_, i)
		{
			X_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"X:" + species_[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						write_option
					),
					mesh,
					dimensionedScalar("X", dimensionSet(0, 0, 0, 0, 0), 0.)
				)
			);
		}
	}


	Info << " * Allocating memory for Dmix of species ..." << endl;
	{	
		IOobject::writeOption  write_option = IOobject::NO_WRITE;
		if ( lookupOrDefault<Switch>("outputDiffusionCoefficients", false) == true)
			write_option = IOobject::AUTO_WRITE;

		forAll(species_, i)
		{
			Dmix_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"mix:Dmix:" + species_[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						write_option
					),
					mesh,
					dimensionedScalar("mix:Dmix", dimensionSet(0, 2, -1, 0, 0), 0.)
				)
			);
		}
	}

	Info << " * Allocating memory for Cp of species ..." << endl;
	{
		IOobject::writeOption  write_option = IOobject::NO_WRITE;
		if ( lookupOrDefault<Switch>("outputSpeciesSpecificHeats", false) == true)
			write_option = IOobject::AUTO_WRITE;

		forAll(species_, i)
		{
			CpSpecies_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"mix:Cp:" + species_[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						write_option
					),
					mesh,
					dimensionedScalar("mix:Cp", dimensionSet(0, 2, -2, -1, 0), 0.)
				)
			);
		}
	}

	Info << " * Allocating memory for sensible enthalpy of species ..." << endl;
	{
		IOobject::writeOption  write_option = IOobject::NO_WRITE;
		if ( lookupOrDefault<Switch>("outputSpeciesSensibleEnthalpies", false) == true)
			write_option = IOobject::AUTO_WRITE;

		forAll(species_, i)
		{
			hsSpecies_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"mix:hs:" + species_[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						write_option
					),
					mesh,
					dimensionedScalar("mix:h", dimensionSet(0, 2, -2, 0, 0), 0.)
				)
			);
		}
	}

	Info << " * Allocating memory for fluxes J of species ..." << endl;
	{
		IOobject::writeOption  write_option = IOobject::NO_WRITE;
		if ( lookupOrDefault<Switch>("outputSpeciesFluxes", false) == true)
			write_option = IOobject::AUTO_WRITE;

		forAll (species_,i)
		{
			J_.set
			(
				i,
				new volVectorField
				(
					IOobject
					(
						"mix:J:" + species_[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						write_option
					),
					mesh,
					dimensionedVector("mix:J", dimensionSet(1, -2, -1, 0, 0), vector(0.,0.,0.))
				)
			);
		}
	}

	if (diffusion_via_mole_fraction_gradient_ == true)
	{
		Info << " * Allocating memory for MW correction on diffusion fluxes..." << endl;
		forAll (species_,i)
		{
			mwDiffusionCorrection_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"mix:mwDiffCorr:" + species_[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("mix:mwDiffCorr", dimensionSet(1, -3, -1, 0, 0), 0.)
				)
			);
		}
	}

	if (soret_effect_ == true)
	{
		IOobject::writeOption  write_option = IOobject::NO_WRITE;
		if ( lookupOrDefault<Switch>("outputSoretDiffusionCoefficients", false) == true)
			write_option = IOobject::AUTO_WRITE;

		for(unsigned int i=0;i<transportMap_->iThermalDiffusionRatios().size();i++)
		{ 
			unsigned int index = transportMap_->iThermalDiffusionRatios()[i]-1;
			soret_index_[index] = i;
			soret_list_[index] = true;

			Info << " * Soret effect: species " 	<< i+1 << "/" << transportMap_->iThermalDiffusionRatios().size() << " " 
								<< thermodynamicsMap_->NamesOfSpecies()[index] << " (" << index << ")" << endl;

			Dsoret_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"mix:Dsoret:"+ species_[index],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						write_option
					),
					mesh,
					dimensionedScalar("mix:Dsoret", dimensionSet(0, 2, -1, 0, 0), 0.)
				)
			);
		}
	}

	#if STEADYSTATE == 1
	Info << " * Allocating memory for implicit/explicit chemical source terms..." << endl;
	{
		sourceImplicit_.resize(ns_+1);
		sourceExplicit_.resize(ns_+1);

		forAll (species_,i)
		{
			sourceImplicit_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"mix:SI_" + species_[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("mix:SI", dimensionSet(1, -3, -1, 0, 0), 0.0)
				)
			);

			sourceExplicit_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"mix:SE_" + species_[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("mix:SE", dimensionSet(1, -3, -1, 0, 0), 0.0)
				)
			);
		}

		sourceImplicit_.set
		(
			ns_,
			new volScalarField
			(
				IOobject
				(
					"mix:SI_T",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensionedScalar("mix:SI_T", dimensionSet(1, -1, -3, -1, 0), 0.0)
			)
		);

		sourceExplicit_.set
		(
			ns_,
			new volScalarField
			(
				IOobject
				(
					"mix:SE_T",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensionedScalar("mix:SE_T", dimensionSet(1, -1, -3, 0, 0), 0.0)
			)
		);
	}
	#endif

	// Initialize the enthalpy field
	init();

	// Calculate the properties
	calculate();

	// Switch on saving old time
	this->psi_.oldTime();

	// Read ODE solver options
	if (homogeneousReactions_ == true)
	{
		// Batch reactor homogeneous: ode parameters
		const dictionary& odeHomogeneousDictionary = this->subDict("OdeHomogeneous");

		//- Mass fractions tolerance
		scalar relTolerance = readScalar(odeHomogeneousDictionary.lookup("relTolerance"));
		scalar absTolerance = readScalar(odeHomogeneousDictionary.lookup("absTolerance"));
		odeParameterBatchReactorHomogeneous_.SetRelativeTolerance(relTolerance);
		odeParameterBatchReactorHomogeneous_.SetAbsoluteTolerance(absTolerance);

		//- Full pivoting (only for OpenSMOKE++ solver)
		Switch fullPivoting(odeHomogeneousDictionary.lookup("fullPivoting"));
		odeParameterBatchReactorHomogeneous_.SetFullPivoting(fullPivoting);

		//- Maximum order of integration (only for OpenSMOKE++ solver)
		label maximumOrder = readLabel(odeHomogeneousDictionary.lookup("maximumOrder"));
		odeParameterBatchReactorHomogeneous_.SetMaximumOrder(maximumOrder);

		// Type
		word homogeneousODESolverString(odeHomogeneousDictionary.lookup("odeSolver"));
		if ( homogeneousODESolverString != "OpenSMOKE++" &&  homogeneousODESolverString != "CVODE" )
			FatalError << "Wrong homogeneous ODE Solver: OpenSMOKE++ || CVODE" << ::Foam::exit(FatalError);

		// OpenSMOKE++ Ode solver
		if (homogeneousODESolverString == "OpenSMOKE++")
			odeParameterBatchReactorHomogeneous_.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE);
		
		// CVODE Ode solver
		if (homogeneousODESolverString == "CVODE") 	
		{
			odeParameterBatchReactorHomogeneous_.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE);
			#if OPENSMOKE_USE_SUNDIALS != 1
				FatalError << "The solver was compiled without the CVODE support. Please choose a different ODE solver." << ::Foam::exit(FatalError);
			#endif
		}
	}
	
	// Allocate memory for chemistry
	if (homogeneousReactions_ == true && hmom_ == false && drg_ == false && tio2_ == false)
	{
		// Batch reactor homogeneous
		batchReactorHomogeneousConstantPressure_ = new BatchReactorHomogeneousConstantPressure(*thermodynamicsMap_, *kineticsMap_);

		// ODE Solver (constant pressure)
		odeSolverConstantPressure_ = new OdeSMOKE::MultiValueSolver<methodGearConstantPressure>;
		odeSolverConstantPressure_->SetReactor(batchReactorHomogeneousConstantPressure_);

		// ODE solver (constant pressure, SUNDIALS CVODE)
		#if OPENSMOKE_USE_SUNDIALS == 1
		if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
		{
			cvodeSolverConstantPressure_ = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_CVODE::GetInstance();
			cvodeSolverConstantPressure_->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure_);
		}
		#endif
	}
	else if (homogeneousReactions_ == true && hmom_ == false && drg_ == true && tio2_ == false)
	{
		// Batch reactor homogeneous
		batchReactorHomogeneousConstantPressure_DRG_ = new BatchReactorHomogeneousConstantPressure_DRG( *thermodynamicsMap_, *kineticsMap_ );

		// ODE Solver (constant pressure)
		odeSolverConstantPressure_DRG_ = new OdeSMOKE::MultiValueSolver<methodGearConstantPressure_DRG>;
		odeSolverConstantPressure_DRG_->SetReactor(batchReactorHomogeneousConstantPressure_DRG_);

		// ODE solver (constant pressure, SUNDIALS CVODE)
		#if OPENSMOKE_USE_SUNDIALS == 1
		if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
		{
			cvodeSolverConstantPressure_DRG_ = OpenSMOKE::BatchReactorHomogeneousConstantPressure_DRG_ODE_CVODE::GetInstance();
			cvodeSolverConstantPressure_DRG_->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure_DRG_);
		}
		#endif
	}
	else if (homogeneousReactions_ == true && hmom_ == true && drg_ == false && tio2_ == false)
	{
		// Batch reactor homogeneous
		batchReactorHomogeneousConstantPressure_HMOM_ = new BatchReactorHomogeneousConstantPressure_HMOM(*thermodynamicsMap_, *kineticsMap_, hmomModel_().hmomAnalyzer());

		// ODE Solver (constant pressure)
		odeSolverConstantPressure_HMOM_ = new OdeSMOKE::MultiValueSolver<methodGearConstantPressure_HMOM>;
		odeSolverConstantPressure_HMOM_->SetReactor(batchReactorHomogeneousConstantPressure_HMOM_);

		// ODE solver (constant pressure, SUNDIALS CVODE)
		#if OPENSMOKE_USE_SUNDIALS == 1
		if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
		{
			cvodeSolverConstantPressure_HMOM_ = OpenSMOKE::BatchReactorHomogeneousConstantPressure_HMOM_ODE_CVODE::GetInstance();
			cvodeSolverConstantPressure_HMOM_->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure_HMOM_);
		}
		#endif
	}
	else if (homogeneousReactions_ == true && hmom_ == false && drg_ == false && tio2_ == true)
	{
		// Batch reactor homogeneous
		batchReactorHomogeneousConstantPressure_TIO2_ = new BatchReactorHomogeneousConstantPressure_TIO2(*thermodynamicsMap_, *kineticsMap_, tio2Model_().tio2Model());

		// ODE Solver (constant pressure)
		odeSolverConstantPressure_TIO2_ = new OdeSMOKE::MultiValueSolver<methodGearConstantPressure_TIO2>;
		odeSolverConstantPressure_TIO2_->SetReactor(batchReactorHomogeneousConstantPressure_TIO2_);

		// ODE solver (constant pressure, SUNDIALS CVODE)
		#if OPENSMOKE_USE_SUNDIALS == 1
		if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
		{
			cvodeSolverConstantPressure_TIO2_ = OpenSMOKE::BatchReactorHomogeneousConstantPressure_TIO2_ODE_CVODE::GetInstance();
			cvodeSolverConstantPressure_TIO2_->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure_TIO2_);
		}
		#endif
	}


	/*
	// Additional output
	{
		// Formation rates of selected species
		{
			List<word>  formationRates(outputDictionary.lookup("formationRates"));
			outputFormationRatesIndices_.resize(formationRates.size());
			if (formationRates.size() == 1)
			{
				if ( formationRates[0] == "all" || formationRates[0] == "ALL" )
				{
					outputFormationRatesIndices_.resize(ns_);
					for (int i=0;i<ns_;i++)
						outputFormationRatesIndices_(i) = i;
				}
				else
					outputFormationRatesIndices_(0) = thermodynamicsMap_->IndexOfSpecies(formationRates[0])-1;
			}
			else
			{
				for (int i=0;i<formationRates.size();i++)
					outputFormationRatesIndices_(i) = thermodynamicsMap_->IndexOfSpecies(formationRates[i])-1;
			}
		}
	}
	*/

	// Final checks
	if (polimiSoot_ == true)
		if ( (polimiSootModel_().physicalDiffusion()) && (diffusivityModel_ != DIFFUSIVITY_MODEL_MIX_AVERAGED) )
			FatalError << "The soot physical diffusion model can be used only together with a mix-averaged diffusivity model. Please check your options." << ::Foam::exit(FatalError);
}


Foam::autoPtr<Foam::OpenSMOKEppReactingMixture> Foam::OpenSMOKEppReactingMixture::New(const fvMesh& mesh)
{
	return Foam::autoPtr<Foam::OpenSMOKEppReactingMixture>( new OpenSMOKEppReactingMixture(mesh));
}


Foam::OpenSMOKEppReactingMixture::~OpenSMOKEppReactingMixture()
{}


Foam::tmp<Foam::volScalarField> Foam::OpenSMOKEppReactingMixture::rho() const
{
	return p_*psi_;
}

Foam::tmp<Foam::scalarField> Foam::OpenSMOKEppReactingMixture::rho(const label patchi) const
{
	return p_.boundaryField()[patchi]*psi_.boundaryField()[patchi];
}

Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::p()
{
	return p_;
}

const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::p() const
{
	return p_;
}

const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::T() const
{
	return T_;
}

Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::T()
{
	return T_;
}

const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::psi() const
{
	return psi_;
}

Foam::tmp<Foam::volScalarField> Foam::OpenSMOKEppReactingMixture::Cp() const
{
	return Cp_;
}

Foam::tmp<Foam::scalarField> Foam::OpenSMOKEppReactingMixture::Cp(const label patchi) const
{
	return Cp_.boundaryField()[patchi];
}

Foam::tmp<Foam::volScalarField> Foam::OpenSMOKEppReactingMixture::Cpv() const
{
	return Cpv_;
}

Foam::tmp<Foam::scalarField> Foam::OpenSMOKEppReactingMixture::Cpv(const label patchi) const
{
	return Cpv_.boundaryField()[patchi];
}

const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::Hs() const
{
	return Hs_;
}

Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::Hs()
{
	return Hs_;
}

Foam::tmp<Foam::scalarField> Foam::OpenSMOKEppReactingMixture::Hs(const label patchi) const
{
	return Hs_.boundaryField()[patchi];
}

Foam::tmp<Foam::volScalarField> Foam::OpenSMOKEppReactingMixture::mu() const
{
	return mu_;
}

Foam::tmp<Foam::scalarField> Foam::OpenSMOKEppReactingMixture::mu(const label patchi) const
{
	return mu_.boundaryField()[patchi];
}

Foam::tmp<Foam::volScalarField> Foam::OpenSMOKEppReactingMixture::lambda() const
{
	return lambda_;
}

Foam::tmp<Foam::scalarField> Foam::OpenSMOKEppReactingMixture::lambda(const label patchi) const
{
	return lambda_.boundaryField()[patchi];
}

const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::alpha() const
{
	return alpha_;
}

const Foam::scalarField& Foam::OpenSMOKEppReactingMixture::alpha(const label patchi) const
{
	return alpha_.boundaryField()[patchi];
}

const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::MW() const
{
	return MW_;
}

const Foam::scalarField& Foam::OpenSMOKEppReactingMixture::MW(const label patchi) const
{
	return MW_.boundaryField()[patchi];
}

const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::cTot() const
{
	return cTot_;
}

const Foam::scalarField& Foam::OpenSMOKEppReactingMixture::cTot(const label patchi) const
{
	return cTot_.boundaryField()[patchi];
}

Foam::volVectorField& Foam::OpenSMOKEppReactingMixture::Jc()
{
	return Jc_;
}

Foam::tmp<Foam::vectorField> Foam::OpenSMOKEppReactingMixture::Jc(const label patchi) const
{
	return Jc_.boundaryField()[patchi];
}

Foam::surfaceScalarField& Foam::OpenSMOKEppReactingMixture::phic()
{
	return phic_;
}

Foam::tmp<Foam::scalarField> Foam::OpenSMOKEppReactingMixture::phic(const label patchi) const
{
	return phic_.boundaryField()[patchi];
}

Foam::volVectorField& Foam::OpenSMOKEppReactingMixture::Cps_Js()
{
	return Cps_Js_;
}

Foam::tmp<Foam::vectorField> Foam::OpenSMOKEppReactingMixture::Cps_Js(const label patchi) const
{
	return Cps_Js_.boundaryField()[patchi];
}

Foam::volVectorField& Foam::OpenSMOKEppReactingMixture::hs_Js()
{
	return hs_Js_;
}

Foam::tmp<Foam::vectorField> Foam::OpenSMOKEppReactingMixture::hs_Js(const label patchi) const
{
	return hs_Js_.boundaryField()[patchi];
}

Foam::volVectorField& Foam::OpenSMOKEppReactingMixture::HsConduction()
{
	return HsConduction_;
}

Foam::tmp<Foam::vectorField> Foam::OpenSMOKEppReactingMixture::HsConduction(const label patchi) const
{
	return HsConduction_.boundaryField()[patchi];
}

const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::Q() const
{
	return Q_;
}

const Foam::scalarField& Foam::OpenSMOKEppReactingMixture::Q(const label patchi) const
{
	return Q_.boundaryField()[patchi];
}

Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::Qrad()
{
	return Qrad_;
}

const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::cpuChemistry() const
{
	return cpuChemistry_;
}

const Foam::scalarField& Foam::OpenSMOKEppReactingMixture::cpuChemistry(const label patchi) const
{
	return cpuChemistry_.boundaryField()[patchi];
}


void Foam::OpenSMOKEppReactingMixture::calculate()
{
	Info << " * Properties evaluation... ";

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	// Auxiliary variables
	Eigen::VectorXd mass_fractions(ns_);
	Eigen::VectorXd mole_fractions(ns_);
	Eigen::VectorXd cp_species(ns_);
	Eigen::VectorXd h_species(ns_);
	Eigen::VectorXd Dmix_species(ns_);
	Eigen::VectorXd tetamix_vector(ns_);

	// Const fields
	const scalarField& pCells = this->p_;

	// Non-const fields
	scalarField& TCells = this->T_.primitiveFieldRef();
	scalarField& psiCells = this->psi_.primitiveFieldRef();
	scalarField& cpCells = this->Cp_.primitiveFieldRef();
	scalarField& cpvCells = this->Cpv_.primitiveFieldRef();
	scalarField& hsCells = this->Hs_.primitiveFieldRef();
	scalarField& muCells = this->mu_.primitiveFieldRef();
	scalarField& lambdaCells = this->lambda_.primitiveFieldRef();
	scalarField& alphaCells = this->alpha_.primitiveFieldRef();
	scalarField& mwCells = this->MW_.primitiveFieldRef();
	scalarField& cTotCells = this->cTot_.primitiveFieldRef();


	// Loop over all the internal cells
	forAll(TCells, celli)
	{			
		// Extract the mass fractions
		for(int i=0;i<ns_;++i)
			mass_fractions(i) = Y_[i].internalField()[celli];

		// Calculate the mole fractions from the mass fractions
		thermodynamicsMap_->MoleFractions_From_MassFractions(mole_fractions.data(), mwCells[celli], mass_fractions.data());
		for(int i=0;i<ns_;++i)
			X_[i].ref()[celli] = mole_fractions(i);

		if (energyEquation_ == ENERGY_EQ_ENTHALPY)
		{
			// Calculate total enthalpy [J/kg]
			const double htot = hsCells[celli] + mass_fractions.dot(h0_species_);

			// Reconstruct the temperaure [K] from molar enthalpy [J/kmol]
			const double temp = thermodynamicsMap_->GetTemperatureFromEnthalpyAndMoleFractions(htot*mwCells[celli], pCells[celli], mole_fractions.data(), TCells[celli]);
			if (temp != 0.) 
			{
				TCells[celli] = temp;
				TCells[celli] = std::max(TCells[celli], Tmin_);
				TCells[celli] = std::min(TCells[celli], Tmax_);
			}

			// Set temperature and pressure for thermodynamic and transport maps
			thermodynamicsMap_->SetPressure(pCells[celli]);
			thermodynamicsMap_->SetTemperature(TCells[celli]);
			transportMap_->SetPressure(pCells[celli]);
			transportMap_->SetTemperature(TCells[celli]);
		}
		else if (energyEquation_ == ENERGY_EQ_TEMPERATURE)
		{
			// Set temperature and pressure for thermodynamic and transport maps
			thermodynamicsMap_->SetPressure(pCells[celli]);
			thermodynamicsMap_->SetTemperature(TCells[celli]);
			transportMap_->SetPressure(pCells[celli]);
			transportMap_->SetTemperature(TCells[celli]);

			// Calculate the enthalpy [J/kg]
			hsCells[celli]  = thermodynamicsMap_->hMolar_Mixture_From_MoleFractions(mole_fractions.data())/mwCells[celli];		// total enthalpy [J/kg]
			hsCells[celli] -= mass_fractions.dot(h0_species_);									// sensible enthalpy [J/kg]
		}
			
		// Calculate the total concentration [kmol/m3]
		cTotCells[celli] = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];

		// Calculate the compressibility [s2/m2]
		psiCells[celli]  = cTotCells[celli]*mwCells[celli]/pCells[celli];

		// Calculate the constant pressure specific heat [J/kg/K]
            	{
			cpCells[celli]   = thermodynamicsMap_->cpMolar_Mixture_From_MoleFractions(mole_fractions.data());	//[J/kmol/K]
			cpvCells[celli]  = cpCells[celli] - PhysicalConstants::R_J_kmol;					//[J/kmol/K]
			cpCells[celli]  /= mwCells[celli];
			cpvCells[celli] /= mwCells[celli];
	
			{
				thermodynamicsMap_->cpMolar_Species(cp_species.data());
				for(unsigned int i=0;i<ns_;++i)
		                       CpSpecies_[i].ref()[celli] = cp_species(i) / thermodynamicsMap_->MW(i);
			}

			{
				thermodynamicsMap_->hMolar_Species(h_species.data());						// total ethalpies of species [J/kmol]
				for(unsigned int i=0;i<ns_;++i)
		                       hsSpecies_[i].ref()[celli] = h_species(i)/thermodynamicsMap_->MW(i) - h0_species_(i);	// sensible enthalpies of species [J/kg]
			}
		}

		// Calculate the viscosity [kg/m/s]
		{
			muCells[celli] = transportMap_->DynamicViscosity(mole_fractions.data());
		}

		// Calculate the thermal conductivity [W/m/K]
            	{
			lambdaCells[celli] = transportMap_->ThermalConductivity(mole_fractions.data());
		}

		// Thermal diffusivity [m2/s]
		{
			alphaCells[celli] = lambdaCells[celli]/cpCells[celli]/(pCells[celli]*psiCells[celli]);
		}

		// Mass diffusion coefficients [m2/s]
		if (diffusivityModel_ == DIFFUSIVITY_MODEL_MIX_AVERAGED)
		{
			transportMap_->MassDiffusionCoefficients( Dmix_species.data(), mole_fractions.data(), transportMap_->is_species_bundling() );
			for(int i=0;i<ns_;++i)
				Dmix_[i].ref()[celli] = Dmix_species(i);
		}
		else if (diffusivityModel_ == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
		{
			for(int i=0;i<ns_;++i)
				Dmix_[i].ref()[celli] = alphaCells[celli]/lewis_numbers_(i);
		}

		// Thermal diffusion coefficients [-]
		if (soret_effect_ == true)
		{
			transportMap_->ThermalDiffusionRatios(tetamix_vector.data(), mole_fractions.data());		
			for(int i=0;i<transportMap_->iThermalDiffusionRatios().size();i++)
			{
				unsigned int index = transportMap_->iThermalDiffusionRatios()[i]-1;
				Dsoret_[i].ref()[celli] = Dmix_[index].internalField()[celli]*tetamix_vector(index)*thermodynamicsMap_->MW(index)/mwCells[celli];
			} 
		}
	}

	// Boundaries
	{
		volScalarField::Boundary& pBf 		= this->p_.boundaryFieldRef();
		volScalarField::Boundary& TBf 		= this->T_.boundaryFieldRef();
		volScalarField::Boundary& psiBf 	= this->psi_.boundaryFieldRef();
		volScalarField::Boundary& cpBf 		= this->Cp_.boundaryFieldRef();
		volScalarField::Boundary& cpvBf 	= this->Cpv_.boundaryFieldRef();
		volScalarField::Boundary& hsBf 		= this->Hs_.boundaryFieldRef();

		volScalarField::Boundary& muBf 		= this->mu_.boundaryFieldRef();
		volScalarField::Boundary& lambdaBf 	= this->lambda_.boundaryFieldRef();
		volScalarField::Boundary& alphaBf 	= this->alpha_.boundaryFieldRef();
		
		volScalarField::Boundary& mwBf 		= this->MW_.boundaryFieldRef();
		volScalarField::Boundary& cTotBf 	= this->cTot_.boundaryFieldRef();

		forAll(this->T_.boundaryField(), patchi)
		{
			fvPatchScalarField& pp = pBf[patchi];
			fvPatchScalarField& pT = TBf[patchi];
			fvPatchScalarField& ppsi = psiBf[patchi];
			fvPatchScalarField& pcp = cpBf[patchi];
			fvPatchScalarField& pcpv = cpvBf[patchi];
			fvPatchScalarField& phs = hsBf[patchi];

			fvPatchScalarField& pmu = muBf[patchi];
			fvPatchScalarField& plambda = lambdaBf[patchi];
			fvPatchScalarField& palpha = alphaBf[patchi];
			fvPatchScalarField& pmw = mwBf[patchi];
			fvPatchScalarField& pcTot = cTotBf[patchi];
	
			forAll(pT, facei)
			{
				// Calculate the mole fractions from the mass fractions
				for(int i=0;i<ns_;++i)
					mass_fractions(i) = Y_[i].boundaryField()[patchi][facei];
				thermodynamicsMap_->MoleFractions_From_MassFractions(mole_fractions.data(),pmw[facei],mass_fractions.data());
				for(int i=0;i<ns_;++i)
					X_[i].boundaryFieldRef()[patchi][facei] = mole_fractions(i);

				if (energyEquation_ == ENERGY_EQ_ENTHALPY)
				{
					if (pT.fixesValue())
					{
						// Set temperature and pressure for thermodynamic and transport maps
						thermodynamicsMap_->SetPressure(pp[facei]);
						thermodynamicsMap_->SetTemperature(pT[facei]);
						transportMap_->SetPressure(pp[facei]);
						transportMap_->SetTemperature(pT[facei]);

						// Calculate the enthalpy [J/kg]
						phs[facei]  = thermodynamicsMap_->hMolar_Mixture_From_MoleFractions(mole_fractions.data())/pmw[facei];	// total enthalpy [J/kg]
						phs[facei] -= mass_fractions.dot(h0_species_);								// sensible enthalpy [J/kg]
					}
					else
					{
						// Calculate total enthalpy [J/kg]
						const double htot = phs[facei] + mass_fractions.dot(h0_species_);

						// Reconstruct the temperaure [K] from molar enthalpy [J/kmol]
						const double temp = thermodynamicsMap_->GetTemperatureFromEnthalpyAndMoleFractions(htot*pmw[facei], pp[facei], mole_fractions.data(), pT[facei]);
						if (temp != 0.) 
						{
							pT[facei] = temp;
							pT[facei] = std::max(pT[facei], Tmin_);
							pT[facei] = std::min(pT[facei], Tmax_);
						}

						// Set temperature and pressure for thermodynamic and transport maps
						thermodynamicsMap_->SetPressure(pp[facei]);
						thermodynamicsMap_->SetTemperature(pT[facei]);
						transportMap_->SetPressure(pp[facei]);
						transportMap_->SetTemperature(pT[facei]);
					}
				}
				else if (energyEquation_ == ENERGY_EQ_TEMPERATURE)
				{
					// Set temperature and pressure for thermodynamic and transport maps
					thermodynamicsMap_->SetPressure(pp[facei]);
					thermodynamicsMap_->SetTemperature(pT[facei]);
					transportMap_->SetPressure(pp[facei]);
					transportMap_->SetTemperature(pT[facei]);

					// Calculate the enthalpy [J/kg]
					phs[facei]  = thermodynamicsMap_->hMolar_Mixture_From_MoleFractions(mole_fractions.data())/pmw[facei];	// total enthalpy [J/kg]
					phs[facei] -= mass_fractions.dot(h0_species_);								// sensible enthalpy [J/kg]
				}

				// Calculate the total concentration [kmol/m3]
				pcTot[facei] = pp[facei]/(PhysicalConstants::R_J_kmol*pT[facei]);

				// Calculate the compressibility [s2/m2]
				ppsi[facei]  = pcTot[facei]*pmw[facei]/pp[facei];

				// Calculate the constant pressure specific heat [J/kg/K]
		    		{
					pcp[facei]   = thermodynamicsMap_->cpMolar_Mixture_From_MoleFractions(mole_fractions.data());	//[J/Kmol/K]
					pcpv[facei]  = pcp[facei] - PhysicalConstants::R_J_kmol;					//[J/Kmol/K]
					pcp[facei]  /= pmw[facei];
					pcpv[facei] /= pmw[facei];
		
					{
						thermodynamicsMap_->cpMolar_Species(cp_species.data());
						for(unsigned int i=0;i<ns_;++i)
			                                CpSpecies_[i].boundaryFieldRef()[patchi][facei] = cp_species(i) / thermodynamicsMap_->MW(i);
					}

					{
						thermodynamicsMap_->hMolar_Species(h_species.data());									// total enthalpies of species [J/kg]
						for(unsigned int i=0;i<ns_;++i)
			                                hsSpecies_[i].boundaryFieldRef()[patchi][facei] = h_species(i) / thermodynamicsMap_->MW(i) -h0_species_(i);	// sensible enthalpies of species [J/kg]
					}
				} 
							
				// Calculate the viscosity [kg/m/s]
				{
					pmu[facei] = transportMap_->DynamicViscosity(mole_fractions.data());
				}
					
				// Calculate the thermal conductivity [W/m/K]
		    		{
					plambda[facei] = transportMap_->ThermalConductivity(mole_fractions.data());
				} 

				// Thermal diffusivity [m2/s]
				{
					palpha[facei] = plambda[facei]/pcp[facei]/(pp[facei]*ppsi[facei]);
				}

				// Mass diffusion coefficients [m2/s]
				if (diffusivityModel_ == DIFFUSIVITY_MODEL_MIX_AVERAGED)
				{
					transportMap_->MassDiffusionCoefficients( Dmix_species.data(),mole_fractions.data(), transportMap_->is_species_bundling() );
					for(int i=0;i<ns_;++i)
						Dmix_[i].boundaryFieldRef()[patchi][facei] = Dmix_species(i);
				}
				else if (diffusivityModel_ == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
				{
					for(int i=0;i<ns_;++i)
						Dmix_[i].boundaryFieldRef()[patchi][facei] = palpha[facei]/lewis_numbers_(i);
				}

				// Thermal diffusion coefficients [-]
				if (soret_effect_ == true)
				{
					transportMap_->ThermalDiffusionRatios(tetamix_vector.data(), mole_fractions.data());		
					for(int i=0;i<transportMap_->iThermalDiffusionRatios().size();i++)
					{
						unsigned int index = transportMap_->iThermalDiffusionRatios()[i]-1;
						Dsoret_[i].boundaryFieldRef()[patchi][facei] = Dmix_[index].boundaryFieldRef()[patchi][facei]*tetamix_vector(index)*thermodynamicsMap_->MW(index)/pmw[facei];
					} 
				}

			}
		}
	}

        // Update properties for Polimi Soot Model
	if (polimiSoot_ == true)
	{
		// Correct soot diffusion coefficients
		if ( (polimiSootModel_().physicalDiffusion() == true) && (diffusivityModel_ == DIFFUSIVITY_MODEL_MIX_AVERAGED) )
		{
			Info << "[Physical diffusion model for BINs...] ";
			const volScalarField& DmixReference = Dmix_[polimiSootModel_().physicalDiffusionReferenceIndex()];

			for(int i=0;i<polimiSootModel_().physicalDiffusionCorrectionIndex().size();i++)
			{
				const double teta = polimiSootModel_().physicalDiffusionCorrection()[i];
				if (teta > 0.)
				{
					const unsigned int index = polimiSootModel_().physicalDiffusionCorrectionIndex()[i];
					Dmix_[index] = DmixReference * teta;
				}
			}
		}

		// Update properties
		polimiSootModel_().UpdateProperties(p_, T_, Y_);
	}

	// Update properties for HMOM Model
	if (hmom_ == true)
	{
		hmomModel_().UpdateProperties(p_, T_);
	}

	// Update properties for DRG Model
	if (drg_ == true)
	{
		// nothing to do
	}

	// Update properties for Clustering Model
	if (clustering_ == true)
	{
		// nothing to do
	}

	// Update properties for Repairing Model
	if (repairing_ == true)
	{
		// nothing to do
	}

	// Update properties for TiO2 Model
	if (tio2_ == true)
	{
		const volScalarField rho = this->cTot_*this->MW_;
		tio2Model_().UpdateProperties(p_, T_, rho, mu_);
	}

	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(ncells_)*1000. << " ms per cell)" << endl;
}


void Foam::OpenSMOKEppReactingMixture::update_properties()
{
	// Force the saving of the old-time values
	this->psi_.oldTime();

	// Update the properties
	calculate();
}

#if STEADYSTATE == 1
void Foam::OpenSMOKEppReactingMixture::update_properties_for_steady_state()
{
	// Force the saving of the old-time values
	this->psi_.oldTime();

	if (energyEquation_ == ENERGY_EQ_TEMPERATURE)
	{
		// Update the properties
		if (chemistryLinearModel_().propertiesCounter() == chemistryLinearModel_().propertiesUpdate())
			calculate();
	
		// Update the counter
		chemistryLinearModel_().updatePropertiesCounter();
	}
	else
	{
		calculate();
	}
}
#endif

void Foam::OpenSMOKEppReactingMixture::update_transport_terms(const fvMesh& mesh, const Foam::volScalarField& rho)
{

	Info << "Transport terms evaluation... " ;

	const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	// Evaluation of diffusion fluxes
	if (diffusion_via_mole_fraction_gradient_ == true)
	{
		forAll(Y_, i)
		{
			volScalarField& Xi = X_[i];
			volScalarField& Yi = Y_[i];
			volScalarField& Dmixi = Dmix_[i];
			volVectorField& Ji = J_[i];
			dimensionedScalar MWi = mw_species_[i];

			Ji = -rho*Dmixi*MWi/MW_*fvc::grad(Xi);
		}
	}
	else
	{
		forAll(Y_, i)
		{
			const volScalarField& Yi = Y_[i];
			volScalarField& Dmixi = Dmix_[i];
			volVectorField& Ji = J_[i];
			Ji = -rho*Dmixi*fvc::grad(Yi);
		}
	}


	// Soret effect
	if (soret_effect_ == true)
	{
		forAll (Y_,i)
		{
			if (soret_list_[i] == true)
			{
				volVectorField& Ji = J_[i];
				Ji += -rho*Dsoret_[soret_index_[i]]/T_ *fvc::grad(T_);
			}
		}
	}

	// Thermophoretic effect
	if (thermophoretic_effect_ == true)
	{
		if (polimiSoot_ == true)
		{
			for (int k=0;k<polimiSootModel_().thermophoretic_effect_list().size();k++)
			{
				const int i = polimiSootModel_().thermophoretic_effect_list()[k];

				const volScalarField& Yi = Y_[i];
				volVectorField& Ji = J_[i];

				Ji += -0.538*mu_/T_*Yi * fvc::grad(T_) ;
			}
		}
	}

	// Correction of diffusion fluxes
	if (diffusion_fluxes_correction_ == true)
	{
		// Correction flux
		Jc_ *= 0.;
		forAll (Y_,i)
		{
			volVectorField& Ji = J_[i];
			Jc_ -= Ji;
		}

		// Evaluation of correction flow rate
		phic_ = linearInterpolate(Jc_) & mesh.Sf();

		// Correction of fluxes
		forAll(Y_, i)
		{
			volScalarField& Yi = Y_[i];
			volVectorField& Ji = J_[i];

			Ji += Yi*Jc_;
		}
	}

	// Molecular weight correction
	if (diffusion_via_mole_fraction_gradient_ == true)
	{
		forAll(Y_, i)
		{
			volScalarField& Yi = Y_[i];
			volScalarField& Dmixi = Dmix_[i];
			volScalarField& sum = mwDiffusionCorrection_[i];

			volScalarField rho_Dmixi_MW_Yi  = rho*Dmixi*MW_*Yi;

			sum *= 0.;
			forAll(Y_, k)
			{
				volScalarField& Yk = Y_[k];
				dimensionedScalar MWk = mw_species_[k];

				if (i!=k) sum += fvc::laplacian(rho_Dmixi_MW_Yi/MWk, Yk);
			}
		}
	}

	// Contribution to energy equation from mass diffusion fluxes	
	if (diffusion_in_energy_equation_ == true)
	{
		if (energyEquation_ == ENERGY_EQ_TEMPERATURE)
		{
			Cps_Js_ *= 0.;
			forAll(Y_, i)
			{
				volVectorField& Ji = J_[i];
				Cps_Js_ += CpSpecies_[i]*Ji;
			}
		}

		if (energyEquation_ == ENERGY_EQ_ENTHALPY)
		{
			hs_Js_ *= 0.;
			forAll(Y_, i)
			{
				volVectorField& Ji = J_[i];
				hs_Js_ += hsSpecies_[i]*Ji;
			}

			if (conduction_via_temperature_gradient_ == true)
			{
				HsConduction_ *= 0.;
				forAll(Y_, i)
				{
					volScalarField& Yi = Y_[i];
					HsConduction_ += lambda_/Cp_ * hsSpecies_[i]*fvc::grad(Yi);
				}
			}
		}
	}

	const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;

}

Foam::wordList Foam::OpenSMOKEppReactingMixture::hBoundaryBaseTypes()
{
	const volScalarField::Boundary& tbf = this->T_.boundaryField();

	wordList hbt(tbf.size(), word::null);

	forAll(tbf, patchi)
	{
		if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
		{
			Info << "Enthalpy Patch BaseType " << patchi << ": fixedJumpFvPatchScalarField" << endl;
			const fixedJumpFvPatchScalarField& pf = dynamic_cast<const fixedJumpFvPatchScalarField&>(tbf[patchi]);
			hbt[patchi] = pf.interfaceFieldType();
		}
		else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
		{
			Info << "Enthalpy Patch BaseType " << patchi << ": fixedJumpAMIFvPatchScalarField" << endl;
			const fixedJumpAMIFvPatchScalarField& pf = dynamic_cast<const fixedJumpAMIFvPatchScalarField&>( tbf[patchi] );
			hbt[patchi] = pf.interfaceFieldType();
		}
	}

	return hbt;
}
 
 
Foam::wordList Foam::OpenSMOKEppReactingMixture::hBoundaryTypes()
{
	const volScalarField::Boundary& tbf = this->T_.boundaryField();

	wordList hbt = tbf.types();

	forAll(tbf, patchi)
	{
		if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
		{
			Info << "Enthalpy Patch " << patchi << ": " << fixedEnthalpyFvPatchScalarField::typeName << endl;
			hbt[patchi] = fixedEnthalpyFvPatchScalarField::typeName;
		}
		else if ( isA<zeroGradientFvPatchScalarField>(tbf[patchi]) || isA<fixedGradientFvPatchScalarField>(tbf[patchi]) )
		{
			Info << "Enthalpy Patch " << patchi << ": " << gradientEnthalpyFvPatchScalarField::typeName << endl;
			hbt[patchi] = gradientEnthalpyFvPatchScalarField::typeName;
		}
		else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
		{
			Info << "Enthalpy Patch " << patchi << ": " << mixedEnthalpyFvPatchScalarField::typeName << endl;
			hbt[patchi] = mixedEnthalpyFvPatchScalarField::typeName;
		}
		else if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
		{
			Info << "Enthalpy Patch " << patchi << ": " << enthalpyJumpFvPatchScalarField::typeName << endl;
			hbt[patchi] = enthalpyJumpFvPatchScalarField::typeName;
		}
		else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
		{
			Info << "Enthalpy Patch " << patchi << ": " << enthalpyJumpAMIFvPatchScalarField::typeName << endl;
			hbt[patchi] = enthalpyJumpAMIFvPatchScalarField::typeName;
		}
		else if (tbf[patchi].type() == "energyRegionCoupledFvPatchScalarField")
		{
			Info << "Enthalpy Patch " << patchi << ": " << "energyRegionCoupledFvPatchScalarField" << endl;
			hbt[patchi] = "energyRegionCoupledFvPatchScalarField";
		}
	}

	return hbt;
}

void Foam::OpenSMOKEppReactingMixture::hBoundaryCorrection(Foam::volScalarField& h)
{
	volScalarField::Boundary& hBf = h.boundaryFieldRef();

	forAll(hBf, patchi)
	{
		if (isA<gradientEnthalpyFvPatchScalarField>(hBf[patchi]))
		{
			refCast<gradientEnthalpyFvPatchScalarField>(hBf[patchi]).gradient() = hBf[patchi].fvPatchField::snGrad();
		}
		else if (isA<mixedEnthalpyFvPatchScalarField>(hBf[patchi]))
		{
			refCast<mixedEnthalpyFvPatchScalarField>(hBf[patchi]).refGrad() = hBf[patchi].fvPatchField::snGrad();
		}
	}
}

void Foam::OpenSMOKEppReactingMixture::init()
{
	Info << " * Initialize the enthalpy field..." << endl;

	// Auxiliary variables
	Eigen::VectorXd mass_fractions(ns_);
	Eigen::VectorXd mole_fractions(ns_);

	// Const fields
	const scalarField& pCells = this->p_;
	const scalarField& TCells = this->T_;

	// Non-const fields
	scalarField& hsCells = this->Hs_.primitiveFieldRef();

	// Repairing solution
	if (repairing_ == true)
	{
		repairingModel_().RepairMassFractions(this->T_.mesh(), Y_);
		repairingModel_().RepairMinMaxFields(this->T_.mesh());
		repairingModel_().RepairPatch(this->T_.mesh(), Y_);
	//	repairingModel_().RepairVelocity(this->T_.mesh());
	}

	// Loop over all the internal cells
	Info << "   Internal cells... " << endl;
	Info << "    - Min/Max T: " << min(this->T_).value() << "/" << max(this->T_).value() << endl;
	Info << "    - Min/Max p: " << min(this->p_).value() << "/" << max(this->p_).value() << endl;
	for(int i=0;i<ns_;++i)
		Info << "    - Min/Max " << thermodynamicsMap_->NamesOfSpecies()[i] << ": " << min(Y_[i]).value() << "/" << max(Y_[i]).value() << endl;

	forAll(TCells, celli)
	{			
		// Calculate the mole fractions from the mass fractions
		double mw = 0.;
		for(int i=0;i<ns_;++i)
			mass_fractions(i) = std::max(0., Y_[i].internalField()[celli]);
		
		thermodynamicsMap_->MoleFractions_From_MassFractions(mole_fractions.data(), mw, mass_fractions.data());

		// Set temperature and pressure for thermodynamic and transport maps
		thermodynamicsMap_->SetPressure(pCells[celli]);
		thermodynamicsMap_->SetTemperature(TCells[celli]);

		// Calculate the enthalpy [J/kg]
		hsCells[celli]  = thermodynamicsMap_->hMolar_Mixture_From_MoleFractions(mole_fractions.data())/mw;	// total enthalpy [J/kg]
		hsCells[celli] -= mass_fractions.dot(h0_species_);							// sensible enthalpy [J/kg]
	}

	// Boundaries
	{
		volScalarField::Boundary& pBf 		= this->p_.boundaryFieldRef();
		volScalarField::Boundary& TBf 		= this->T_.boundaryFieldRef();
		volScalarField::Boundary& hsBf 		= this->Hs_.boundaryFieldRef();

		forAll(this->T_.boundaryField(), patchi)
		{
			Info << "   Patch field: " << T_.mesh().boundary()[patchi].name() << endl;

			fvPatchScalarField& pp = pBf[patchi];
			fvPatchScalarField& pT = TBf[patchi];
			fvPatchScalarField& phs = hsBf[patchi];
			
			forAll(pT, facei)
			{
				// Calculate the mole fractions from the mass fractions
				double mw = 0.;
				for(int i=0;i<ns_;++i)
					mass_fractions(i) = Y_[i].boundaryField()[patchi][facei];
				thermodynamicsMap_->MoleFractions_From_MassFractions(mole_fractions.data(),mw,mass_fractions.data());
				
				// Set temperature and pressure for thermodynamic and transport maps
				thermodynamicsMap_->SetPressure(pp[facei]);
				thermodynamicsMap_->SetTemperature(pT[facei]);

				// Calculate the enthalpy [J/kg]
				phs[facei]  = thermodynamicsMap_->hMolar_Mixture_From_MoleFractions(mole_fractions.data())/mw;	// total enthalpy [J/kg]
				phs[facei] -= mass_fractions.dot(h0_species_);							// sensible enthalpy [J/kg]
			}
		}
	}

	this->hBoundaryCorrection(Hs_);
}


Foam::tmp<Foam::scalarField> Foam::OpenSMOKEppReactingMixture::hs( const Foam::scalarField& p, const Foam::scalarField& T, const label patchi ) const
{
     Foam::tmp<Foam::scalarField> ths(new Foam::scalarField(T.size()));
     Foam::scalarField& phs = ths.ref();

     Eigen::VectorXd mass_fractions(ns_);
     Eigen::VectorXd mole_fractions(ns_);
 
     forAll(T, facei)
     {
		// Mole fractions and molecular weight from mass fractions
		double mw = 0.;
		for(int i=0;i<ns_;++i)
			mass_fractions(i) = Y_[i].boundaryField()[patchi][facei];
		thermodynamicsMap_->MoleFractions_From_MassFractions(mole_fractions.data(),mw,mass_fractions.data());

		// Set temperature and pressure for thermodynamic and transport maps
		thermodynamicsMap_->SetPressure(p[facei]);
		thermodynamicsMap_->SetTemperature(T[facei]);

		// Calculate the enthalpy [J/kg]
		phs[facei]  = thermodynamicsMap_->hMolar_Mixture_From_MoleFractions(mole_fractions.data())/mw;	// total enthalpy [J/kg]
		phs[facei] -= mass_fractions.dot(h0_species_);							// sensible enthalpy [J/kg]
     }
 
     return ths;
}

Foam::tmp<Foam::scalarField> Foam::OpenSMOKEppReactingMixture::hs( const Foam::scalarField& p, const Foam::scalarField& T, const labelList& cells ) const
{
     Foam::tmp<Foam::scalarField> ths(new Foam::scalarField(T.size()));
     Foam::scalarField& hsCells = ths.ref();
     
     Eigen::VectorXd mass_fractions(ns_);
     Eigen::VectorXd mole_fractions(ns_);

     forAll(T, celli)
     {
		// Mole fractions and molecular weight from mass fractions
		double mw = 0.;
		for(int i=0;i<ns_;++i)
			mass_fractions(i) = Y_[i].internalField()[celli];
		thermodynamicsMap_->MoleFractions_From_MassFractions(mole_fractions.data(), mw, mass_fractions.data());

		// Set temperature and pressure for thermodynamic and transport maps
		thermodynamicsMap_->SetPressure(p[celli]);
		thermodynamicsMap_->SetTemperature(T[celli]);

		// Calculate the enthalpy [J/kg]
		hsCells[celli]  = thermodynamicsMap_->hMolar_Mixture_From_MoleFractions(mole_fractions.data())/mw;	// total enthalpy [J/kg]
		hsCells[celli] -= mass_fractions.dot(h0_species_);							// sensible enthalpy [J/kg]
     }
 
     return ths;
}

#if STEADYSTATE == 1
void Foam::OpenSMOKEppReactingMixture::update_chemical_source_terms()
{
	if (homogeneousReactions_ == true && (solveForSpeciesEquations_ == true || solveForEnergyEquation_ == true))
	{
		if (chemistryLinearModel_().implicitSourceTerm() == true && chemistryLinearModel_().jacobianUpdate() == 1)
		{
			Info<< " * Implicit chemical source terms evaluation (continuous)... " << endl;
		
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
			{
				const scalarField& TCells = T_.internalField();
				const scalarField& pCells = p_.internalField(); 
				scalarField& QCells = Q_.ref();
	
				Eigen::VectorXd J(ns_+1);
				Eigen::VectorXd source(ns_+1);	
				Eigen::VectorXd y(ns_+1);

				forAll(TCells, celli)
				{
					for(int i=0;i<ns_;i++)
						y(i) = Y_[i].internalField()[celli];
					y(ns_) = TCells[celli];

					chemistryLinearModel_().reactionSourceTerms( y, pCells[celli], source);

					if (chemistryLinearModel_().sparseJacobian() == false)
						chemistryLinearModel_().reactionJacobian( y, pCells[celli], J );
					else
						chemistryLinearModel_().reactionJacobianSparse( y, pCells[celli], J, solveForEnergyEquation_, solveForSpeciesEquations_ );

					for(int i=0;i<ns_+1;i++)
						sourceImplicit_[i].ref()[celli] = J(i);
			
					for(int i=0;i<ns_+1;i++)
						sourceExplicit_[i].ref()[celli] = source(i) - J(i)*y(i);

					// Heat release rate (W/m3)
					QCells[celli] = source(ns_);
				}
			}
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(ncells_)*1000. << " ms per cell)" << endl;
		}

		else if (chemistryLinearModel_().implicitSourceTerm() == true && chemistryLinearModel_().jacobianUpdate() != 1)
		{
			Info<< " * Implicit chemical source terms evaluation (discrete)... " << endl;

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
			{
				const scalarField& TCells = T_.internalField();
				const scalarField& pCells = p_.internalField(); 
				scalarField& QCells = Q_.ref();
	
				Eigen::VectorXd J(ns_+1);
				Eigen::VectorXd source(ns_+1);	
				Eigen::VectorXd y(ns_+1);

				forAll(TCells, celli)
				{
					for(int i=0;i<ns_;i++)
						y(i) = Y_[i].internalField()[celli];
					y(ns_) = TCells[celli];

					chemistryLinearModel_().reactionSourceTerms( y, pCells[celli], source );

					if (chemistryLinearModel_().jacobianCounter() == chemistryLinearModel_().jacobianUpdate())
					{
						if (chemistryLinearModel_().sparseJacobian() == false)
							chemistryLinearModel_().reactionJacobian( y, pCells[celli], J );
						else
							chemistryLinearModel_().reactionJacobianSparse( y, pCells[celli], J, solveForEnergyEquation_, solveForSpeciesEquations_);

						Jstored_[celli] = J;
					}

					for(int i=0;i<ns_+1;i++)
						sourceImplicit_[i].ref()[celli] = Jstored_[celli](i);
		
					for(int i=0;i<ns_+1;i++)
						sourceExplicit_[i].ref()[celli] = source(i) - Jstored_[celli](i)*y(i);

					// Heat release rate (W/m3)
					QCells[celli] = source(ns_);
				}
			}
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(ncells_)*1000. << " ms per cell)" << endl;

			// Update the counter
			chemistryLinearModel_().updateJacobianCounter();
		}

		else
		{
			Info<< " * Explicit source terms evaluation... " << endl;

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
			{
				const scalarField& TCells = T_.internalField();
				const scalarField& pCells = p_.internalField(); 
	
				Eigen::VectorXd source(ns_+1);	
				Eigen::VectorXd y(ns_+1);

				forAll(TCells, celli)
				{
					for(int i=0;i<ns_;i++)
						y(i) = Y_[i].internalField()[celli];
					y(ns_) = TCells[celli];

					chemistryLinearModel_().reactionSourceTerms( y, pCells[celli], source );

					for(int i=0;i<ns_+1;i++)
						sourceImplicit_[i].ref()[celli] = 0.;
		
					for(int i=0;i<ns_+1;i++)
						sourceExplicit_[i].ref()[celli] = source(i);
				}
			}
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(ncells_)*1000. << " ms per cell)" << endl;
		}
	}
}

void Foam::OpenSMOKEppReactingMixture::update_species_order_policy() 
{ 
	chemistryLinearModel_().updateSpeciesOrderPolicy(); 
}
		
label Foam::OpenSMOKEppReactingMixture::species_order(const label i) const 
{ 
	return chemistryLinearModel_().species_order(i); 
}

#endif


Foam::tmp<Foam::scalarField> Foam::OpenSMOKEppReactingMixture::Cpv( const Foam::scalarField& p, const Foam::scalarField& T, const label patchi ) const
{
     Foam::tmp<Foam::scalarField> tCpv(new Foam::scalarField(T.size()));
     Foam::scalarField& pcp = tCpv.ref();

     Eigen::VectorXd mass_fractions(ns_);
     Eigen::VectorXd mole_fractions(ns_);
 
     forAll(T, facei)
     {
		// Mole fractions and molecular weight from mass fractions
		double mw = 0.;
		for(int i=0;i<ns_;++i)
			mass_fractions(i) = Y_[i].boundaryField()[patchi][facei];
		thermodynamicsMap_->MoleFractions_From_MassFractions(mole_fractions.data(),mw,mass_fractions.data());

		// Set temperature and pressure for thermodynamic and transport maps
		thermodynamicsMap_->SetPressure(p[facei]);
		thermodynamicsMap_->SetTemperature(T[facei]);

		pcp[facei] = thermodynamicsMap_->cpMolar_Mixture_From_MoleFractions(mole_fractions.data());	// [J/kmol/K]
		pcp[facei] = pcp[facei]/mw;									// [J/kg/K]
     }
 
     return tCpv;
}

double Foam::OpenSMOKEppReactingMixture::normalizeMassFractions(Eigen::VectorXd& omega_plus_temperature, const label celli)
{
	double sumMassFractions = 0.;
	for(int i=0; i < ns_; i++)
		sumMassFractions += omega_plus_temperature(i);

	for(int i=0; i < ns_; i++)
		omega_plus_temperature(i) /= sumMassFractions;

	return std::fabs(sumMassFractions-1.);
}

void Foam::OpenSMOKEppReactingMixture::onthefly_postprocessing(const fvMesh& mesh, const Foam::volScalarField& rho)
{
	if (polimiSoot_ == true)
	{
		polimiSootModel_().OnTheFlyAnalysis(mesh, rho, Y_);
		polimiSootModel_().FieldsReconstruction(mesh, p_, T_, Y_);
	}

	if (hmom_ == true)
	{
		hmomModel_().OnTheFlyAnalysis(mesh, rho, Y_);
		hmomModel_().FieldsReconstruction(mesh, p_, T_);
	}

	if (drg_ == true)
	{
		// nothing to do
	}

	if (clustering_ == true)
	{
		// nothing to do
	}

	if (repairing_ == true)
	{
		// nothing to do
	}

	if (tio2_ == true)
	{
		tio2Model_().OnTheFlyAnalysis(mesh, rho, Y_);
		tio2Model_().FieldsReconstruction(mesh, rho, T_);
	}
}

void Foam::OpenSMOKEppReactingMixture::solve_additional_equations(const fvMesh& mesh, const Foam::volScalarField& rho, const surfaceScalarField& phi)
{
	if (polimiSoot_ == true)
	{
		// nothing to do
	}

	if (hmom_ == true)
	{
		hmomModel_().SolveEquations(mesh, p_, T_, phi, rho, mu_, thermophoretic_effect_);
	}

	if (drg_ == true)
	{
		// nothing to do
	}

	if (clustering_ == true)
	{
		// nothing to do
	}

	if (repairing_ == true)
	{
		// nothing to do
	}

	if (tio2_ == true)
	{
		tio2Model_().SolveEquations(mesh, p_, T_, phi, rho, mu_);
	}
}

#include "chemistry_DI.H"


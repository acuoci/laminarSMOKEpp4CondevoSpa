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
const Foam::word Foam::NOxPostProcessor::dictName("NOxPostProcessor");

const double Foam::NOxPostProcessor::MW_NO_ = 30.006000;	// molecular weight of NO (kg/kmol)
const double Foam::NOxPostProcessor::MW_N2_ = 28.014000;	// molecular weight of N2 (kg/kmol)
const double Foam::NOxPostProcessor::MW_O2_ = 31.997999;	// molecular weight of O2 (kg/kmol)
const double Foam::NOxPostProcessor::MW_H2O_ = 18.015000;	// molecular weight of H2O (kg/kmol)
const double Foam::NOxPostProcessor::MW_O_ = 15.999000;		// molecular weight of O (kg/kmol)
const double Foam::NOxPostProcessor::MW_OH_ = 17.007000;	// molecular weight of OH (kg/kmol)

// Constructor
Foam::NOxPostProcessor::NOxPostProcessor(const Foam::fvMesh& mesh, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap)
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

	YNO_
	(
	    IOobject
	    (
		"nox:NO",
		mesh.time().timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	),

	gamma_
	(
	    IOobject
	    (
		"nox:gamma",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("nox:gamma", dimensionSet(1., -1., -1., 0., 0.), 0.)
	),

	SNO_thermal_
	(
	    IOobject
	    (
		"nox:SNO_thermal",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("nox:S", dimensionSet(1., -3., -1., 0., 0.), 0.)
	),

	SNO_prompt_
	(
	    IOobject
	    (
		"nox:SNO_prompt",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("nox:S", dimensionSet(1., -3., -1., 0., 0.), 0.)
	),

	SNO_intermediate_n2o_
	(
	    IOobject
	    (
		"nox:SNO_int_n2o",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("nox:S", dimensionSet(1., -3., -1., 0., 0.), 0.)
	),

	sourceExplicitNO_
	(
	    IOobject
	    (
		"nox:SE_NO",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("nox:S", dimensionSet(1., -3., -1., 0., 0.), 0.)
	),

	sourceImplicitNO_
	(
	    IOobject
	    (
		"nox:SI_NO",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("nox:S", dimensionSet(1., -3., -1., 0., 0.), 0.)
	)
{
	// Total number of equations
	ne_ = 1;

	// Kinetic mechanism
	NS_ = thermodynamicsMap_.NumberOfSpecies();
	NR_ = kineticsMap_.NumberOfReactions();

	// User-defined options
	{
		word name_FUEL(this->lookup("FUEL"));
		word name_N2(this->lookup("N2"));
		word name_O2(this->lookup("O2"));
		word name_H2O(this->lookup("H2O"));
		word name_O(this->lookup("O"));
		word name_OH(this->lookup("OH"));
		word name_CH(this->lookup("CH"));

		index_fuel_ = thermodynamicsMap_.IndexOfSpecies(name_FUEL)-1;
		index_n2_ = thermodynamicsMap_.IndexOfSpecies(name_N2)-1;	
		index_o2_ = thermodynamicsMap_.IndexOfSpecies(name_O2)-1;
		index_h2o_ = thermodynamicsMap_.IndexOfSpecies(name_H2O)-1;
		index_o_ = thermodynamicsMap_.IndexOfSpeciesWithoutError(name_O)-1;	
		index_oh_ = thermodynamicsMap_.IndexOfSpeciesWithoutError(name_OH)-1;
		index_ch_ = thermodynamicsMap_.IndexOfSpeciesWithoutError(name_CH)-1;	

		MW_Fuel_ = thermodynamicsMap_.MWs()[index_fuel_];
	}

	// Activation of processes
	{
		is_thermal_ = Switch(this->lookup(word("thermalNOx")));
		is_prompt_ = Switch(this->lookup(word("promptNOx")));
		is_fuel_nox_ = Switch(this->lookup(word("fuelNOx")));
		is_intermediate_n2o_ = Switch(this->lookup(word("intermediateN2O")));
		is_reduction_by_reburning_ = Switch(this->lookup(word("reductionReburning")));
		is_reduction_by_sncr_ = Switch(this->lookup(word("reductionSNCR")));
	}

	// Numerical parameters
	{
		is_implicit_mode_ = Switch(this->lookup(word("implicitMode")));
	}

	// On-the-fly analysis
	{
		ontheflyAnalysis_ = false;
		frequency_ = 10;
		frequencyCounter_ = 10;	
		ontheflyAnalysis_ = Switch(this->lookup(word("onTheFlyAnalysis")));
		frequency_  = readLabel(this->lookup("frequency"));
		frequencyCounter_ = frequency_;
	}


	// Setup of NOx post-processor model
	{
		model_NOxPostProcessor_ = new ModelNOxPostProcessor();

		// Activation of subprocesses
		model_NOxPostProcessor_->SetThermalNOx(is_thermal_);
		model_NOxPostProcessor_->SetPromptNOx(is_prompt_);
		model_NOxPostProcessor_->SetIntermediateN2O(is_intermediate_n2o_);

		// O radical mode
		{	
			Foam::string option = this->lookup("radicalModelO");
			if (option == "equilibrium")			model_NOxPostProcessor_->SetORadicalMode(ModelNOxPostProcessor::O_RADICAL_EQUILIBRIUM);
			else if (option == "partial-equilibrium")	model_NOxPostProcessor_->SetORadicalMode(ModelNOxPostProcessor::O_RADICAL_PARTIAL_EQUILIBRIUM);
			else if (option == "kinetics")			model_NOxPostProcessor_->SetORadicalMode(ModelNOxPostProcessor::O_RADICAL_KINETICS);
			else FatalError << "Wrong radicalModelO option. Allowed options: equilibrium | partial-equilibrium | kinetics" << ::Foam::exit(FatalError); 
		}

		// OH radical mode
		{	
			Foam::string option = this->lookup("radicalModelOH");
			if (option == "exclusion")			model_NOxPostProcessor_->SetOHRadicalMode(ModelNOxPostProcessor::OH_RADICAL_EXCLUSION);
			else if (option == "partial-equilibrium")	model_NOxPostProcessor_->SetOHRadicalMode(ModelNOxPostProcessor::OH_RADICAL_PARTIAL_EQUILIBRIUM);
			else if (option == "kinetics")			model_NOxPostProcessor_->SetOHRadicalMode(ModelNOxPostProcessor::OH_RADICAL_KINETICS);
			else FatalError << "Wrong radicalModelOH option. Allowed options: exclusion | partial-equilibrium | kinetics" << ::Foam::exit(FatalError); 
		}
	}


	// NOx integrals
	if (ontheflyAnalysis_ == true)
	{
	}

	
	// Summary on the screen
	{
		std::cout << "------------------------------------------------------------------------------------------" << std::endl; 
		std::cout << "                         NOxPostProcessor Solver Summary                                  " << std::endl;
		std::cout << "------------------------------------------------------------------------------------------" << std::endl; 

		Info << "FUEL species:    " << thermodynamicsMap_.NamesOfSpecies()[index_fuel_] << " (" << index_fuel_ << ")" << endl;
		Info << "N2 species:      " << thermodynamicsMap_.NamesOfSpecies()[index_n2_] << " (" << index_n2_ << ")" << endl;
		Info << "O2 species:      " << thermodynamicsMap_.NamesOfSpecies()[index_o2_] << " (" << index_o2_ << ")" << endl;
		Info << "H2O species:     " << thermodynamicsMap_.NamesOfSpecies()[index_h2o_] << " (" << index_h2o_ << ")" << endl;

		if (index_o_ >= 0)
			Info << "O species:       " << thermodynamicsMap_.NamesOfSpecies()[index_o_] << " (" << index_o_ << ")" << endl;
		if (index_oh_ >= 0)
			Info << "OH species:      " << thermodynamicsMap_.NamesOfSpecies()[index_oh_] << " (" << index_oh_ << ")" << endl;
		if (index_ch_ >= 0)
			Info << "CH species:      " << thermodynamicsMap_.NamesOfSpecies()[index_ch_] << " (" << index_ch_ << ")" << endl;
	
		model_NOxPostProcessor_->SummaryOnScreen();

		std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	}
}

Foam::NOxPostProcessor::~NOxPostProcessor()
{}

const Foam::volScalarField& Foam::NOxPostProcessor::YNO() const
{
	return YNO_;
}

Foam::volScalarField& Foam::NOxPostProcessor::YNO()
{
	return YNO_;
}

void Foam::NOxPostProcessor::Setup(const Foam::fvMesh& mesh)
{
	
}

void Foam::NOxPostProcessor::reactionSourceTerms( const Eigen::VectorXd& Y, const double T0, const double P0, const double rho0, const Eigen::VectorXd& Y0, Eigen::VectorXd& S )
{
	const double cTot = P0/PhysicalConstants::R_J_kmol/T0;	// (kmol/m3)
	const double MWmix = rho0/cTot;				// (kg/kmol)

	// Evaluation of concentration of relevant species (kmol/m3)
	const double cNO = cTot * (Y(0) * MWmix / MW_NO_);

	// Evaluation of concentration of relevant species (kmol/m3)
	const double cN2 = cTot * (Y0(0) * MWmix / MW_N2_);
	const double cO2 = cTot * (Y0(1) * MWmix / MW_O2_);
	const double cH2O = cTot * (Y0(2) * MWmix / MW_H2O_);
	const double cFuel = cTot * (Y0(3) * MWmix / MW_Fuel_);

	double cO = 0.;
	if (index_o_ >= 0)
		cO = cTot * (Y0(4) * MWmix / MW_O_);

	double cOH = 0.;
	if (index_oh_ >= 0)
		cOH = cTot * (Y0(5) * MWmix / MW_OH_);

			
	// Source terms
	double OmegaThermal, OmegaPrompt, OmegaIntermediateN2O;
	model_NOxPostProcessor_->ThermalNOx( T0, cN2, cO2, cH2O, cO, cOH, cNO, OmegaThermal);
	model_NOxPostProcessor_->PromptNOx( T0, P0, cN2, cO2, cFuel, OmegaPrompt);
	model_NOxPostProcessor_->IntermediateN2O( T0, P0, cN2, cO2, cO, cNO, OmegaIntermediateN2O);

	S(0) = OmegaThermal + OmegaPrompt + OmegaIntermediateN2O;
}

void Foam::NOxPostProcessor::reactionJacobian( const Eigen::VectorXd& y, const double T0, const double P0, const double rho0, const Eigen::VectorXd& y0, Eigen::VectorXd &J)
{
	// Calculated as suggested by Guido Buzzi-Ferraris (private communication)
	const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
	const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
	const double BETA = 1.e3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;
	const double TOLR = 1.e-7;
	const double TOLA = 1.e-12;

	// Memory allocation
	Eigen::VectorXd dy_original(ne_);
	Eigen::VectorXd y_plus(ne_);
	Eigen::VectorXd dy_plus(ne_);

	// Start calculations
	y_plus = y;

	// Call equations
	reactionSourceTerms(y, T0, P0, rho0, y0, dy_original);     

	// Derivatives with respect to y(kd)
	for(int kd=0;kd<ne_;kd++)
	{
		const double hf = 1.e0;
		const double error_weight = 1./(TOLA+TOLR*std::fabs(y(kd)));
		      double hJ = ETA2 * std::fabs(std::max(y(kd), 1./error_weight));
		const double hJf = hf/error_weight;
		hJ = std::max(hJ, hJf);
		hJ = std::max(hJ, ZERO_DER);

		// This is what is suggested by Guido Buzzi-Ferraris (private communication)
		const double dy = std::min(hJ, 1.e-3 + 1e-3*std::fabs(y(kd)));
		const double udy = 1. / dy;
		y_plus(kd) += dy;

		reactionSourceTerms(y_plus, T0, P0, rho0, y0, dy_plus);

		J(kd) = (dy_plus(kd)-dy_original(kd)) * udy;

		y_plus(kd) = y(kd);
	}
}


void Foam::NOxPostProcessor::SolveEquations(const Foam::fvMesh& mesh, const volScalarField& p, const volScalarField& T, const surfaceScalarField& phi, const volScalarField& rho, const PtrList<volScalarField>& Y)
{

#if STEADYSTATE == 1

	Eigen::VectorXd S(ne_);
	Eigen::VectorXd J(ne_);
	Eigen::VectorXd y(ne_);
	Eigen::VectorXd y0(6);
	y0.setZero();

	const scalarField& TCells = T.internalField();
	const scalarField& rhoCells = rho.internalField();
	const scalarField& pCells = p.internalField();
	const scalarField& YNOCells  = YNO_.internalField();

	scalarField& SExplicitCells  = sourceExplicitNO_.ref();
	scalarField& SImplicitCells  = sourceImplicitNO_.ref();

	forAll(TCells, celli)
	{
		y(0) = YNOCells[celli];

		y0(0) = Y[index_n2_].internalField()[celli];	
		y0(1) = Y[index_o2_].internalField()[celli];
		y0(2) = Y[index_h2o_].internalField()[celli];
		y0(3) = Y[index_fuel_].internalField()[celli];
		if (index_o_ >= 0)
			y0(4) = Y[index_o_].internalField()[celli];

		if (index_oh_ >= 0)
			y0(5) = Y[index_oh_].internalField()[celli];

		// Source terms
		reactionSourceTerms( y, TCells[celli], pCells[celli], rhoCells[celli], y0, S );

		if (is_implicit_mode_ == true)
		{
			reactionJacobian( y, TCells[celli], pCells[celli], rhoCells[celli], y0, J);
			SImplicitCells[celli] = J(0);
			SExplicitCells[celli] = S(0)-J(0)*y(0);
		}
		else
		{
			SExplicitCells[celli] = S(0);
		}
	}

	if (is_implicit_mode_ == true)
	{
		fvScalarMatrix YNOEqn
		(
			  fvm::div(phi, YNO_)
			- fvm::laplacian(gamma_, YNO_)
			==
			sourceExplicitNO_ + fvm::Sp(sourceImplicitNO_, YNO_)
		);

		// Addiitonal operations
		YNOEqn.relax();
		YNOEqn.solve("nox:NO");

		// Correct
		YNO_.max(0.);
		YNO_.min(1.);
	}
	else
	{	
		fvScalarMatrix YNOEqn
		(
			  fvm::div(phi, YNO_)
			- fvm::laplacian(gamma_, YNO_)
			==
			sourceExplicitNO_
		);

		// Addiitonal operations
		YNOEqn.relax();
		YNOEqn.solve("nox:NO");

		// Correct
		YNO_.max(0.);
		YNO_.min(1.);
	}

#else
		
	// Equations
	{
		fvScalarMatrix YEqn
		(
			fvm::ddt(rho,YNO_)
		      + fvm::div(phi, YNO_) 
		      - fvm::laplacian(gamma_, YNO_) 
		);

		// Solve
		YEqn.relax();
		YEqn.solve("YNO");
		
		// Correct
		YNO_.max(0.);
		YNO_.min(1.);
	}

#endif
		
}

void Foam::NOxPostProcessor::UpdateProperties(const volScalarField& p, const volScalarField& T, const volScalarField& rho, const PtrList<volScalarField>& Dmix)
{
	gamma_ = rho*Dmix[index_n2_];
}

void Foam::NOxPostProcessor::OnTheFlyAnalysis(const Foam::fvMesh& mesh, const volScalarField& rho, const PtrList<volScalarField>& Y)
{
	if (ontheflyAnalysis_ == true && frequencyCounter_ == frequency_)
	{
		frequencyCounter_ = 0;
	}

	frequencyCounter_++;
}

void Foam::NOxPostProcessor::FieldsReconstruction(const Foam::fvMesh& mesh, const volScalarField& rho, const volScalarField& p, const volScalarField& T, const PtrList<volScalarField>& Y)
{
	if (mesh.time().outputTime() == true)
	{
		const scalarField& TCells = T.internalField();
		const scalarField& rhoCells = rho.internalField();
		const scalarField& pCells = p.internalField();

		const scalarField& YNOCells  = YNO_.internalField();

		scalarField& SNOThermalCells  = SNO_thermal_.ref();
		scalarField& SNOPromptCells  = SNO_prompt_.ref();
		scalarField& SNOIntermediateN2OCells  = SNO_intermediate_n2o_.ref();

		forAll(TCells, celli)
		{
			const double cTot = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];	// (kmol/m3)
			const double MWmix = rhoCells[celli]/cTot;					// (kg/kmol)

			// Evaluation of concentration of relevant species (kmol/m3)
			const double cNO = cTot * (YNOCells[celli] * MWmix / MW_NO_);
			const double cN2 = cTot * (Y[index_n2_].internalField()[celli] * MWmix / MW_N2_);
			const double cO2 = cTot * (Y[index_o2_].internalField()[celli] * MWmix / MW_O2_);
			const double cH2O = cTot * (Y[index_h2o_].internalField()[celli] * MWmix / MW_H2O_);
			const double cFuel = cTot * (Y[index_fuel_].internalField()[celli] * MWmix / MW_Fuel_);

			double cO = 0.;
			if (index_o_ >= 0)
				cO = cTot * (Y[index_o_].internalField()[celli] * MWmix / MW_O_);

			double cOH = 0.;
			if (index_oh_ >= 0)
				cOH = cTot * (Y[index_oh_].internalField()[celli] * MWmix / MW_OH_);

			
			// Source terms
			model_NOxPostProcessor_->ThermalNOx( TCells[celli], cN2, cO2, cH2O, cO, cOH, cNO, SNOThermalCells[celli]);
			model_NOxPostProcessor_->PromptNOx( TCells[celli], pCells[celli], cN2, cO2, cFuel, SNOPromptCells[celli]);
			model_NOxPostProcessor_->IntermediateN2O( TCells[celli], pCells[celli], cN2, cO2, cO, cNO, SNOIntermediateN2OCells[celli]);

		}

		forAll(T.boundaryField(), patchi)
		{
			const fvPatchScalarField& pT  = T.boundaryField()[patchi];
			const fvPatchScalarField& prho  = rho.boundaryField()[patchi];
			const fvPatchScalarField& pp  = p.boundaryField()[patchi];

			const fvPatchScalarField& pyno  = YNO_.boundaryField()[patchi];

			fvPatchScalarField& pSNOThermal  = SNO_thermal_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pSNOPrompt  = SNO_prompt_.boundaryFieldRef()[patchi];
			fvPatchScalarField& pSNOIntermediateN2O  = SNO_intermediate_n2o_.boundaryFieldRef()[patchi];

			forAll(pT, facei)
			{
				const double cTot = pp[facei]/PhysicalConstants::R_J_kmol/pT[facei];	// (kmol/m3)
				const double MWmix = prho[facei]/cTot;					// (kg/kmol)

				// Evaluation of concentration of relevant species (kmol/m3)
				const double cNO = cTot * (pyno[facei] * MWmix / MW_NO_);
				const double cN2 = cTot * (Y[index_n2_].boundaryField()[patchi][facei] * MWmix / MW_N2_);
				const double cO2 = cTot * (Y[index_o2_].boundaryField()[patchi][facei] * MWmix / MW_O2_);
				const double cH2O = cTot * (Y[index_h2o_].boundaryField()[patchi][facei] * MWmix / MW_H2O_);
				const double cFuel = cTot * (Y[index_fuel_].boundaryField()[patchi][facei] * MWmix / MW_Fuel_);

				double cO = 0.;
				if (index_o_ >= 0)
					cO = cTot * (Y[index_o_].boundaryField()[patchi][facei] * MWmix / MW_O_);

				double cOH = 0.;
				if (index_oh_ >= 0)
					cOH = cTot * (Y[index_oh_].boundaryField()[patchi][facei] * MWmix / MW_OH_);

			
				// Source terms
				model_NOxPostProcessor_->ThermalNOx( pT[facei], cN2, cO2, cH2O, cO, cOH, cNO, pSNOThermal[facei]);
				model_NOxPostProcessor_->PromptNOx( pT[facei], pp[facei], cN2, cO2, cFuel, pSNOPrompt[facei]);
				model_NOxPostProcessor_->IntermediateN2O( pT[facei], pp[facei], cN2, cO2, cO, cNO, pSNOIntermediateN2O[facei]);
			}
		}
	}
}

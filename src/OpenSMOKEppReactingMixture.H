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

#ifndef OpenSMOKEppReactingMixture_H
#define OpenSMOKEppReactingMixture_H

// OpenSMOKE++ Framework
#include "OpenSMOKEpp"
#include "maps/Maps_CHEMKIN"

// OpenSMOKE++ Dictionaries
#include "dictionary/OpenSMOKE_Dictionary"

// Include standard OpenFOAM files
#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"
#include "hashedWordList.H"

// Detailed chemistry
#if CHEMKINBASED == 1
#include "chemistry/chemkin/BatchReactorHomogeneousConstantPressure.H"
#else
#include "chemistry/BatchReactorHomogeneousConstantPressure.H"
#endif
#include "chemistry/BatchReactorHomogeneousConstantPressure_ODE_Interface.H"

// Detailed chemistry + HMOM
#include "chemistry/BatchReactorHomogeneousConstantPressure_HMOM.H"
#include "chemistry/BatchReactorHomogeneousConstantPressure_HMOM_ODE_Interface.H"

// Detailed chemistry + DRG
#include "chemistry/BatchReactorHomogeneousConstantPressure_DRG.H"
#include "chemistry/BatchReactorHomogeneousConstantPressure_DRG_ODE_Interface.H"

// Detailed chemistry + TIO2
#include "chemistry/BatchReactorHomogeneousConstantPressure_TIO2.H"
#include "chemistry/BatchReactorHomogeneousConstantPressure_TIO2_ODE_Interface.H"

// Include additional sub-models
#include "soot/polimiSoot/PolimiSootModel.H"
#include "soot/hmom/HMOMModel.H"
#include "drg/DRG.H"
#include "clustering/Clustering.H"
#include "repairing/Repairing.H"
#include "materialSynthesis/TiO2.H"

// Include chemistry linear model for steady state calculations
#if STEADYSTATE == 1
#include "ChemistryLinearModel.H"
#endif

// ODE solvers
#include "math/native-ode-solvers/MultiValueSolver"
#include "math/external-ode-solvers/ODE_Parameters.h"

namespace Foam
{
	typedef OdeSMOKE::KernelDense<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_OpenSMOKE> denseOdeConstantPressure;
	typedef OdeSMOKE::MethodGear<denseOdeConstantPressure> methodGearConstantPressure;

	typedef OdeSMOKE::KernelDense<OpenSMOKE::BatchReactorHomogeneousConstantPressure_HMOM_ODE_OpenSMOKE> denseOdeConstantPressure_HMOM;
	typedef OdeSMOKE::MethodGear<denseOdeConstantPressure_HMOM> methodGearConstantPressure_HMOM;

	typedef OdeSMOKE::KernelDense<OpenSMOKE::BatchReactorHomogeneousConstantPressure_DRG_ODE_OpenSMOKE> denseOdeConstantPressure_DRG;
	typedef OdeSMOKE::MethodGear<denseOdeConstantPressure_DRG> methodGearConstantPressure_DRG;

	typedef OdeSMOKE::KernelDense<OpenSMOKE::BatchReactorHomogeneousConstantPressure_TIO2_ODE_OpenSMOKE> denseOdeConstantPressure_TIO2;
	typedef OdeSMOKE::MethodGear<denseOdeConstantPressure_TIO2> methodGearConstantPressure_TIO2;

	enum EnergyEquationFormulation { ENERGY_EQ_ENTHALPY, ENERGY_EQ_TEMPERATURE } ;
	enum DiffusivityModel { DIFFUSIVITY_MODEL_MIX_AVERAGED, DIFFUSIVITY_MODEL_LEWIS_NUMBERS } ;

	class OpenSMOKEppReactingMixture
		:	public IOdictionary
	{

	protected:

		// OpenSMOKE++ Maps

		//- Thermodynamics map
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*		thermodynamicsMap_; 

		//- Kinetics map
		OpenSMOKE::KineticsMap_CHEMKIN* 		kineticsMap_;

		//- Transport properties map
		OpenSMOKE::TransportPropertiesMap_CHEMKIN* 	transportMap_;

		//- Polimi Soot Model (based on the discrete sectional method)
		autoPtr<PolimiSootModel> polimiSootModel_;

		//- HMOM Model
		autoPtr<HMOMModel> hmomModel_;

		//- DRG Model
		autoPtr<DRG> drgModel_;

		//- Clustering Model
		autoPtr<Clustering> clusteringModel_;

		//- Repairing Model
		autoPtr<Repairing> repairingModel_;

		//- Repairing Model
		autoPtr<TiO2> tio2Model_;

		//- Chemistry Linear Model
		#if STEADYSTATE == 1
		autoPtr<ChemistryLinearModel> chemistryLinearModel_;
		#endif

		//- OpenSMOKE++ ODE solver for stiff problems
		OdeSMOKE::MultiValueSolver<methodGearConstantPressure>* odeSolverConstantPressure_;
		OdeSMOKE::MultiValueSolver<methodGearConstantPressure_HMOM>* odeSolverConstantPressure_HMOM_;
		OdeSMOKE::MultiValueSolver<methodGearConstantPressure_DRG>* odeSolverConstantPressure_DRG_;
		OdeSMOKE::MultiValueSolver<methodGearConstantPressure_TIO2>* odeSolverConstantPressure_TIO2_;

		//- CVODE solver for stiff problems
		#if OPENSMOKE_USE_SUNDIALS == 1
		OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_CVODE *cvodeSolverConstantPressure_;
		OpenSMOKE::BatchReactorHomogeneousConstantPressure_HMOM_ODE_CVODE *cvodeSolverConstantPressure_HMOM_;
		OpenSMOKE::BatchReactorHomogeneousConstantPressure_DRG_ODE_CVODE *cvodeSolverConstantPressure_DRG_;
		OpenSMOKE::BatchReactorHomogeneousConstantPressure_TIO2_ODE_CVODE *cvodeSolverConstantPressure_TIO2_;
		#endif
		
		//- Batch reactor with constant pressure
		BatchReactorHomogeneousConstantPressure* batchReactorHomogeneousConstantPressure_;
		BatchReactorHomogeneousConstantPressure_HMOM* batchReactorHomogeneousConstantPressure_HMOM_;
		BatchReactorHomogeneousConstantPressure_DRG* batchReactorHomogeneousConstantPressure_DRG_;
		BatchReactorHomogeneousConstantPressure_TIO2* batchReactorHomogeneousConstantPressure_TIO2_;

		//- ODE parameters for solution of ODE systems
		OpenSMOKE::ODE_Parameters odeParameterBatchReactorHomogeneous_;


		// Protected Data

		//- Temperature [K]
		volScalarField T_;

		//- Pressure [Pa]
		volScalarField p_;

		//- Compressibility [s2/m2]
		volScalarField psi_;

		//- Mass fractions
		PtrList<volScalarField> Y_;

		//- Mass fractions
		PtrList<volScalarField> X_;

		//- Constant pressure specific heat [J/kg/K]
		volScalarField Cp_;

		//- Constant volume specific heat [J/kg/K]
		volScalarField Cpv_;

		//- Constant pressure specific heat of species [J/kg/K]
		PtrList<volScalarField> CpSpecies_;

		//- Sensible enthalpy of species [J/kg]
		PtrList<volScalarField> hsSpecies_;

		//- Sensible enthalpy [J/kg]
		volScalarField Hs_;

		//- Dynamic viscosity [kg/m/s]
		volScalarField mu_;

		//- Thermal conductivity [W/m/K]
		volScalarField lambda_;

		//- Mass diffusion coefficients [m2/s]
		PtrList<volScalarField> Dmix_;

		//- Thermal diffusion (Soret) coefficients [m2/s]
		PtrList<volScalarField> Dsoret_;

		//- Mass diffusion fluxes [kg/m2/s]
		PtrList<volVectorField> J_;

		//- Additional terms for correcting the diffusion fluxes in case the gradient of mole fractions is considered [kg/m3/s]
		PtrList<volScalarField> mwDiffusionCorrection_;

		//- Thermal diffusivity [m2/s]
		volScalarField alpha_;

		//- Molecular weight [kg/kmol]
		volScalarField MW_;

		//- Total concentration [kmol/m3]
		volScalarField cTot_;

		//- Correction diffusion flux [kg/m2/s]
		volVectorField Jc_;

		//- Correction diffusion flow [kg/s]
		surfaceScalarField phic_;

		//- Contribution to energy equation from mass diffusion fluxes [kg/s3/K]
		volVectorField Cps_Js_;

		//- Contribution to energy equation from mass diffusion fluxes [kg/s3]
		volVectorField hs_Js_;

		//- Explicit part of heat conduction in terms of enthalpy [kg/s3]
		volVectorField HsConduction_;

		//- Heat release [kg/m/s3]
		volScalarField Q_;

		//- Radiative loss [kg/m/s3]
		volScalarField Qrad_;

		//- CPU time for chemical step [ms]
		volScalarField cpuChemistry_;

		//- Number of species
		int ns_;

		//- Number of reactions
		int nr_;

		//- Number of cells
		int ncells_;

		//- Minimum allowed temperature [K]
		scalar Tmin_;

		//- Maximum allowed temperature [K]
		scalar Tmax_;

		//- Index of inert species (zero-based)
		int inertIndex_;

		//- List of species
		hashedWordList species_;

		//- List of species active flags
		List<bool> active_;

		//- List of molecular weights as dimensined scalars [kg/kmol]
		List<dimensionedScalar> mw_species_;

		//- Should the Soret effect be included in the species equations?
		Switch soret_effect_;

		//- Has the species a non-zero Soret diffusion coefficient?
		List<Switch> soret_list_;

		//- Transfer indices between global and Soret numbering
		List<label> soret_index_;

		//- Should the thermophoretic effect be included in the species equations?
		Switch thermophoretic_effect_;

		//- Should the dpdt term be included in the enthalpy equation?
		Switch dpdt_;

		//- Should the diffusion flux be evaulated via gradient of mole fraction?
		Switch diffusion_via_mole_fraction_gradient_;

		//- Should the correction of diffusion fluxes be applied?
		Switch diffusion_fluxes_correction_;

		//- Should the mass fraction closure be based on inert species?
		Switch closure_on_inert_;

		//- Should the enthalpy flux associated to the mass diffusion fluxes be included in the energy equation?
		Switch diffusion_in_energy_equation_;

		//- Should the conductive heat flux contribution in the enthalpy equation be written via gradient of temperature? 
		Switch conduction_via_temperature_gradient_;

		//- Should the viscous heating source_term be included in the energy equation?
		Switch viscous_heating_source_term_;

		//- Should the transport equation for local residence time be solved?
		Switch solveForTauEquation_;

		//- Should the transport equation for passive scalar be solved?
		Switch solveForCsiEquation_;

		//- Should the momentum equation be solved?
		Switch solveForMomentumEquation_;

		//- Should the energy equation be solved?
		Switch solveForEnergyEquation_;

		//- Should the species equations be solved?
		Switch solveForSpeciesEquations_;

		//- Should the chemical reactions in the gaseous phase be taken into account?
		Switch homogeneousReactions_;

		//- Minimum temperature below which no chemistry integration is carried out
		scalar direct_integration_minimum_temperature_for_chemistry_;

		//- Formation enthalpy of species
		Eigen::VectorXd h0_species_;

		//- List of Lewis numbers
		Eigen::VectorXd lewis_numbers_;

		//- Energy equation formulation
		EnergyEquationFormulation energyEquation_;

		//- Diffusivity model
		DiffusivityModel diffusivityModel_; 

		//- Is the HMOM Model turned on?
		Switch	hmom_;

		//- Is the Polimi Soot Model turned on?
		Switch polimiSoot_;

		//- Is the DRG Model turned on?
		Switch drg_;

		//- Is the Clustering Model turned on?
		Switch clustering_;

		//- Is the Repairing Model turned on?
		Switch repairing_;

		//- Is the TiO2 Model turned on?
		Switch tio2_;

		#if STEADYSTATE == 1

		//- Chemical source terms: implicit contribution
		PtrList<volScalarField> sourceImplicit_;

		//- Chemical source terms: explicit contribution
		PtrList<volScalarField> sourceExplicit_;

		//- Saved Jacobian matrix
		std::vector<Eigen::VectorXd> Jstored_;

		#endif


		// Protected Member Functions

		//- Construct as copy (not implemented)
		OpenSMOKEppReactingMixture(const OpenSMOKEppReactingMixture&);

		//- Calculate the transport and thermodynamic properties
		void calculate();

		//- Return the enthalpy/internal energy field boundary types by interrogating the temperature field boundary types
		wordList hBoundaryBaseTypes();

		//- Return the enthalpy/internal energy field boundary base types by interrogating the temperature field boundary types
		wordList hBoundaryTypes();

		//- Initialize the fields
		void init();

		//- Correct the enthalpy field boundaries
		void hBoundaryCorrection(volScalarField& h);

		//- Normalize the mass fractions after the chemical step and return theumbalance		
		double normalizeMassFractions(Eigen::VectorXd& omega_plus_temperature, const label celli);


	public:

		//- Construct from mesh
		OpenSMOKEppReactingMixture(const fvMesh& mesh);


		//- Destructor
		virtual ~OpenSMOKEppReactingMixture();

		//- Thermodynamics map
		OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap() { return *thermodynamicsMap_; }; 

		//- Kinetics map
		OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap() { return *kineticsMap_; }; 

		//- Transport properties map
		OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap() { return *transportMap_; }; 

        	//- Return false because the equation of state is compressible
		bool incompressible() const { return false; }


		//- Name of main dictionary
		static const word dictName;

		//- New function
		static autoPtr<Foam::OpenSMOKEppReactingMixture> New(const fvMesh& mesh);

		//- Update the thermodynamic and transport properties
		void update_properties();

		//- Update fluxes and transport terms in species/energy equations
		void update_transport_terms(const fvMesh& mesh, const volScalarField& rho);

		//- Solve additional equations
		void solve_additional_equations(const fvMesh& mesh, const Foam::volScalarField& rho, const Foam::surfaceScalarField& phi);

		//- Carry out the chemical step via Direct Integration
		void chemistry_direct_integration(const double t0, const double tf, const fvMesh& mesh, const Foam::volScalarField& rho);

		//- Post-processing operations (on-th-fly)
		void onthefly_postprocessing(const fvMesh& mesh, const Foam::volScalarField& rho);

		#if STEADYSTATE == 1
		
		//- Update the chemical source terms
		void update_chemical_source_terms();

		//- Update the mixture properties
		void update_properties_for_steady_state();

		//- Update the sequence of species to be solved according to the defined policy
		void update_species_order_policy();
		
		//- Returns the index of species to be solved according to the defined policy
		label species_order(const label i) const;
		
		#endif

		//- Correct the density
		void correctRho(const volScalarField& deltaRho) { }

		//- Return the number of species
		int ns() const { return ns_; }

		//- Returns the number of reactions
		int nr() const { return nr_; }

		//- Returns the minimum allowed temperature [K]
		int Tmin() const { return Tmin_; }

		//- Returns the maximum allowed temperature [K]
		int Tmax() const { return Tmax_; }

		//- Return the index of inert species (zero-based)
		int inertIndex() const { return inertIndex_; }

		//- Should the Soret effect be included in the species equations?
		Switch soret_effect() const { return soret_effect_; }

		//- Has the species a non-zero Soret diffusion coefficient?
		Switch soret_list(const label i) const { return soret_list_[i]; }

		//- Transfer indices between global and Soret numbering
		label soret_index(const label i) const { return soret_index_[i]; }

		//- Should the thermophoretic effect be included in the species equations?
		Switch thermophoretic_effect() const { return thermophoretic_effect_; }

		//- Should the dpdt term be included in the enthalpy equation?
		Switch dpdt() const { return dpdt_; }

		//- Should the diffusion flux be evaulated via gradient of mole fraction?
		Switch diffusion_via_mole_fraction_gradient() const { return diffusion_via_mole_fraction_gradient_; }

		//- Should the correction of diffusion fluxes be applied?
		Switch diffusion_fluxes_correction() const { return diffusion_fluxes_correction_; }

		//- Should the mass fraction closure be based on inert species?
		Switch closure_on_inert() const { return closure_on_inert_; }

		//- Should the enthalpy flux associated to the mass diffusion fluxes be included in the energy equation?
		Switch diffusion_in_energy_equation() const { return diffusion_in_energy_equation_; }

		//- Should the conductive heat flux contribution in the enthalpy equation be written via gradient of temperature? 
		Switch conduction_via_temperature_gradient() const { return conduction_via_temperature_gradient_; }

		//- Should the viscous heating source_term be included in the energy equation?
		Switch viscous_heating_source_term() const { return viscous_heating_source_term_; }

		//- Should the transport equation for local residence time be solved?
		Switch solveForTauEquation() const { return solveForTauEquation_; }

		//- Should the transport equation for passive scalar be solved?
		Switch solveForCsiEquation() const { return solveForCsiEquation_; }

		//- Should the transport equation for momentum be solved?
		Switch solveForMomentumEquation() const { return solveForMomentumEquation_; }

		//- Should the transport equation for energy be solved?
		Switch solveForEnergyEquation() const { return solveForEnergyEquation_; }

		//- Should the transport equations for species be solved?
		Switch solveForSpeciesEquations() const { return solveForSpeciesEquations_; }

		//- Energy equation formulation
		EnergyEquationFormulation energyEquation() const { return energyEquation_; }

		//- Diffusivity model
		DiffusivityModel diffusivityModel() const { return diffusivityModel_; }

		//- Is the HMOM model turned on?
		Switch hmom() const { return hmom_; }

		//- Is the Polimi Soot Model turned on?
		Switch polimiSoot() const { return polimiSoot_; }

		//- Is the DRG Model turned on?
		Switch drg() const { return drg_; }

		//- Is the Clustering Model turned on?
		Switch clustering() const { return clustering_; }

		//- Is the Repairing Model turned on?
		Switch repairing() const { return repairing_; }

		//- Is the TiO2 model turned on?
		Switch tio2() const { return tio2_; }

		// Member Functions

		//- Return the list of species
		inline const hashedWordList& species() const;

		//- Does the mixture include this specie?
		inline bool contains(const word& specieName) const;

		//- Return true for active species
		inline bool active(label speciei) const;

		//- Return the bool list of active species
		inline const List<bool>& active() const;

		//- Set speciei active
		inline void setActive(label speciei);

		//- Set speciei inactive
		inline void setInactive(label speciei);

		
		//- Return the mass fraction fields
		inline PtrList<volScalarField>& Y();

		//- Return the const mass fraction fields
		inline const PtrList<volScalarField>& Y() const;

		//- Return the mass fraction field for a specie given by index
		inline volScalarField& Y(const label i);

		//- Return the const mass fraction field for a specie given by index
		inline const volScalarField& Y(const label i) const;

		//- Return the mass fraction field for a specie given by name
		inline volScalarField& Y(const word& specieName);

		//- Return the const mass fraction field for a specie given by name
		inline const volScalarField& Y(const word& specieName) const;

		
		//- Return the mole fraction fields
		inline PtrList<volScalarField>& X();

		//- Return the const mole fraction fields
		inline const PtrList<volScalarField>& X() const;

		//- Return the mole fraction field for a specie given by index
		inline volScalarField& X(const label i);

		//- Return the const mole fraction field for a specie given by index
		inline const volScalarField& X(const label i) const;

		//- Return the mole fraction field for a specie given by name
		inline volScalarField& X(const word& specieName);

		//- Return the const mole fraction field for a specie given by name
		inline const volScalarField& X(const word& specieName) const;


		//- Return the mass diffusion coefficient fields
		inline PtrList<volScalarField>& Dmix();

		//- Return the const mass diffusion coefficient fields
		inline const PtrList<volScalarField>& Dmix() const;

		//- Return the mass diffusion coefficient field for a specie given by index
		inline volScalarField& Dmix(const label i);

		//- Return the const mass diffusion coefficient field for a specie given by index
		inline const volScalarField& Dmix(const label i) const;

		//- Return the mass diffusion coefficient field for a specie given by name
		inline volScalarField& Dmix(const word& specieName);

		//- Return the const mass diffusion coefficient field for a specie given by name
		inline const volScalarField& Dmix(const word& specieName) const;


		//- Return the thermal diffusion (Soret) coefficient fields
		inline PtrList<volScalarField>& Dsoret();

		//- Return the const mass diffusion (Soret) coefficient fields
		inline const PtrList<volScalarField>& Dsoret() const;


		//- Return the mass diffusion flux fields
		inline PtrList<volVectorField>& J();

		//- Return the const mass diffusion flux fields
		inline const PtrList<volVectorField>& J() const;

		//- Return the mass diffusion flux for a specie given by index
		inline volVectorField& J(const label i);

		//- Return the const mass diffusion flux for a specie given by index
		inline const volVectorField& J(const label i) const;

		//- Return the mass diffusion flux for a specie given by name
		inline volVectorField& J(const word& specieName);

		//- Return the const mass diffusion flux field for a specie given by name
		inline const volVectorField& J(const word& specieName) const;

		#if STEADYSTATE == 1

		//- Return the implicit chemical source term fields
		inline PtrList<volScalarField>& sourceImplicit();

		//- Return the const implicit chemical source term fields
		inline const PtrList<volScalarField>& sourceImplicit() const;

		//- Return the implicit chemical source term field by index
		inline volScalarField& sourceImplicit(const label i);

		//- Return the const implicit chemical source term field by index
		inline const volScalarField& sourceImplicit(const label i) const;

		//- Return the explicit chemical source term fields
		inline PtrList<volScalarField>& sourceExplicit();

		//- Return the const explicit chemical source term fields
		inline const PtrList<volScalarField>& sourceExplicit() const;

		//- Return the explicit chemical source term field by index
		inline volScalarField& sourceExplicit(const label i);

		//- Return the const explicit chemical source term field by index
		inline const volScalarField& sourceExplicit(const label i) const;

		#endif

		//- Return the constant pressure specific heat of species fields
		inline PtrList<volScalarField>& CpSpecies();

		//- Return the const constant pressure specific heat of species fields
		inline const PtrList<volScalarField>& CpSpecies() const;

		//- Return the constant pressure specific heat of species field for a specie given by index
		inline volScalarField& CpSpecies(const label i);

		//- Return the const constant pressure specific heat of species field for a specie given by index
		inline const volScalarField& CpSpecies(const label i) const;

		//- Return the constant pressure specific heat of species field for a specie given by name
		inline volScalarField& CpSpecies(const word& specieName);

		//- Return the const constant pressure specific heat of species field for a specie given by name
		inline const volScalarField& CpSpecies(const word& specieName) const;


		//- Return the enthalpy of species fields
		inline PtrList<volScalarField>& hsSpecies();

		//- Return the enthalpy species fields
		inline const PtrList<volScalarField>& hsSpecies() const;

		//- Return the enthalpy field for a specie given by index
		inline volScalarField& hsSpecies(const label i);

		//- Return the enthalpy of species field for a specie given by index
		inline const volScalarField& hsSpecies(const label i) const;

		//- Return the enthalpy of species field for a specie given by name
		inline volScalarField& hsSpecies(const word& specieName);

		//- Return the enthalpy of species field for a specie given by name
		inline const volScalarField& hsSpecies(const word& specieName) const;


		//- Return the correction on the diffusion fluxes associated to the molecular weight
		inline PtrList<volScalarField>& mwDiffusionCorrection();

		//- Return the const correction on the diffusion fluxes associated to the molecular weight
		inline const PtrList<volScalarField>& mwDiffusionCorrection() const;

		//- Return the correction on the diffusion fluxes associated to the molecular weight
		inline volScalarField& mwDiffusionCorrection(const label i);

		//- Return the const correction on the diffusion fluxes associated to the molecular weight
		inline const volScalarField& mwDiffusionCorrection(const label i) const;

		//- Return the correction on the diffusion fluxes associated to the molecular weight
		inline volScalarField& mwDiffusionCorrection(const word& specieName);

		//- Return the const correction on the diffusion fluxes associated to the molecular weight
		inline const volScalarField& mwDiffusionCorrection(const word& specieName) const;


		//- Density [kg/m3] - uses current value of pressure
		virtual tmp<volScalarField> rho() const;

		//- Density for patch [kg/m^3]
		virtual tmp<scalarField> rho(const label patchi) const;

		//- Pressure [Pa]
		//  Non-const access allowed for transport equations
		virtual volScalarField& p();

		//- Pressure [Pa]
		virtual const volScalarField& p() const;

		//- Temperature [K]
		virtual const volScalarField& T() const;

		//- Temperature [K]
		//  Non-const access allowed for transport equations
		virtual volScalarField& T();

		//- Compressibility [s^2/m^2]
		virtual const volScalarField& psi() const;

		//- Constant pressure specific heat of mixture [J/kg/K]
		virtual tmp<volScalarField> Cp() const;

		//- Constant pressure specific heat of mixture for patch [J/kg/K]
		virtual tmp<scalarField> Cp(const label patchi) const;

		//- Constant volume specific heat of mixture [J/kg/K]
		virtual tmp<volScalarField> Cpv() const;

		//- Constant volume specific heat of mixture for patch [J/kg/K]
		virtual tmp<scalarField> Cpv(const label patchi) const;

		//- Sensible enthalpy of mixture [J/kg]
		virtual const volScalarField& Hs() const;

		//- Sensible enthalpy of mixture [J/kg]
		//  Non-const access allowed for transport equations
		virtual volScalarField& Hs();

		//- Sensible enthalpy of mixture for patch [J/kg]
		virtual tmp<scalarField> Hs(const label patchi) const;

		//- Dynamic viscosity of mixture [kg/m/s]
		virtual tmp<volScalarField> mu() const;

		//- Dynamic viscosity of mixture for patch [kg/m/s]
		virtual tmp<scalarField> mu(const label patchi) const;

		//- Thermal conductivity of mixture [W/m/K]
		virtual tmp<volScalarField> lambda() const;

		//- Thermal conductivity of mixture for patch [W/m/K]
		virtual tmp<scalarField> lambda(const label patchi) const;

		//- Thermal diffusivity of mixture [m2/s]
		virtual const volScalarField& alpha() const;

		//- Thermal diffusivity of mixture for patch [m2/s]
		virtual const scalarField& alpha ( const label patchi ) const;

		//- Molecular weight of mixture [kg/kmol]
		virtual const volScalarField& MW() const;

		//- Molecular weight of mixture for patch [kg/kmol]
		virtual const scalarField& MW ( const label patchi ) const;

		//- Total concentration of mixture [kmol/m3]
		virtual const volScalarField& cTot() const;

		//- Total concentration of mixture for patch [kmol/m3]
		virtual const scalarField& cTot ( const label patchi ) const;

		//- Mass diffusion correction flux [kg/m2/s]
		//  Non-const access allowed for transport equations
		virtual volVectorField& Jc();

		//- Mass diffusion correction flux for patch [kg/m2/s]
		virtual tmp<vectorField> Jc(const label patchi) const;

		//- Correction diffusion flow [kg/s]
		//  Non-const access allowed for transport equations
		virtual surfaceScalarField& phic();

		//- Correction diffusion flow [kg/s]
		virtual tmp<scalarField> phic(const label patchi) const;

		//- Contribution to energy equation from mass diffusion fluxes [kg/s3/K]
		//  Non-const access allowed for transport equations
		virtual volVectorField& Cps_Js();

		//- Contribution to energy equation from mass diffusion fluxes [kg/s3/K]
		virtual tmp<vectorField> Cps_Js(const label patchi) const;

		//- Contribution to energy equation from mass diffusion fluxes [kg/s3]
		//  Non-const access allowed for transport equations
		virtual volVectorField& hs_Js();

		//- Contribution to energy equation from mass diffusion fluxes [kg/s3]
		virtual tmp<vectorField> hs_Js(const label patchi) const;

		//- Explicit part of heat conduction in terms of enthalpy [kg/s3]
		//  Non-const access allowed for transport equations
		virtual volVectorField& HsConduction();

		//- Explicit part of heat conduction in terms of enthalpy [kg/s3]
		virtual tmp<vectorField> HsConduction(const label patchi) const;

		//- Reaction heat [kg/m/s3]
		virtual const volScalarField& Q() const;

		//- Reaction heat for patch [kg/m/s3]
		virtual const scalarField& Q ( const label patchi ) const;

		//- Radiation heat loss [kg/m/s3]
		virtual volScalarField& Qrad();

		//- CPU time for integration of chemical step [s]
		virtual const volScalarField& cpuChemistry() const;

		//- CPU time for integration of chemical step [s]
		virtual const scalarField& cpuChemistry ( const label patchi ) const;


		//- Return the mass fraction fields
		const dimensionedScalar& mw_species(const label i) { return mw_species_[i]; }


		//- Sensible enthalpy for cell-set [J/kg]
             	virtual tmp<scalarField> hs ( const scalarField& p, const scalarField& T, const labelList& cells ) const;

		//- Sensible enthalpy for patch [J/kg]
             	virtual tmp<scalarField> hs ( const scalarField& p, const scalarField& T, const label patchi ) const;

		//- Heat capacity at constant pressure/volume for patch [J/kg/K]
             	virtual tmp<scalarField> Cpv ( const scalarField& p, const scalarField& T, const label patchi ) const;
	};


} // End namespace Foam

#include "OpenSMOKEppReactingMixtureI.H"

#include "OpenSMOKEppReactingMixture.HPP"

#endif // OpenSMOKEppReactingMixture_H


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
|   Copyright(C) 2021 Alberto Cuoci                                       |
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
const Foam::word Foam::PCA::dictName("PCA");

// Constructor
Foam::PCA::PCA(const Foam::fvMesh& mesh)
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
	)
{
	// Set accuracy coefficients (mandatory)
	n_ = readLabel(this->lookup("numberPCs"));

	// Enable calculation of PCs on boundary faces (for graphical reasons only)
	boundaryFaces_ = lookupOrDefault<Switch>("boundaryFaces", false);

	// Set the algorithm to be used (mandatory)
	{
		word samplingType = lookup("sampling");
		if (samplingType == "SimpleStep")
			samplingType_ = SAMPLING_SIMPLE_STEP;
		else FatalError << "PCA: wrong sampling argument. Available options: SimpleStep" << ::Foam::exit(FatalError);
	}

	// Set sampling step
	samplingStep_ = readLabel(this->lookup("samplingStep"));

	// Set the verbosity level to be used (optional, default: 1)
	verbosity_ = lookupOrDefault<label>("verbosity", 1);

	// Set the list of enabled species
	{
		List<word> listSpecies(this->lookup("species"));
		enabled_species_names_ = listSpecies;
	}
}

Foam::PCA::~PCA()
{}

void Foam::PCA::Setup(const Foam::fvMesh& mesh, const OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap)
{
	// List of species enabled for PCA
	if (enabled_species_names_.size() == 1 && ( enabled_species_names_[0] == "ALL" || enabled_species_names_[0] == "all") )
	{
		enabled_species_indices_.resize(thermodynamicsMap.NumberOfSpecies());
		for (unsigned int i=0;i<thermodynamicsMap.NumberOfSpecies();i++)
			enabled_species_indices_[i] = i;
	}
	else
	{
		enabled_species_indices_.resize(enabled_species_names_.size());
		for (unsigned int i=0;i<enabled_species_indices_.size();i++)
			enabled_species_indices_[i] = thermodynamicsMap.IndexOfSpecies(enabled_species_names_[i])-1;
	}

	// Summary on the screen
	{
		Info << endl;
		Info << "List of species included in PCA:" << endl;
		for (unsigned int i=0;i<enabled_species_indices_.size();i++)
			Info << " * " << enabled_species_names_[i] << " (" << enabled_species_indices_[i] << ")" << endl;
		Info << endl;
	}

	// Allocate memory for PCs fields
	{
		Info << "Creating PC fields..." << endl;

		pc_.resize(n_);
		for (unsigned int i=0;i<n_;i++)
		{
			std::stringstream label;
			label << i+1;
			const std::string name = "pc:" + label.str();

			pc_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						name,
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					mesh,
					dimensionedScalar("pc", dimensionSet(0, 0, 0, 0, 0), 0.)
				)
			);
		}
	}
}

void Foam::PCA::Analysis(const volScalarField& T, const PtrList<volScalarField>& Y, Eigen::VectorXd& explained, Eigen::MatrixXd& cosines)
{
	const scalarField& TCells = T.internalField();

	const label ns = Y.size();				// total number of species
	const label nc = enabled_species_indices_.size();	// total number of enabled species

	// Calculation of observations
	const label ncells = T.size();
	const int modulus = ncells % samplingStep_;
	const int nobs = (ncells-modulus)/samplingStep_ + (modulus == 0? 0 : 1);

	Eigen::MatrixXd X(nobs, nc+1);				// observation matrix
	Eigen::VectorXd x(nc+1);				// single observation

	unsigned int iobs=0;
	forAll(TCells, celli)
	{	
		if (celli % samplingStep_ == 0)
		{
			x(0) = TCells[celli];
			for(unsigned int i=0;i<nc;i++)
				x(i+1) = Y[enabled_species_indices_[i]].internalField()[celli];

			X.row(iobs) = x;

			iobs++;
		}
	}

	if (nobs != iobs)
		FatalError << "PCAModel: something wrong in the application of the sampling procedure" << ::Foam::exit(FatalError);

	// PCA Class
	PCAEigen* pca = new PCAEigen();

	// Perform PCA using the SVD approach
	double tStartPCA = OpenSMOKE::OpenSMOKEGetCpuTime();
	{
		const int result = pca->Calculate(X);

		if (result != 0) 
			FatalError << "PCAModel: There is an error during PCA calculation!" << ::Foam::exit(FatalError);
	}
	double tEndPCA = OpenSMOKE::OpenSMOKEGetCpuTime();

	// Explained variance
	explained = pca->explained();

	// Cosines 
	pca->CalculateCosines();
	cosines = pca->cosines();
}

void Foam::PCA::Analysis(const volScalarField& T, const PtrList<volScalarField>& Y)
{		
	const scalarField& TCells = T.internalField();

	const label ns = Y.size();				// total number of species
	const label nc = enabled_species_indices_.size();	// total number of enabled species

	// Calculation of observations
	const label ncells = T.size();
	const int modulus = ncells % samplingStep_;
	const int nobs = (ncells-modulus)/samplingStep_ + (modulus == 0? 0 : 1);

	//Info << endl;
	//Info << "   PCA: " << "ncells: " << ncells << " sampling: " << samplingStep_ << " observations: " << nobs << " modulus: " << modulus << endl;

	Eigen::MatrixXd X(nobs, nc+1);				// observation matrix
	Eigen::VectorXd x(nc+1);				// single observation

	unsigned int iobs=0;
	forAll(TCells, celli)
	{	
		if (celli % samplingStep_ == 0)
		{
			x(0) = TCells[celli];
			for(unsigned int i=0;i<nc;i++)
				x(i+1) = Y[enabled_species_indices_[i]].internalField()[celli];

			X.row(iobs) = x;

			iobs++;
		}
	}

	if (nobs != iobs)
		FatalError << "PCAModel: something wrong in the application of the sampling procedure" << ::Foam::exit(FatalError);

	// PCA Class
	PCAEigen* pca = new PCAEigen();

	// Perform PCA using the SVD approach
	double tStartPCA = OpenSMOKE::OpenSMOKEGetCpuTime();
	{
		const int result = pca->Calculate(X);

		if (result != 0) 
			FatalError << "PCAModel: There is an error during PCA calculation!" << ::Foam::exit(FatalError);
	}
	double tEndPCA = OpenSMOKE::OpenSMOKEGetCpuTime();

	// Explained variance
	Info << "Principal Component Analysis..." << endl;

	Info << " * Explained variance (n=1):        " << pca->explained()(0) << endl;
	Info << " * Explained variance (n=2):        " << pca->explained()(1) << endl;
	Info << " * Explained variance (n=3):        " << pca->explained()(2) << endl;
	Info << " * Explained variance (n=4):        " << pca->explained()(3) << endl;

	Info << " * Explained variance (95%):        " << pca->explained_095() << endl;
	Info << " * Explained variance (98%):        " << pca->explained_098() << endl;
	Info << " * Explained variance (99%):        " << pca->explained_099() << endl;

	Info << " * CPU time (s):                    " << tEndPCA-tStartPCA << endl;
	Info << " * CPU time per cell (ms):          " << (tEndPCA-tStartPCA)/double(ncells)*1000. << endl;

	// Evaluation of scores
	double tStartScores = OpenSMOKE::OpenSMOKEGetCpuTime();
	{
		Eigen::VectorXd xt(nc+1);		// single observation
		Eigen::VectorXd scores(n_);		// single score vector
	
		// Internal cells
		forAll(TCells, celli)
		{	
			xt(0) = TCells[celli];
			for(unsigned int i=0;i<nc;i++)
				xt(i+1) = Y[enabled_species_indices_[i]].internalField()[celli];
		
			pca->Scores(xt, n_, scores);

			for(unsigned int i=0;i<n_;i++)
				pc_[i].ref()[celli] = scores(i);
		}

		// Boundary faces (for graphical purposes only)
		if (boundaryFaces_ == true)
		{
			const volScalarField::Boundary& TBf = T.boundaryField();

			forAll(T.boundaryField(), patchi)
			{
				
				const fvPatchScalarField& pT = TBf[patchi];
	
				forAll(pT, facei)
				{
					xt(0) = pT[facei];
					for(unsigned int i=0;i<nc;i++)
						xt(i+1) = Y[enabled_species_indices_[i]].boundaryField()[patchi][facei];
		
					pca->Scores(xt, n_, scores);

					for(unsigned int i=0;i<n_;i++)
						pc_[i].boundaryFieldRef()[patchi][facei] = scores(i);
				}
			}
		}
	}
	double tEndScores = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << " * CPU time (scores)(s):            " << tEndScores-tStartScores << endl;
	Info << " * CPU time per cell (scores) (ms): " << (tEndScores-tStartScores)/double(ncells)*1000. << endl;
}



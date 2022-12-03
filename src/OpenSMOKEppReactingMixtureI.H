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
|   Copyright(C) 2020, 2021 Alberto Cuoci                                 |
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

inline Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::Y()
{
	return Y_;
}
  
inline const Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::Y() const
{
	return Y_;
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::Y(const label i)
{
	return Y_[i];
}
  
inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::Y( const label i ) const
{
	return Y_[i];
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::Y( const word& speciesName )
{
	return Y_[species_[speciesName]];
}

inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::Y ( const word& speciesName ) const
{
	return Y_[species_[speciesName]];
}


inline Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::X()
{
	return X_;
}
  
inline const Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::X() const
{
	return X_;
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::X(const label i)
{
	return X_[i];
}
  
inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::X( const label i ) const
{
	return X_[i];
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::X( const word& speciesName )
{
	return X_[species_[speciesName]];
}

inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::X ( const word& speciesName ) const
{
	return X_[species_[speciesName]];
}


inline Foam::PtrList<Foam::volVectorField>& Foam::OpenSMOKEppReactingMixture::J()
{
	return J_;
}
  
inline const Foam::PtrList<Foam::volVectorField>& Foam::OpenSMOKEppReactingMixture::J() const
{
	return J_;
}

inline Foam::volVectorField& Foam::OpenSMOKEppReactingMixture::J(const label i)
{
	return J_[i];
}
  
inline const Foam::volVectorField& Foam::OpenSMOKEppReactingMixture::J( const label i ) const
{
	return J_[i];
}

inline Foam::volVectorField& Foam::OpenSMOKEppReactingMixture::J( const word& speciesName )
{
	return J_[species_[speciesName]];
}

inline const Foam::volVectorField& Foam::OpenSMOKEppReactingMixture::J ( const word& speciesName ) const
{
	return J_[species_[speciesName]];
}

#if STEADYSTATE == 1

inline Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::sourceImplicit()
{
	return sourceImplicit_;
}
  
inline const Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::sourceImplicit() const
{
	return sourceImplicit_;
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::sourceImplicit(const label i)
{
	return sourceImplicit_[i];
}
  
inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::sourceImplicit( const label i ) const
{
	return sourceImplicit_[i];
}

inline Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::sourceExplicit()
{
	return sourceExplicit_;
}
  
inline const Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::sourceExplicit() const
{
	return sourceExplicit_;
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::sourceExplicit(const label i)
{
	return sourceExplicit_[i];
}
  
inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::sourceExplicit( const label i ) const
{
	return sourceExplicit_[i];
}

#endif

inline Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::Dmix()
{
	return Dmix_;
}

inline const Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::Dmix() const
{
	return Dmix_;
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::Dmix(const label i)
{
	return Dmix_[i];
}

inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::Dmix( const label i ) const
{
	return Dmix_[i];
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::Dmix( const word& speciesName )
{
	return Dmix_[species_[speciesName]];
}

inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::Dmix( const word& speciesName ) const
{
	return Dmix_[species_[speciesName]];
}


inline Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::Dsoret()
{
	return Dsoret_;
}

inline const Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::Dsoret() const
{
	return Dsoret_;
}


inline Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::CpSpecies()
{
	return CpSpecies_;
}

inline const Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::CpSpecies() const
{
	return CpSpecies_;
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::CpSpecies(const label i)
{
	return CpSpecies_[i];
}

inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::CpSpecies( const label i ) const
{
	return CpSpecies_[i];
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::CpSpecies( const word& speciesName )
{
	return CpSpecies_[species_[speciesName]];
}

inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::CpSpecies( const word& speciesName ) const
{
	return CpSpecies_[species_[speciesName]];
}


inline Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::hsSpecies()
{
	return hsSpecies_;
}

inline const Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::hsSpecies() const
{
	return hsSpecies_;
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::hsSpecies(const label i)
{
	return hsSpecies_[i];
}

inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::hsSpecies( const label i ) const
{
	return hsSpecies_[i];
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::hsSpecies( const word& speciesName )
{
	return hsSpecies_[species_[speciesName]];
}

inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::hsSpecies( const word& speciesName ) const
{
	return hsSpecies_[species_[speciesName]];
}


inline Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::mwDiffusionCorrection()
{
	return mwDiffusionCorrection_;
}

inline const Foam::PtrList<Foam::volScalarField>& Foam::OpenSMOKEppReactingMixture::mwDiffusionCorrection() const
{
	return mwDiffusionCorrection_;
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::mwDiffusionCorrection(const label i)
{
	return mwDiffusionCorrection_[i];
}

inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::mwDiffusionCorrection( const label i ) const
{
	return mwDiffusionCorrection_[i];
}

inline Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::mwDiffusionCorrection( const word& speciesName )
{
	return mwDiffusionCorrection_[species_[speciesName]];
}

inline const Foam::volScalarField& Foam::OpenSMOKEppReactingMixture::mwDiffusionCorrection( const word& speciesName ) const
{
	return mwDiffusionCorrection_[species_[speciesName]];
}


inline const Foam::hashedWordList& Foam::OpenSMOKEppReactingMixture::species() const
{
	return species_;
}

inline bool Foam::OpenSMOKEppReactingMixture::contains( const word& specieName ) const
{
	return species_.found(specieName);
}

inline bool Foam::OpenSMOKEppReactingMixture::active(label speciei) const
{
	return active_[speciei];
}

inline const Foam::List<bool>& Foam::OpenSMOKEppReactingMixture::active() const
{
	return active_;
}

inline void Foam::OpenSMOKEppReactingMixture::setActive(label speciei)
{
	active_[speciei] = true;
}

inline void Foam::OpenSMOKEppReactingMixture::setInactive(label speciei)
{
	active_[speciei] = false;
}


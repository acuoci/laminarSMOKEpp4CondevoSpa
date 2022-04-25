#pragma once

#include "ClusteringAlgorithms.H"

class DMZAlgorithm : public ClusteringAlgorithms
{
public:

	/**
		*@brief Default constructor
	*/
	DMZAlgorithm();

	/**
	*@brief Set the number of features
	*@param D number of features
	*/
	virtual void SetNumberOfFeatures(const unsigned int D);

	/**
	*@brief Set the vector of data corresponding to one of the features on which to cluster
	*@param index index of the feature (from 0)
	*@param psi vector of data
	*@param epsilon accuracy level (absolute value)
	*@param relative flag to specify if accuracy is given in relative (true) or absolute (false) units
	*/
	virtual void SetData(const unsigned int index, const std::vector<double>& psi, const double epsilon, const bool relative);

	/**
	*@brief Cluster the data provided by the user via SetData function
	*/
	virtual void ClusterData();


protected:

	/**
	*@brief Check the data provided by the user
	*/
	virtual void CheckInputData();

	/**
	*@brief Fatal error message
	*/
	virtual void FatalErrorMessage(const std::string message);


private:

	unsigned int nc_init_;						//!< initial number of clusters along the first feature (default 5)
	std::vector<std::vector<double>> y0_;		//!< values first (main feature)
	std::vector<std::vector<int>> ljloc_;		//!< local index (temporary)
};


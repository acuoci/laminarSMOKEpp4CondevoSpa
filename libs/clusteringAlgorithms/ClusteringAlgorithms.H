#pragma once

#include <string>
#include <vector>

/**
*@brief Calculation of centroid of a vector of values (i.e. the arithmetic mean)
*@param y values
*@return the centroid
*/
double centroid(const std::vector<double>& y);

/**
*@brief Calculation of dispersion coefficient of a vector of values (i.e. the maximum distance between 2 elements)
*@param y values
*@return the dispersion coefficient
*/
double dispersion_coefficient(const std::vector<double>& y);

class ClusteringAlgorithms
{
public:

	/**
	*@brief Default constructor
	*/
	ClusteringAlgorithms();

	/**
	*@brief Set the number of features
	*@param D number of features
	*/
	void SetNumberOfFeatures(const unsigned int D);

	/**
	*@brief Set the verbosity level
	*@param verbosity level: 0=minimum
	*/
	void SetVerbosityLevel(const unsigned int verbosity);

	/**
	*@brief Set the vector of data corresponding to one of the features on which to cluster
	*@param index index of the feature (from 0)
	*@param psi vector of data
	*@param epsilon accuracy level (absolute value)
	*@param relative flag to specify if accuracy is given in relative (true) or absolute (false) units
	*/
	virtual void SetData(const unsigned int index, const std::vector<double>& psi, const double epsilon, const bool relative) = 0;

	/**
	*@brief Cluster the data provided by the user via SetData function
	*/
	virtual void ClusterData() = 0;

	/**
	*@brief Returns the clustering's results
	*@param nc number of clusters
	*@param gj vector of global indices
	*@param lj vector of local indices
	*/
	void GetData(unsigned int& nc, std::vector<int>& gj, std::vector< std::vector<int> >& lj);

	/**
	*@brief Calculate the clusters' centroids
	*/
	void CalculateCentroids();

	/**
	*@brief Calculate the clusters' dispersion coefficients
	*/
	void CalculateDispersion();

	/**
	*@brief Print a summary of relevant data about the clustering
	*/
	void Print();


protected:

	/**
	*@brief Check the data provided by the user
	*/
	virtual void CheckInputData() = 0;

	/**
	*@brief Fatal error message
	*/
	virtual void FatalErrorMessage(const std::string message) = 0;


protected:

	unsigned int D_;							//!< number of features
	unsigned int n_;							//!< number of points to be clustered
	unsigned int nc_;							//!< number of clusters
	std::vector<double>	epsilon_;				//!< tolerance along each direction (absolute value)
	std::vector<int> ne_;						//!< number of elements in a cluster

	std::vector<std::vector<double>> C_;		//!< centroids for each cluster
	std::vector<std::vector<double>> sigma_;	//!< dispersion coefficients for each cluster
	std::vector<double>	C_avg_;					//!< average centroid
	std::vector<double>	sigma_avg_;				//!< average dispersion coefficient

	std::vector<int> gj_;						//!< global index
	std::vector<std::vector<int>> lj_;			//!< local index

	std::vector< std::vector<double> > psi_;	//!< features' values

	double cpu_clustering_;						//!< cpu time for clustering (in s)
	double cpu_centroids_;						//!< cpu time for centroid calculation (in s)
	double cpu_dispersion_;						//!< cpu time for dispersion coefficient calculation (in s)

	int verbosity_;								//!< verbosity level
};

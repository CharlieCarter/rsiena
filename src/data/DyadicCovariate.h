/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: DyadicCovariate.h
 * 
 * Description: This file contains the implementation of the
 * DyadicCovariate class.
 *****************************************************************************/

#ifndef DYADICCOVARIATE_H_
#define DYADICCOVARIATE_H_

#include "utils/NamedObject.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ActorSet;


// ----------------------------------------------------------------------------
// Section: DyadicCovariate class
// ----------------------------------------------------------------------------

/**
 * The base class for dyadic covariates. A dyadic covariate assigns double
 * values to pairs of actors, where the first and the second actor in a pair
 * may belong to different actor sets.
 */
class DyadicCovariate : public NamedObject
{
public:
	DyadicCovariate(std::string name,
		const ActorSet * pFirstActorSet,
		const ActorSet * pSecondActorSet);
	const ActorSet * pFirstActorSet() const;
	const ActorSet * pSecondActorSet() const;
	
	inline double mean() const;
	void mean(double value);
	inline bool hasLogTransformedMean() const;
	inline double logTransformedMean() const;
	void logTransformedMean(double value);
	inline bool hasAsinhTransformedMean() const;
	inline double asinhTransformedMean() const;
	void asinhTransformedMean(double value);
	
private:
	// One of the involved actor sets
	const ActorSet * lpFirstActorSet;
	
	// The other of the involved actor sets
	const ActorSet * lpSecondActorSet;
	
	// The average covariate value
	double lmean {};

	// Optional externally supplied means on named transformed scales. Keeping
	// these separate prevents a fitted log centre from being reused by an
	// asinh-transformed effect (or vice versa).
	bool lhasLogTransformedMean {};
	double llogTransformedMean {};
	bool lhasAsinhTransformedMean {};
	double lasinhTransformedMean {};
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the average value of this covariate.
 */
double DyadicCovariate::mean() const
{
	return this->lmean;
}

/**
 * Returns whether an externally supplied transformed-scale mean is present.
 */
bool DyadicCovariate::hasLogTransformedMean() const
{
	return this->lhasLogTransformedMean;
}

/**
 * Returns the externally supplied transformed-scale mean.
 */
double DyadicCovariate::logTransformedMean() const
{
	return this->llogTransformedMean;
}

/**
 * Returns whether an externally supplied asinh-scale mean is present.
 */
bool DyadicCovariate::hasAsinhTransformedMean() const
{
	return this->lhasAsinhTransformedMean;
}

/**
 * Returns the externally supplied asinh-scale mean.
 */
double DyadicCovariate::asinhTransformedMean() const
{
	return this->lasinhTransformedMean;
}

}

#endif /*DYADICCOVARIATE_H_*/

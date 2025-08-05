/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalAlterWeightedContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * TotalAlterWeightedContinuousEffect class.
 *****************************************************************************/

// modelled after DyadicCovariateAndNetworkBehaviorEffect.h and AverageAlterContinuousEffect.h

#ifndef TOTALINALTERWEIGHTEDCONTINUOUSEFFECT_H_
#define TOTALINALTERWEIGHTEDCONTINUOUSEFFECT_H_

#include "NetworkDependentContinuousEffect.h"

namespace siena
{


// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

//class Network;
//class BehaviorVariable;
class BehaviorLongitudinalData;
class ConstantDyadicCovariate;
class ChangingDyadicCovariate;
class DyadicCovariateValueIterator;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------


/**
 * Total alter effect defined as the total of an ego's neighbors (with
 * respect to a certain network).
 */
class TotalInAlterWeightedContinuousEffect : public NetworkDependentContinuousEffect
{
public:
    enum TransformType { NONE, LOG, ASINH };

    // Constructor with optional transform
    TotalInAlterWeightedContinuousEffect(const EffectInfo * pEffectInfo,
                                         TransformType transformType = NONE);

    virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double calculateChangeContribution(int actor);
	virtual double egoStatistic(int ego, double * currentValues);

protected:
	double dycoValue(int i, int j) const;
	bool missingDyCo(int i, int j) const;
	DyadicCovariateValueIterator rowValues(int i) const;
	DyadicCovariateValueIterator columnValues(int j) const;
	bool constantDyadicCovariate() const;
virtual void initializeStatisticCalculation();
virtual void cleanupStatisticCalculation();

private:
	// The constant covariate this effect depends on or 0, if the
	// effect depends on a changing covariate.
	ConstantDyadicCovariate * lpConstantDyadicCovariate;

	// The changing covariate this effect depends on or 0, if the
	// effect depends on a constant covariate.
	ChangingDyadicCovariate * lpChangingDyadicCovariate;

	BehaviorLongitudinalData * lpBehaviorData;

	// flag to control exclusion of missing values
	bool lexcludeMissings;

    // Selected transformation type
    TransformType lTransformType;
};

}

#endif /*TOTALINALTERWEIGHTEDCONTINUOUSEFFECT_H_*/

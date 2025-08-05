/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalDyadicCovariateInAltersEffect.h
 *
 * Description: This file contains the definition of the
 * TotalDyadicCovariateInAltersEffect class.
 *****************************************************************************/

#ifndef TOTALDYADICCOVARIATEINALTERSEFFECT_H_
#define TOTALDYADICCOVARIATEINALTERSEFFECT_H_

#include "NetworkDependentContinuousEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class BehaviorLongitudinalData;
class ConstantDyadicCovariate;
class ChangingDyadicCovariate;
class DyadicCovariateValueIterator;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------


/**
 * Total weighted in-alter effect defined as the product of an ego's behavior
 * and the sum of a dyadic covariate over all of an ego's in-alters (with
 * respect to a certain network).
 */
class TotalDyadicCovariateInAltersEffect : public NetworkDependentContinuousEffect
{
public:
    enum TransformType { NONE, LOG, ASINH };

    // Constructor with optional transformation type (default NONE)
    TotalDyadicCovariateInAltersEffect(const EffectInfo * pEffectInfo,
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

    // Selected transformation type for this effect
    TransformType lTransformType;

    // Mean of the TRANSFORMED dyadic covariate values (computed once per period)
    double lTransformedDycoMean;
};

}

#endif /*TOTALDYADICCOVARIATEINALTERSEFFECT_H_*/

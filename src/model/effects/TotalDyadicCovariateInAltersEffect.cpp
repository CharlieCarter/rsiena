/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalDyadicCovariateInAltersEffect.cpp
 *
 * Description: This file contains the implementation of the
 * TotalDyadicCovariateInAltersEffect class.
 *****************************************************************************/

#include <cmath>
#include <stdexcept>
#include "TotalDyadicCovariateInAltersEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

#include "model/variables/NetworkVariable.h"
#include "model/variables/ContinuousVariable.h"

#include "data/ConstantDyadicCovariate.h"
#include "data/ChangingDyadicCovariate.h"
#include "data/DyadicCovariateValueIterator.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/EpochSimulation.h"

// for unix-alike machines only
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
#include <R_ext/Print.h>
#endif

namespace siena
{

/**
 * Constructor.
 */
TotalDyadicCovariateInAltersEffect::TotalDyadicCovariateInAltersEffect(
    const EffectInfo * pEffectInfo,
    TransformType transformType) :
        NetworkDependentContinuousEffect(pEffectInfo),
        lTransformType(transformType),
        lTransformedDycoMean(0)
{
    this->lpConstantDyadicCovariate = 0;
    this->lpChangingDyadicCovariate = 0;
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void TotalDyadicCovariateInAltersEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkDependentContinuousEffect::initialize(pData, pState, period, pCache);
	string name1 = this->pEffectInfo()->interactionName1();
	string name2 = this->pEffectInfo()->interactionName2();

	this->lpConstantDyadicCovariate = pData->pConstantDyadicCovariate(name2);
	this->lpChangingDyadicCovariate = pData->pChangingDyadicCovariate(name2);
	    this->lpBehaviorData = pData->pBehaviorData(name1);
    this->lexcludeMissings = false;

    if (!this->lpConstantDyadicCovariate && !this->lpChangingDyadicCovariate)
    {
        throw logic_error(
            "Dyadic covariate variable '" + name2 + "' expected.");
    }

    // ------------------------------------------------------------------
    // Compute mean of TRANSFORMED dyadic covariate values for centering
    // ------------------------------------------------------------------
    if (this->lTransformType == NONE)
    {
        this->lTransformedDycoMean = 0; // not used when NONE
    }
    else
    {
        const Network * pNet = this->pNetwork();
        int nActors = pNet->n();
        double sumT = 0;
        long   count = 0;

        for (int i = 0; i < nActors; ++i)
        {
            if (this->lpConstantDyadicCovariate)
            {
                DyadicCovariateValueIterator iter = this->lpConstantDyadicCovariate->rowValues(i);
                for (; iter.valid(); iter.next())
                {
                    int j = iter.actor();
                    double raw = iter.value();
                    if (this->missingDyCo(i, j)) continue;
                    double tval = (this->lTransformType == LOG) ? std::log1p(raw) : std::asinh(raw);
                    sumT += tval;
                    ++count;
                }
            }
            else // changing dyadic covariate: average over ALL observations
            {
                int nObs = pData->observationCount();
                for (int obs = 0; obs < nObs; ++obs)
                {
                    DyadicCovariateValueIterator iter = this->lpChangingDyadicCovariate->rowValues(i, obs, /*excludeMissings*/false);
                    for (; iter.valid(); iter.next())
                    {
                        int j = iter.actor();
                        double raw = iter.value();
                        if (this->lpChangingDyadicCovariate->missing(i, j, obs)) continue;
                        double tval = (this->lTransformType == LOG) ? std::log1p(raw) : std::asinh(raw);
                        sumT += tval;
                        ++count;
                    }
                }
            }
        }
        // Compute grand mean of transformed values
        this->lTransformedDycoMean = (count > 0) ? sumT / static_cast<double>(count) : 0;
    }
}

/**
 * Returns the covariate value for the given pair of actors.
 */
double TotalDyadicCovariateInAltersEffect::dycoValue(int i, int j) const
{
	double value = 0;

    if (this->lpConstantDyadicCovariate)
    {
        value = this->lpConstantDyadicCovariate->value(i, j);
    }
    else
    {
        value = this->lpChangingDyadicCovariate->value(i, j, this->period());
    }

    // Perform centering only when no non-linear transformation is requested
    if (this->lTransformType == NONE)
    {
        const double mean = this->lpConstantDyadicCovariate
                              ? this->lpConstantDyadicCovariate->mean()
                              : this->lpChangingDyadicCovariate->mean();
        value -= mean;
    }

    return value;
}


/**
 * Returns if the covariate value for the given pair of actors is missing.
 */
bool TotalDyadicCovariateInAltersEffect::missingDyCo(int i, int j) const
{
	bool missing = false;

	if (this->lpConstantDyadicCovariate)
	{
		missing = this->lpConstantDyadicCovariate->missing(i, j);
	}
	else
	{
		missing = this->lpChangingDyadicCovariate->missing(i, j, this->period());
	}

	return missing;
}


/**
 * Returns if the associated covariate is a constant covariate or not
 */
bool TotalDyadicCovariateInAltersEffect::constantDyadicCovariate() const
{
	if (this->lpConstantDyadicCovariate)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/**
 * Returns an iterator over non-zero non-missing values of the given row
 * of the covariate.
 */
DyadicCovariateValueIterator
	TotalDyadicCovariateInAltersEffect::rowValues(int i) const
{
	if (this->lpConstantDyadicCovariate)
	{
		return this->lpConstantDyadicCovariate->rowValues(i);
	}
	else
	{
		return this->lpChangingDyadicCovariate->rowValues(i, this->period(),
			this->lexcludeMissings);
	}
}

/**
 * Returns an iterator over non-zero non-missing values of the given column
 * of the covariate.
 */
DyadicCovariateValueIterator
	TotalDyadicCovariateInAltersEffect::columnValues(int j) const
{
	if (this->lpConstantDyadicCovariate)
	{
		return this->lpConstantDyadicCovariate->columnValues(j);
	}
	else
	{
		return this->lpChangingDyadicCovariate->columnValues(j, this->period(),
			this->lexcludeMissings);
	}
}

/**
 * This method is called at the start of the calculation of the
 * evaluationStatistic, endowmentStatistic, and creationStatistic
 */
void TotalDyadicCovariateInAltersEffect::initializeStatisticCalculation()
{
		this->lexcludeMissings = true;
// Prevents having to check missingness in egoStatistic()
}

/**
 * This method is called at the end of the calculation of the
 * evaluationStatistic, endowmentStatistic, and creationStatistic
 */
void TotalDyadicCovariateInAltersEffect::cleanupStatisticCalculation()
{
		this->lexcludeMissings = false;
}


/**
 * Returns the sum of a certain actor's in-alters' dyadic covariate values,
 * and thus how much this effect contributes to the change in the continuous behavior.
 */
double TotalDyadicCovariateInAltersEffect::calculateChangeContribution(int actor)
{
        double contribution = 0;
    const Network * pNetwork = this->pNetwork();

    int indeg = pNetwork->inDegree(actor);
    if (indeg == 0)
    {
        return 0;
    }

    // (global transformed mean already pre-computed in initialize())



    // 1. Accumulate dyadic covariate values using helper
    for (IncidentTieIterator iter = pNetwork->inTies(actor);
         iter.valid();
         iter.next())
    {
        int j = iter.actor();
        contribution += this->dycoValue(j, actor);
    }

    // 2. Apply the requested non-linear transformation, if any
    if (this->lTransformType == LOG)
    {
        if (contribution < 0)
        {
            throw domain_error("TotalDyadicCovariateInAltersEffect LOG transform encountered negative value");
        }
        contribution = std::log1p(contribution);
    }
    else if (this->lTransformType == ASINH)
    {
        contribution = std::asinh(contribution);
    }

    // 3. Center the statistic (after transformation)
    if (this->lTransformType != NONE)
    {
        contribution -= this->lTransformedDycoMean;
    }

    return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the continuous behavior variable.
 */
double TotalDyadicCovariateInAltersEffect::egoStatistic(int ego, double * currentValues)
{
        double statistic = 0;
    const Network * pNetwork = this->pNetwork();

    int indeg = pNetwork->inDegree(ego);
    if (indeg == 0)
    {
        return 0;
    }



    // 1. Accumulate raw (un-centered) dyadic covariate values of incoming alters
    for (IncidentTieIterator iter = pNetwork->inTies(ego);
         iter.valid();
         iter.next())
    {
        int j = iter.actor();

        // When lexcludeMissings is false we must exclude missing values
        if (!this->lexcludeMissings)
        {
            if (this->missing(this->period(), j) ||
                this->missing(this->period() + 1, j) ||
                this->missingDyCo(j, ego))
            {
                continue;
            }
        }

        statistic += this->dycoValue(j, ego);
    }

    // 2. Apply transformation
    if (this->lTransformType == LOG)
    {
        if (statistic < 0)
        {
            throw domain_error("TotalDyadicCovariateInAltersEffect LOG transform encountered negative statistic");
        }
        statistic = std::log1p(statistic);
    }
    else if (this->lTransformType == ASINH)
    {
        statistic = std::asinh(statistic);
    }

    // 3. Center
    if (this->lTransformType != NONE)
    {
        statistic -= this->lTransformedDycoMean;
    }

    // 4. Multiply by ego's current value
    statistic *= currentValues[ego];

    return statistic;
}

}

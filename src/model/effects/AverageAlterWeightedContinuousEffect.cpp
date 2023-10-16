/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageAlterContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageAlterWeightedContinuousEffect class.
 *****************************************************************************/

#include <cmath>
#include "AverageAlterWeightedContinuousEffect.h"
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
AverageAlterWeightedContinuousEffect::AverageAlterWeightedContinuousEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentContinuousEffect(pEffectInfo)
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
void AverageAlterWeightedContinuousEffect::initialize(const Data * pData,
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
}

/**
 * Returns the covariate value for the given pair of actors.
 */
double AverageAlterWeightedContinuousEffect::dycoValue(int i, int j) const
{
	double value = 0;

	if (this->lpConstantDyadicCovariate)
	{
		value = this->lpConstantDyadicCovariate->value(i, j) -
			this->lpConstantDyadicCovariate->mean();
	}
	else
	{
		value = this->lpChangingDyadicCovariate->value(i, j, this->period()) -
			this->lpChangingDyadicCovariate->mean();
	}

	return value;
}


/**
 * Returns if the covariate value for the given pair of actors is missing.
 */
bool AverageAlterWeightedContinuousEffect::missingDyCo(int i, int j) const
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
bool AverageAlterWeightedContinuousEffect::constantDyadicCovariate() const
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
	AverageAlterWeightedContinuousEffect::rowValues(int i) const
{
	if (this->lpConstantDyadicCovariate)
	{
		return this->lpConstantDyadicCovariate->rowValues(i);
	}
	else
	{
		//	Rprintf("%d %d effect \n", i, this->lexcludeMissings);
		return this->lpChangingDyadicCovariate->rowValues(i, this->period(),
			this->lexcludeMissings);
	}
}

/**
 * Returns an iterator over non-zero non-missing values of the given column
 * of the covariate.
 */
DyadicCovariateValueIterator
	AverageAlterWeightedContinuousEffect::columnValues(int j) const
{
	if (this->lpConstantDyadicCovariate)
	{
		return this->lpConstantDyadicCovariate->columnValues(j);
	}
	else
	{
		//	Rprintf("%d effect \n", this->lexcludeMissings);
		return this->lpChangingDyadicCovariate->columnValues(j, this->period(),
			this->lexcludeMissings);
	}
}

/**
 * This method is called at the start of the calculation of the
 * evaluationStatistic, endowmentStatistic, and creationStatistic
 */
void AverageAlterWeightedContinuousEffect::initializeStatisticCalculation()
{
		this->lexcludeMissings = true;
// Prevents having to check missingness in egoStatistic()
}

/**
 * This method is called at the end of the calculation of the
 * evaluationStatistic, endowmentStatistic, and creationStatistic
 */
void AverageAlterWeightedContinuousEffect::cleanupStatisticCalculation()
{
		this->lexcludeMissings = false;
}


/**
 * Returns the average of a certain actor's alters, and thus how much
 * this effect contributes to the change in the continuous behavior.
 */
double AverageAlterWeightedContinuousEffect::calculateChangeContribution(int actor)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{
		double totalAlterValue = 0;
		double totalWeightValue = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
            int j = iter.actor();                // identifies alter
            double dycova = this->dycoValue(actor, j);
			double alterValue = this->centeredValue(iter.actor()) * dycova;  // for simstudy: value
			totalAlterValue += alterValue;
            totalWeightValue += (double) dycova;
		}

		contribution = totalAlterValue / totalWeightValue;
	}

	return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the continuous behavior variable.
 */
double AverageAlterWeightedContinuousEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int neighborCount = 0;
    double totalWeightValue = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		 iter.valid();
		 iter.next())
	{
		int j = iter.actor();

		if (!this->missing(this->period(), j) &&
			!this->missing(this->period() + 1, j) &&
            !this->missingDyCo(ego, j))
		{
            double dycova = this->dycoValue(ego, j);
			statistic += currentValues[j] * dycova;
			neighborCount++;
            totalWeightValue += (double) dycova;
		}
	}

    // may present issue if weight total is negative, though not anticipated
    // in current use-case
	if (neighborCount > 0)
	{
		statistic *= currentValues[ego] / totalWeightValue;
	}

	return statistic;
}

}

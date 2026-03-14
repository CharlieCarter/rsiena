/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageInAlterContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageInAlterContinuousEffect class.
 *****************************************************************************/

#include <cmath>
#include "AverageInAlterContinuousEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

#include "model/variables/NetworkVariable.h"
#include "model/variables/ContinuousVariable.h"
namespace siena
{

/**
 * Constructor.
 */
AverageInAlterContinuousEffect::AverageInAlterContinuousEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentContinuousEffect(pEffectInfo)
{
	this->ldivide = true;
	// Default behavior is to divide by indegree (average)
}

/**
 * Constructor.
 */
AverageInAlterContinuousEffect::AverageInAlterContinuousEffect(
	const EffectInfo * pEffectInfo, bool divide) :
		NetworkDependentContinuousEffect(pEffectInfo)
{
	this->ldivide = divide;
	// Indicates whether there will be division by the indegree of ego
}

/**
 * Precomputes the in-alter contribution for the given ego.
 * Called once per ego before calculateChangeContribution().
 */
void AverageInAlterContinuousEffect::preprocessEgo(int ego)
{
	ContinuousEffect::preprocessEgo(ego);
	this->lCachedContribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->inDegree(ego) > 0)
	{
		double totalInAlterValue = 0;

		for (IncidentTieIterator iter = pNetwork->inTies(ego);
			iter.valid();
			iter.next())
		{
			totalInAlterValue += this->centeredValue(iter.actor());
		}

		this->lCachedContribution = totalInAlterValue;
		if (this->ldivide)
		{
			this->lCachedContribution /= pNetwork->inDegree(ego);
		}
	}
}

/**
 * Returns the precomputed in-alter contribution (O(1)).
 */
double AverageInAlterContinuousEffect::calculateChangeContribution(int actor)
{
	return this->lCachedContribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the continuous behavior variable.
 */
double AverageInAlterContinuousEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int neighborCount = 0;

	for (IncidentTieIterator iter = pNetwork->inTies(ego);
		 iter.valid();
		 iter.next())
	{
		int j = iter.actor();

		if (!this->missing(this->period(), j) &&
			!this->missing(this->period() + 1, j))
		{
			statistic += currentValues[j];
			neighborCount++;
		}
	}

	if (neighborCount > 0)
	{
		if (this->ldivide)
		{
			statistic *= currentValues[ego] / neighborCount;
		}
		else
		{
			statistic *= currentValues[ego];
		}
	}

	return statistic;
}

}

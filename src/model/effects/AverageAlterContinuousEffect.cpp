/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageAlterContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageAlterContinuousEffect class.
 *****************************************************************************/

#include <cmath>
#include "AverageAlterContinuousEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

#include "model/variables/NetworkVariable.h"
#include "model/variables/ContinuousVariable.h"
namespace siena
{

/**
 * Constructor.
 */
AverageAlterContinuousEffect::AverageAlterContinuousEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentContinuousEffect(pEffectInfo)
{
}


/**
 * Precomputes the average alter value for the given ego.
 * Called once per ego before calculateChangeContribution().
 */
void AverageAlterContinuousEffect::preprocessEgo(int ego)
{
	ContinuousEffect::preprocessEgo(ego);
	this->lCachedContribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(ego) > 0)
	{
		double totalAlterValue = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(ego);
			iter.valid();
			iter.next())
		{
			totalAlterValue += this->centeredValue(iter.actor());
		}

		this->lCachedContribution = totalAlterValue / pNetwork->outDegree(ego);
	}
}


/**
 * Returns the precomputed average alter contribution (O(1)).
 */
double AverageAlterContinuousEffect::calculateChangeContribution(int actor)
{
	return this->lCachedContribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the continuous behavior variable.
 */
double AverageAlterContinuousEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int neighborCount = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
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
		statistic *= currentValues[ego] / neighborCount;
	}

	return statistic;
}

}

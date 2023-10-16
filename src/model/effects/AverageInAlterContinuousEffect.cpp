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
}


/**
 * Returns the average of a certain actor's alters, and thus how much
 * this effect contributes to the change in the continuous behavior.
 */
double AverageInAlterContinuousEffect::calculateChangeContribution(int actor)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->inDegree(actor) > 0)
	{
		double totalInAlterValue = 0;

		for (IncidentTieIterator iter = pNetwork->inTies(actor);
			iter.valid();
			iter.next())
		{
			double alterValue = this->centeredValue(iter.actor());  // for simstudy: value
			totalInAlterValue += alterValue;
		}

		contribution = totalInAlterValue / pNetwork->inDegree(actor);
	}

	return contribution;
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
		statistic *= currentValues[ego] / neighborCount;
	}

	return statistic;
}

}

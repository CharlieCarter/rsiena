/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageInAlterContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * AverageInAlterContinuousEffect class.
 *****************************************************************************/

#ifndef AVERAGEINALTERCONTINUOUSEFFECT_H_
#define AVERAGEINALTERCONTINUOUSEFFECT_H_

#include "NetworkDependentContinuousEffect.h"

namespace siena
{

/**
 * Average in-alter effect defined as the average of an ego's in-neighbors (with
 * respect to a certain network).
 */
class AverageInAlterContinuousEffect : public NetworkDependentContinuousEffect
{
public:
	AverageInAlterContinuousEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*AVERAGEINALTERCONTINUOUSEFFECT_H_*/

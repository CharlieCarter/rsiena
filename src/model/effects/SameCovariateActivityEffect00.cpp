/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateActivityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * SameCovariateActivityEffect class.
 *****************************************************************************/

#include <cmath>
#include "SameCovariateActivityEffect.h"
#include "utils/SqrtTable.h"
#include "utils/Utils.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"


namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
SameCovariateActivityEffect::SameCovariateActivityEffect(
		const EffectInfo * pEffectInfo, bool same, bool recip, bool straight) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lsame = same;
	this->lrecip = recip;
	this->lstraight = straight;
	this->lsqrt = (pEffectInfo->internalEffectParameter() == 2);
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Determines the condition combining everything
 */

bool SameCovariateActivityEffect::lcondition1(int theAlter, double theOwnValue) const
{
	return ((fabs(this->value(theAlter) - theOwnValue) < EPSILON) &&
			(!lrecip || (this->inTieExists(theAlter))));
}

bool SameCovariateActivityEffect::lcondition2(int theAlter, double theOwnValue) const
{
	return ((fabs(this->value(theAlter) - theOwnValue) >= EPSILON) &&
			(!lrecip || (this->inTieExists(theAlter))));
}

double SameCovariateActivityEffect::changeStat(double d) const
{
	if (this->lsqrt)
	{
		return(((d+1)*this->lsqrtTable->sqrt(d+1)) - (d * this->lsqrtTable->sqrt(d)));
	}
	else
	{
		return((2*d) + 1);
	}	
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
// For recip this is not very efficient, because each alter has the same contribution;
// could be made more efficient by preprocessEgo?
double SameCovariateActivityEffect::calculateContribution(int alter) const
{
	double myvalue = this->value(this->ego());
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if ((lsame) && 
		(lrecip || ((fabs(this->value(alter) - myvalue) < EPSILON) == lstraight)))
	{
		for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
		{
			// Get the receiver of the outgoing tie.
			int h = iter.actor();
			if (lcondition1(h, myvalue)) // (myvalue == this->value(h))
			{
				contribution++;
			}
		}
		if ((this->outTieExists(alter)) && lstraight)
		{
			contribution--;
		}
	}

	if ((!lsame) && 
		(lrecip || ((fabs(this->value(alter) - myvalue) >= EPSILON) == lstraight)))  
	{
		for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
		{
			// Get the receiver of the outgoing tie.
			int h = iter.actor();
			if (lcondition2(h, myvalue)) //  (myvalue != this->value(h))
			{
				contribution++;
			}
		}
		if ((this->outTieExists(alter)) && lstraight)
		{
			contribution--;
		}
	}
	return changeStat(contribution);
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double SameCovariateActivityEffect::tieStatistic(int alter)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (!(this->missing(alter) || this->missing(this->ego())))
	{
		double myvalue = this->value(this->ego());

		if (lsame)
		{
			if (lrecip || ((fabs(this->value(alter) - myvalue) < EPSILON) == this->lstraight))
			{
				for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
					iter.valid();
					iter.next())
				{
					// Get the receiver of the outgoing tie.
					int h = iter.actor();
					if ((!this->missing(h)) &&
							(lcondition1(h, myvalue))) // (myvalue == this->value(h))
					{
						contribution++;
					}
				}
			}
		}
		else
		{
		if (lrecip || ((fabs(this->value(alter) - myvalue) >= EPSILON) == this->lstraight))
			{
				for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
					iter.valid();
					iter.next())
				{
					// Get the receiver of the outgoing tie.
					int h = iter.actor();
					if ((!this->missing(h)) &&
							(lcondition2(h, myvalue)))	//  (myvalue != this->value(h))
					{
						contribution++;
					}
				}
			}
		}
	}
	if (this->lsqrt)
	{
		return this->lsqrtTable->sqrt(contribution);
	}
	else
	{
		return contribution;
	}
}

}

/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalInAlterWeightedContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the
 * TotalInAlterWeightedContinuousEffect class.
 *****************************************************************************/

 #include <cmath>
 #include "TotalInAlterWeightedContinuousEffect.h"
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
 TotalInAlterWeightedContinuousEffect::TotalInAlterWeightedContinuousEffect(
	 const EffectInfo * pEffectInfo,
	 TransformType transformType) :
		 NetworkDependentContinuousEffect(pEffectInfo),
		 lTransformType(transformType)
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
 void TotalInAlterWeightedContinuousEffect::initialize(const Data * pData,
	 State * pState,
	 int period,
	 Cache * pCache)
 {
	 // Rprintf("initializing total in alter weighted continuous effect");
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
 double TotalInAlterWeightedContinuousEffect::dycoValue(int i, int j) const
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
 bool TotalInAlterWeightedContinuousEffect::missingDyCo(int i, int j) const
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
 bool TotalInAlterWeightedContinuousEffect::constantDyadicCovariate() const
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
	 TotalInAlterWeightedContinuousEffect::rowValues(int i) const
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
	 TotalInAlterWeightedContinuousEffect::columnValues(int j) const
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
 void TotalInAlterWeightedContinuousEffect::initializeStatisticCalculation()
 {
		 this->lexcludeMissings = true;
 // Prevents having to check missingness in egoStatistic()
 }

 /**
  * This method is called at the end of the calculation of the
  * evaluationStatistic, endowmentStatistic, and creationStatistic
  */
 void TotalInAlterWeightedContinuousEffect::cleanupStatisticCalculation()
 {
		 this->lexcludeMissings = false;
 }


 /**
  * Returns the total of a certain actor's alters, and thus how much
  * this effect contributes to the change in the continuous behavior.
  */
 double TotalInAlterWeightedContinuousEffect::calculateChangeContribution(int actor)
 {
	 double contribution = 0;
	 const Network * pNetwork = this->pNetwork();

	 if (pNetwork->inDegree(actor) > 0)
	 {
		 for (IncidentTieIterator iter = pNetwork->inTies(actor);
			 iter.valid();
			 iter.next())
		 {
			 int j = iter.actor(); // identifies in-alter (sends tie to focal actor, despite inversion of 'i' and 'j')
			 double dycova = this->dycoValue(j, actor);
			 double inAlterValue = this->centeredValue(iter.actor());  // for simstudy: value
			 contribution += inAlterValue * dycova;
		 }
		 }

	 // Apply transformation if requested
	 switch (this->lTransformType)
	 {
		 case LOG:
			 if (contribution < 0)
			 {
				 throw domain_error("TotalInAlterWeightedContinuousEffect LOG transform encountered negative value");
			 }
			 contribution = std::log1p(contribution);
			 break;
		 case ASINH:
			 contribution = std::asinh(contribution);
			 break;
		 case NONE:
		 default:
			 break;
	 }

	 return contribution;
 }


 /**
  * Returns the statistic corresponding to the given ego with respect to the
  * given values of the continuous behavior variable.
  */
 double TotalInAlterWeightedContinuousEffect::egoStatistic(int ego, double * currentValues)
 {
	 double statistic = 0;
	 const Network * pNetwork = this->pNetwork();

	 for (IncidentTieIterator iter = pNetwork->inTies(ego);
		  iter.valid();
		  iter.next())
	 {
		 int j = iter.actor();

		 if (!(this->missing(this->period(), j) ||
			   this->missing(this->period() + 1, j) ||
			   this->missingDyCo(j, ego)))
		 {
			 double dycova = this->dycoValue(j, ego);
			 statistic += currentValues[j] * dycova;
		 }
	 }

	 if (pNetwork->inDegree(ego))
	 {
		 // Apply transformation if requested before multiplying by ego's value
		 switch (this->lTransformType)
		 {
			 case LOG:
				 if (statistic < 0)
				 {
					 throw domain_error("TotalInAlterWeightedContinuousEffect LOG transform encountered negative statistic");
				 }
				 statistic = std::log1p(statistic);
				 break;
			 case ASINH:
				 statistic = std::asinh(statistic);
				 break;
			 case NONE:
			 default:
				 break;
		 }

		 statistic *= currentValues[ego];
	 }

	 return statistic;
 }

 }

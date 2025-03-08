/********************************************************************/
/*                  SOFTWARE COPYRIGHT NOTIFICATION                 */
/*                             Cardinal                             */
/*                                                                  */
/*                  (c) 2021 UChicago Argonne, LLC                  */
/*                        ALL RIGHTS RESERVED                       */
/*                                                                  */
/*                 Prepared by UChicago Argonne, LLC                */
/*               Under Contract No. DE-AC02-06CH11357               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*             Prepared by Battelle Energy Alliance, LLC            */
/*               Under Contract No. DE-AC07-05ID14517               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*                 See LICENSE for full restrictions                */
/********************************************************************/

#ifdef ENABLE_OPENMC_COUPLING

#include "FETallyBase.h"

InputParameters
FETallyBase::validParams()
{
    auto params = TallyBase::validParams();
    params.addRequiredParam<std::vector<unsigned>>("orders",
                            "The orders for the expansions in each dimension. "
                            "These must be (x, y, z) for Legendre and "
                            "(rtheta, z) for Zernike.");
    params.addParam<std::string>("function_suffix","_function",
                            "The suffix to append to the score as the name of the"
                            "function object holding the Legendre expansion.");
    params.set<MultiMooseEnum>("output") = "UNRELAXED_TALLY";
    return params;
}

FETallyBase::FETallyBase(const InputParameters & parameters)
: TallyBase(parameters),
  _orders(getParam<std::vector<unsigned>>("orders")),
  _function_suffix(getParam<std::string>("function_suffix"))
{

    //overriding auxvariable names, don't want to create any
    if (isParamValid("name"))
    {
      mooseWarning(this->_name + " does not have any ElementalAuxVariables " 
                   "associated with it! "+this->_name+" creates functions, " 
                   "the names of which are controllable by \"function_suffix\". "
                   "Clearing \"name\" parameter...");
    }
    _tally_name.clear();

    // initializing number of functions and coefficients
    _functions.resize(_tally_score.size());
    
    /**
     * OpenMC spatial FETs only support the collision estimator
     */
    if (isParamValid("estimator"))
    {
        if (_estimator != openmc::TallyEstimator::COLLISION)
        paramError("estimator",
                    "Collision estimators are currently the only compatible "
                    "estimator type for Spatial Legendre expansion tallies!");
    }
    else
      _estimator = openmc::TallyEstimator::COLLISION;

    // initializing the first moments
    _first_moments.resize(_tally_score.size());

}

void 
FETallyBase::initializeTally()
{
    // Clear cached results.
  _local_sum_tally.clear();
  _local_sum_tally.resize(_tally_score.size(), 0.0);
  _local_mean_tally.clear();
  _local_mean_tally.resize(_tally_score.size(), 0.0);

  _current_tally.resize(_tally_score.size());
  _current_raw_tally.resize(_tally_score.size());
  _current_raw_tally_rel_error.resize(_tally_score.size());
  _current_raw_tally_std_dev.resize(_tally_score.size());
  _previous_tally.resize(_tally_score.size());

  auto [_filter_index, spatial_filters] = this->spatialFilters();

  std::vector<openmc::Filter *> filters;
  for (auto & filter : _ext_filters)
    filters.push_back(filter->getWrappedFilter());
  /**
   * We add the functional expansion filters last 
   * to minimize the number of cache misses during 
   * the OpenMC -> Cardinal transfer.
   */
  for (auto & filter : spatial_filters)
  {
    filters.push_back(filter);
  }

  // Create the tally, assign the required filters and apply the triggers.
  _local_tally_index = openmc::model::tallies.size();
  _local_tally = openmc::Tally::create();
  _local_tally->set_scores(_tally_score);
  _local_tally->estimator_ = _estimator;
  _local_tally->set_filters(filters);
  applyTriggersToLocalTally(_local_tally);
}

void
FETallyBase::resetTally()
{
  // Erase the tally.
  openmc::model::tallies.erase(openmc::model::tallies.begin() + _local_tally_index);

  for (int i = 0; i < getNumOrders(); i++)
  {
    openmc::model::tally_filters.erase(openmc::model::tally_filters.begin() + _filter_index + i);
  };
}

void
FETallyBase::computeSumAndMean()
{ 

  for (unsigned int score = 0; score < _tally_score.size(); ++score)
  {
    auto openmc_coeffs = _openmc_problem.tallySum(_local_tally, score);

    

    _first_moments.at(score) = openmc_coeffs(0);

    _local_sum_tally[score] = _first_moments.at(score);
    _local_mean_tally[score] = _first_moments.at(score) / getVolume();

    std::vector<Real> moose_coeffs;
    for (std::size_t i = 0; i < openmc_coeffs.size(); i++)
    {
      moose_coeffs.push_back(openmc_coeffs(i));
    }
    
    std::cout<<"First Moments:\t"<<openmc_coeffs(0)<<"\t"<<moose_coeffs[0]<<std::endl;

    _functions[score]->setCoefficients(moose_coeffs);
  }
}

Real
FETallyBase::storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                 unsigned int local_score,
                                 unsigned int global_score,
                                 std::vector<xt::xtensor<double, 1>> tally_vals,
                                 bool norm_by_src_rate)
{

    if (norm_by_src_rate)
    {
      Real norm_factor = _openmc_problem.tallyMultiplier(global_score);
      _first_moments.at(local_score) = this->normalizeCoefficients(local_score, norm_factor);
    }

    return this->_first_moments.at(local_score);
}

Real
FETallyBase::normalizeCoefficients( unsigned int local_score, Real factor)
{
  std::vector<Real> coeffs = _functions.at(local_score)->getCoefficients();
  for (std::size_t i = 0; i < coeffs.size(); i++)
  {
    coeffs.at(i) /= factor;
  }
  
  _functions.at(local_score)->setCoefficients(coeffs);

  return coeffs.at(0);
}

#endif

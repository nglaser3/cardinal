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
    params.registerBase("FETally");
    params.registerSystemAttributeName("FETally");
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

    //Verifying number of orders passed is correct
    mooseAssert(_orders.size() == getNumOrders(), 
              "Cardinal only supports 3-D, the length of \"orders\" "
              "for "+this->_name+" must be equal to "
              +std::to_string(getNumOrders())+".");

    // initializing number of functions and coefficients
    _functions.resize(_tally_score.size());
    _coefficients.resize(_tally_score.size());

    // getting number of terms for the coefficients
    size_t _size = 1.;
    for (int i = 0; i < _orders.size(); i++)
    {
        _size *= _orders.at(i);
    }
    
    // initializing functions and coefficients
    for (int index; index < _tally_score.size(); ++index)
    {
      _functions.at(index) = this->getFunctionSeries(_tally_score.at(index) + _function_suffix);
      _coefficients.at(index) = std::vector<Real>(_size);
    }
    
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

  auto [_filter_index, spatial_filters] = spatialFilters();

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
    _local_sum_tally[score] = _first_moments.at(score);
    _local_mean_tally[score] = _first_moments.at(score) / getVolume();
  }
}

Real
FETallyBase::storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                 unsigned int local_score,
                                 unsigned int global_score,
                                 std::vector<xt::xtensor<double, 1>> tally_vals,
                                 bool norm_by_src_rate)
{
  /**
   * TODO: local_score to the index of the function
   * to pass to the setCoefficients? 
   * DONE: Don't care about var_numbers
   * 
   */
    unsigned score_id = local_score;

    this->setCoefficients(tally_vals, local_score);

    if (norm_by_src_rate)
    {
      Real norm_factor = _openmc_problem.tallyMultiplier(global_score);
      this->normalizeCoefficients(score_id, norm_factor);
    }

    _first_moments.at(score_id) = _coefficients.at(score_id).at(0);
    
    return this->_first_moments.at(score_id);
}

void
FETallyBase::normalizeCoefficients( unsigned int score_id, Real factor)
{
  for (size_t i = 0; i < _coefficients.at(score_id).size(); i++)
  {
    _coefficients.at(score_id).at(i) *= factor;
  }
  _functions.at(score_id)->setCoefficients(_coefficients.at(score_id));
}

void
FETallyBase::save(unsigned score_id, size_t index, Real coefficient)
{
  
  if (index < _size && index >= 0)
  {
    _coefficients.at(score_id).at(index) = coefficient;
  }
  else mooseError("Hmm, something went wrong. The index trying to be saved is "
    + std::to_string(index) + " but the size of the coefficients vector is "
    + std::to_string(_size) +".");
}

#endif

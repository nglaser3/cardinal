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

registerMooseObject("CardinalApp", FETallyBase);

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

    verifyOrders();

    _functions.resize(_tally_score.size());
    _coefficients.resize(_tally_score.size());

    Real _size = 1.;
    for (int i = 0; i < _orders.size(); i++)
    {
        _size *= _orders.at(i);
    }
    

    for (int index; index < _tally_score.size(); ++index)
    {
      _functions.at(index) = this->getFunctionSeries(_tally_score.at(index) + _function_suffix);
      //initializing coefficients for functions, each with shape x * y * z (orders)
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
   * We add the three spatial legendre filters last 
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


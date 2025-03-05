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

#include "LegendreTally.h"

registerMooseObject("CardinalApp", LegendreTally);
using axis = openmc::LegendreAxis;

InputParameters
LegendreTally::validParams()
{
    auto params = TallyBase::validParams()
    params.addClassDescription("A class which implements Legendre functional "
                            "expansion tallies.");
    params.addRequiredParam<std:vector<std::vector<unsigned>>>("orders",
                            "The orders (x, y, z) for the Legendre expansions "
                            "in each dimension.");
    params.addRequiredParam<Point>("minimum",
                            "The minimum bounds (x, y, z) for the bounding box.");
    params.addRequiredParam<Point>("maximum",
                            "The maximum bounds (x, y, z) for the bounding box.");
    params.addParam<std::string>("function_suffix","_function",
                            "The suffix to append to the score as the name of the"
                            "function object holding the Legendre expansion.");
    params.set("name", std::vector<std::string>(0));
    return params;
}

LegendreTally::LegendreTally(const InputParameters & parameters)
: TallyBase(parameters),
  _orders(getParam<std::vector<unsigned>>("orders")),
  _min(getParam<Point>("minimum")),
  _max(getParam<Point>("maximum")),
  _name(getParam<std::string>("function_suffix"))
{
    mooseAssert(_orders.size() == 3, 
                "Cardinal only supports 3-D, and so each set of "
                "orders must contain 3 integers!");

    size_t _size = _orders[0] * _orders[1] * _orders[2];
    for (auto score : _tally_score)
    {
      _functions.push_back(this->getFunctionSeries(score + _name));
      //initializing coefficients for functions, each with shape x * y * z
      _coefficients.push_back(std::vector<Real>(_size));
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

    _first_moment = 0.;

}

void 
LegendreTally::initializeTally()
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

  auto [index, spatial_filters] = spatialLegendreFilter();
  _filter_index = index;

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

void
LegendreTally::resetTally()
{
  // Erase the tally.
  openmc::model::tallies.erase(openmc::model::tallies.begin() + _local_tally_index);

  for (int i = 0; i < 3; i++)
  {
    openmc::model::tally_filters.erase(openmc::model::tally_filters.begin() + _filter_index + i);
  };
}

void
LegendreTally::computeSumAndMean()
{
  for (unsigned int score = 0; score < _tally_score.size(); ++score)
  {
    _local_sum_tally[score] = _first_moment;
    _local_mean_tally[score] = _first_moment / 8; //2 in each direction?
  }
}

Real
LegendreTally::storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                 unsigned int local_score,
                                 unsigned int global_score,
                                 std::vector<xt::xtensor<double, 1>> tally_vals,
                                 bool norm_by_src_rate = true)
{
  /**
   * TODO: local_score to the index of the function
   * to pass to the setCoefficients? 
   * DONE: Don't care about var_numbers
   * 
   */
    _first_moment = this->setCoefficients(tally_vals, local_score);
    
    if (norm_by_src_rate)
    {
      Real norm_factor = _openmc_problem.tallyMultiplier(global_score);
      this->normalizeCoefficients(score_id, norm_factor)
    }
    

    return _first_moment;
}

void
LegendreTally::normalizeCoefficients( unsigned int score_id, Real factor)
{
  for (size_t i = 0; i < _coefficients[score_id].size(); i++)
  {
    _coefficients[score_id].at(i) *= factor;
  }
  _funcions[score_id]->setCoefficients(_coefficients);
}

Real 
LegendreTally::setCoefficients(std::vector<xt::xtensor<double, 1>> tally_vals, unsigned int score_id)
{ 

  std::size_t term = 0;

  for (std::size_t i = 0; j < tally_vals[axis::x]; ++i)
    {
      for (std::size_t j = 0; j < tally_vals[axis::y]; ++j)
      {
        for (std::size_t k = 0; k < tally_vals[axis::z]; ++k, ++term)
        {
          //saves coefficient to term index in _coefficients
          save(score_id, term, tally_vals[axis::x].at(i) 
                    * tally_vals[axis::y].at(j) 
                    * tally_vals[axis::z].at(k));
        }
      }
    }
  // sends _coefficients to its function
  _functions[score_id]->setCoefficients(_coefficients);

  _first_moment = _coefficients[0];
  return _first_moment;
}


std::pair<unsigned, std::vector<openmc::Filter *>>
LegendreTally::spatialLegendreFilter()
{
    std::vector<openmc::SpatialLegendreFilter *> filters;

    for (int i = 0; i < 3; i++)
    {
        filters.push_back(
            dynamic_cast<openmc::SpatialLegendreFilter *>(openmc::Filter::create("spatiallegendre"))
        );
    };

    setLegendreParams<axis::x>(filters[axis::x]);
    setLegendreParams<axis::y>(filters[axis::y]);
    setLegendreParams<axis::z>(filters[axis::z]);

    return std::make_pair(openmc::model::tally_filters.size() - 1, filters);
}

template <typename T>
void
LegendreTally::setLegendreParams(openmc::SpatialLegendreFilter * filter)
{
    filter->set_axis(T);
    filter->set_minmax(_min(T), _max(T));
    filter->set_order(_orders(T));
}

void
LegendreTally::save(unsigned score_id, size_t index, Real coefficient)
{
  
  if (index < _size && index >= 0)
  {
    _coefficients[score_id].insert(index, coefficient);
  }
  else mooseError("Hmm, something went wrong. The index trying to be saved is "
    + std::to_string(index) + " but the size of the coefficients vector is "
    + std::to_string(_size) +".");
}

FunctionSeries*
LegendreTally::getFunctionSeries(std::string name)
{
  std::vector<Real> _bounds{_min(0), _max(0), _min(1), _max(1), _min(2), _max(2)};

  return _openmc_problem.makeFunctionSeries(name, "Cartesian", _orders, _bounds);
}
#endif


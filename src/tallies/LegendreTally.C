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
    params.addRequiredParam<std:vector<unsigned int>>("orders",
                            "The orders (x, y, z) for the Legendre expansions "
                            "in each dimension.")
    params.addRequiredParam<Point>("minimum",
                            "The minimum bounds (x, y, z) for the bounding box.");
    params.addRequiredParam<Point>("maximum",
                            "The maximum bounds (x, y, z) for the bounding box.");
    return params;
}

LegendreTally::LegendreTally(const InputParameters & parameters)
: TallyBase(parameters),
  _orders(getParam<std::vector<unsigned int>>("orders")),
  _min(getParam<Point>("minimum")),
  _max(getParam<Point>("maximum"))
{
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
    return;
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

std::pair<unsigned, std::vector<openmc::Filter *>>
LegendreTally::spatialLegendreFilter()
{
    std::vector<openmc::SpatialLegendreFilter *> filters;
    std::vector<unsigned> ids;

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
#endif


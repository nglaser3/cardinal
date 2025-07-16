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

#include "SpatialLegendreTally.h"
#include <iostream>

registerMooseObject("CardinalApp", SpatialLegendreTally);

InputParameters
SpatialLegendreTally::validParams()
{
  auto params = FETallyBase::validParams();
  params.addRequiredParam<std::vector<unsigned>>("orders", "The expansion orders in x, y, and z.");
  params.addRequiredParam<Point>("lower_left", "The minimum bound for the x, y, and z diretions.");
  params.addRequiredParam<Point>("upper_right", "The maximum bound for the x, y, and z diretions.");
  return params;
}

SpatialLegendreTally::SpatialLegendreTally(const InputParameters & parameters)
 : FETallyBase(parameters),
   _orders(getParam<std::vector<unsigned int>>("orders")),
   _min(getParam<Point>("lower_left")),
   _max(getParam<Point>("upper_right"))
  {

    _tally_name.clear();

    std::vector<Real> bounds{_min(0), _max(0), _min(1), _max(1), _min(2), _max(2)};
    std::string func_name = _tally_score.at(0) + "_" + _function_suffix;
    _function = _openmc_problem.makeFunctionSeries(func_name, "Cartesian", _orders, bounds);
    std::cout<<"\n\n\n\n\033[92mCONSTRUCTOR\033[0m\n\n\n\n"<<std::endl;
  };

std::pair<unsigned int, std::vector<openmc::Filter *>>
SpatialLegendreTally::spatialFilter()
{
  auto first_id = openmc::model::tally_filters.size();
  std::vector<openmc::Filter*> filters;
  std::vector<openmc::LegendreAxis> axes{openmc::LegendreAxis::x,
                                         openmc::LegendreAxis::y,
                                         openmc::LegendreAxis::z};
  for (auto i = 0; i < axes.size(); i++)
  {
    auto filter = dynamic_cast<openmc::SpatialLegendreFilter *>(openmc::Filter::create("spatiallegendre"));
    filter->set_order(_orders[i]);
    filter->set_axis(axes[i]);
    filter->set_minmax(_min(i), _max(i));
    filters.push_back(filter);
  }
  std::cout<<"\n\n\n\n\033[92mSPTLFILTER\033[0m\n\n\n\n"<<std::endl;
  return std::make_pair(first_id, filters);
};
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

    std::vector<Real> bounds{_min(0), _max(0), _min(1), _max(1), _min(2), _max(2)};
    for (int i = 0; i < _functions.size(); i++)
    {
      for (int j = 0; j < _ext_filters.size(); j++)
      {
        std::string func_name = _tally_score.at(i) + "_" + _ext_filters.at(j)->name()
                                + "_" + _function_suffix;
        _functions.at(i).at(j) = _openmc_problem.makeFunctionSeries(func_name, "Cartesian",
                                                                    _orders, bounds);
      };
    };
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
  return std::make_pair(first_id, filters);
};

Real
SpatialLegendreTally::zerothMoment(xt::xtensor<double, 1> coefficients)
{
  return coefficients.at(0) * (_max(0) - _min(0))
                            * (_max(1) - _min(1))
                            * (_max(2) - _min(2));
};

Real
SpatialLegendreTally::firstMoment(xt::xtensor<double, 1> coefficients)
{
  return 2 * coefficients.at(1) / ( (_max(0) - _min(0))
                                   *(_max(1) - _min(1))
                                   *(_max(2) - _min(2)));
};
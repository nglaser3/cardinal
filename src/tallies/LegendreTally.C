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

#include "openmc/tallies/filter_sptl_legendre.h"

InputParameters
LegendreTally::validParams()
{
    auto params = FETallyBase::validParams();
    params.addClassDescription("A class which implements Legendre expansion tallies.");
    params.addRequiredParam<Point>("minimum",
                                   "The minimum bounds (x, y, z) for the bounding box.");
    params.addRequiredParam<Point>("maximum",
                                   "The maximum bounds (x, y, z) for the bounding box.");
    return params;
}

LegendreTally::LegendreTally(const InputParameters & parameters)
: FETallyBase(parameters),
  _min(getParam<Point>("minimum")),
  _max(getParam<Point>("maximum"))
{
}

std::pair<long unsigned int, std::vector<openmc::Filter *>> 
LegendreTally::spatialFilters()
{

    std::vector<openmc::Filter*> filters;

    for (int i = 0; i < 3; ++i)
    {
        auto filter = dynamic_cast<openmc::SpatialLegendreFilter *>(openmc::Filter::create("spatiallegendre"));
        filter->set_minmax(_min(i), _max(i));
        filter->set_order(_orders.at(i));
        filter->set_axis(static_cast<openmc::LegendreAxis>(i));
        filters.push_back(static_cast<openmc::Filter *>(filter));
    }

    return std::make_pair(openmc::model::tally_filters.size() - 3, filters);
}

void 
LegendreTally::setCoefficients(std::vector<xt::xtensor<double, 1>> tally_vals,
                               unsigned int score_id)
{
    std::size_t term = 0;

  for (std::size_t i = 0; i < tally_vals.at(0).size(); ++i)
    {
      for (std::size_t j = 0; j < tally_vals.at(1).size(); ++j)
      {
        for (std::size_t k = 0; k < tally_vals.at(2).size(); ++k, ++term)
        {
          //saves coefficient to term index in _coefficients
          save(score_id, term, tally_vals.at(0)(i) 
                             * tally_vals.at(1)(j) 
                             * tally_vals.at(2)(k) );
        }
      }
    }
  // sends _coefficients to its function
  _functions.at(score_id)->setCoefficients(_coefficients.at(score_id));
}

FunctionSeries*
LegendreTally::getFunctionSeries(std::string name)
{
    std::vector<Real> _bounds{_min(0), _max(0), _min(1), _max(1), _min(2), _max(2)};

    return _openmc_problem.makeFunctionSeries(name, "Cartesian", _orders, _bounds);
}

#endif
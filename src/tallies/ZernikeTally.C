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

#include "ZernikeTally.h"
#include <iostream>

registerMooseObject("CardinalApp", ZernikeTally);

InputParameters
ZernikeTally::validParams()
{
  auto params = FETallyBase::validParams();
  params.addRequiredParam<std::vector<unsigned>>("orders", "The expansion orders in r-theta and z");
  params.addRequiredParam<Real>("radius", "The radius of the circle to expand over");
  params.addRequiredParam<Point>("centroid", "The center point for the Zernike / axial Legendre expansion.");
  params.addRequiredParam<Real>("max_z", "The upper z bound for the axial Legendre expansion. The lower z bound is inferred from the centroid and the maximum.");
  return params;
}

ZernikeTally::ZernikeTally(const InputParameters & parameters)
 : FETallyBase(parameters),
   _orders(getParam<std::vector<unsigned int>>("orders")),
   _radius(getParam<Real>("radius")),
   _centroid(getParam<Point>("centroid")),
   _z_max(getParam<Real>("max_z"))
{
  // zmin, zmax, x0, y0, r
  std::vector<Real> bounds{ 2 * _centroid(2) - _z_max, _z_max, _centroid(0), _centroid(1), _radius + 1e-6}; //floating point shenannigans

  std::string func_name = _tally_score.at(0) + "_" + _function_suffix;
  _function = _openmc_problem.makeFunctionSeries(func_name, "CylindricalDuo", _orders, bounds);
  std::cout<<"\n\n\n\n\033[92mCONSTRUCTOR\033[0m\n\n\n\n"<<std::endl;
};

std::pair<unsigned int, std::vector<openmc::Filter *>>
ZernikeTally::spatialFilter()
{
  auto first_id = openmc::model::tally_filters.size();
  std::vector<openmc::Filter*> filters;
  
  // axial Legendre filter
  auto legendre_filter = dynamic_cast<openmc::SpatialLegendreFilter *>(openmc::Filter::create("spatiallegendre"));
  legendre_filter->set_order(_orders[0]);
  legendre_filter->set_axis(openmc::LegendreAxis::z);
  legendre_filter->set_minmax(2 * _centroid(2) - _z_max, _z_max);
  filters.push_back(legendre_filter);

  // Zernike Filter
  auto zernike_filter = dynamic_cast<openmc::ZernikeFilter *>(openmc::Filter::create("zernike"));
  zernike_filter->set_order(_orders[1]);
  zernike_filter->set_r(_radius);
  zernike_filter->set_x(_centroid(0));
  zernike_filter->set_y(_centroid(1));
  filters.push_back(zernike_filter);

  std::cout<<"\n\n\n\n\033[92mSPTLFILTER\033[0m\n\n\n\n"<<std::endl;
  return std::make_pair(first_id, filters);
}
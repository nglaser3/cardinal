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
#endif

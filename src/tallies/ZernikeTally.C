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

#include "ZernikeTally.h"

registerMooseObject("CardinalApp", ZernikeTally);

InputParameters
ZernikeTally::validParams()
{
    auto params = TallyBase::validParams();
    params.addClassDescription("A class which implements Zernike functional "
                            "expansion tallies.");
    params.addRequiredParam<MooseEnum>("normal_axis", getFilterAxisEnum(),
                            "The axis normal to the Zernike expansion. "
                            "Used to define the axis for the associated " 
                            "Legendre expansion tally.");
    params.addRequiredParam<unsigned int>("zernike_order",
                            "The order for the Zernike expansion.");
    params.addRequiredParam<unsigned int>("legendre_order",
                            "The order for the Legendre expansion.")
    params.addRequiredParam<Point>("centroid",
                            "The centroid (x0, y0) for the expansion.");
    params.addRequiredParam<Real>("radius",
                            "The radius for the Zernike expansion.");
    params.addRequiredParam<std::vector<Real>>("legendre_minmax",
                            "The lower and upper bounds (min, max) for "
                            "the Legendre expansion.");
    return params;
}

ZernikeTally::ZernikeTally(const InputParameters & parameters)
  : TallyBase(parameters),
    _axis(getParam<MooseEnum>("normal_axis").getEnum<LegendreAxis>()),
    _zernike_order(getParam<unsigned int>("zernike_order")),
    _legendre_order(getParam<unsigned int>("legendre_order")),
    _centroid(getParam<std::vector<Real>>("centroid")),
    _radius(getParam<Real>("radius")),
    _minmax(getParam<Real>("legendre_minmax"))
{
    /**
     * OpenMC Zernike filters currently only support expansion in the xy plane.
     * The user provided axis must be z, otherwise ZernikeTally throws an error.
     */
    if (_axis != LegendreAxis::z)
    {
        paramError("normal_axis",
            "OpenMC Zernike filters currently only support expansion in the xy plane, "
            "and so the normal_axis must be set to z.");
    }
    
    /**
     * OpenMC spatial FETs only support the collision estimator
     */
    if (isParamValid("estimator"))
    {
        if (_estimator != openmc::TallyEstimator::COLLISION)
        paramError("estimator",
                    "Collision estimators are currently the only compatible "
                    "estimator type for Zernike expansion tallies!");
    }

    return;
}

#endif


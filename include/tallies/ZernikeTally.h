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

#pragma once

#include "TallyBase.h"
#include "OpenMCCellAverageProblem.h"
#include "openmc/tallies/filter_sptl_legendre.h"
#include "openmc/tallies/filter_zernike.h"
#include "CardinalEnums.h"


class ZernikeTally : public TallyBase
{
    public:
        static InputParameters validParams();

        ZernikeTally(const InputParameters & parameters)

        virtual std::pair<unsigned int, openmc::Filter *> spatialFilter() override;

        virtual void resetTally() override;

    protected:
        LegendreAxis _axis;
        unsigned int _zernike_order;
        unsigned int _legendre_order;
        Point _centroid;
        Real _radius;
        Real _minmax;
}

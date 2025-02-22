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

class LegendreTally : public TallyBase
{
    public:
        static InputParameters validParams();

        LegendreTally(const InputParameters & params)

        virtual std::pair<unsigned int, openmc::Filter *> spatialFilter() override;

        virtual void resetTally() override;

    protected:
        std::vector<unsigned int> _orders;
        Point _min;
        Point _max;
}
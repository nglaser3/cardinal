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

class LegendreTally : public TallyBase
{
    public:
        static InputParameters validParams();

        LegendreTally(const InputParameters & params)

        virtual void initializeTally() override;

        virtual void TallyBase::resetTally() override;

        /**
         * spatialFilter is marked as pure virtual, so must be overriden
         * return 0 and null ptr just to satisfy override
         */
        virtual std::pair<unsigned int, openmc::Filter *> spatialFilter() override
        {return std::make_pair(0, nullptr)};

        virtual void resetTally() override;

    protected:
        std::pair<unsigned, std::vector<openmc::Filter *>> spatialLegendreFilter();

        template <typename T> 
        void setLegendreParams(openmc::SpatialLegendreFilter * filter);

        std::vector<unsigned int> _orders;
        Point _min;
        Point _max;

        std::vector<openmc::SpatialLegendreFilter *> _legendre_filter;
        std::vector<unsigned> _filter_ids;
}
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
#include "FunctionSeries.h"

class LegendreTally : public TallyBase
{
    public:
        static InputParameters validParams();

        LegendreTally(const InputParameters & params)

        virtual void initializeTally() override;

        virtual void resetTally() override;

        virtual void computeSumAndMean() override;



        /**
         * spatialFilter is marked as pure virtual, so must be overriden
         * return 0 and null ptr just to satisfy override
         */
        virtual std::pair<unsigned int, openmc::Filter *> spatialFilter() override
        {return std::make_pair(0, nullptr)};

        virtual void resetTally() override;

    protected:
        virtual Real storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                       unsigned int local_score,
                                       unsigned int global_score,
                                       std::vector<xt::xtensor<double, 1>> tally_vals,
                                       bool norm_by_src_rate = true) override;

        std::pair<unsigned, std::vector<openmc::Filter *>> spatialLegendreFilter();

        static template <typename T> 
        void setLegendreParams(openmc::SpatialLegendreFilter * filter);

        std:vector<FunctionSeries*> _functions;
        std::vector<unsigned int> _orders;
        Point _min;
        Point _max;

        std::vector<Real> _coefficients;


}
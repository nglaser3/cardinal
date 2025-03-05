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

        LegendreTally(const InputParameters & params);

        virtual void initializeTally() override;

        virtual void resetTally() override;

        virtual void computeSumAndMean() override;

        /**
         * spatialFilter is marked as pure virtual, so must be overriden
         * return 0 and null ptr just to satisfy override
         */
        virtual std::pair<unsigned int, openmc::Filter *> spatialFilter() override
        {return std::make_pair(0, nullptr);};

    protected:
        virtual Real storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                       unsigned int local_score,
                                       unsigned int global_score,
                                       std::vector<xt::xtensor<double, 1>> tally_vals,
                                       bool norm_by_src_rate) override;


        void normalizeCoefficients(unsigned int score_id, Real factor);
        
        void setCoefficients(std::vector<xt::xtensor<double, 1>> tally_vals, unsigned int score_id);

        std::pair<unsigned int, std::vector<openmc::SpatialLegendreFilter *>> spatialLegendreFilter();

        void setLegendreParams(openmc::LegendreAxis axis,
                               openmc::SpatialLegendreFilter * filter);

        void save(unsigned score_id, size_t index, Real coefficient);

        FunctionSeries* getFunctionSeries(std::string name);
        
        std::vector<unsigned int> _orders;
        Point _min;
        Point _max;

        size_t _size;

        std::vector<FunctionSeries*> _functions;
        std::vector<std::vector<Real>> _coefficients;

        std::vector<Real> _first_moments;

       std::string _function_suffix;
};
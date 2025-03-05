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

class FETallyBase : public TallyBase
{
    public:

        static InputParameters validParams();

        FETallyBase(const InputOParameters & parameters);

        virtual void initializeTally() override;

        virtual void resetTally() override;

        virtual std::pair<unsigned int, openmc::Filter *> spatialFilter() override
        {return std::make_pair(0, nullptr);};

        virtual std::pair<unsigned int, std::vector<openmc::Filter *>> spatialFilters() = 0;  

    protected:

        virtual Real storeResultsInner(const std::vector<unsigned int> & var_numbers,
            unsigned int local_score,
            unsigned int global_score,
            std::vector<xt::xtensor<double, 1>> tally_vals,
            bool norm_by_src_rate) override;

        void normalizeCoefficients(unsigned int score_id, Real factor);

        virtual void setCoefficients(std::vector<xt::xtensor<double, 1>> tally_vals,
                                     unsigned int score_id) = 0;

        void save(unsigned score_id, size_t index, Real coefficient);

        virtual FunctionSeries* getFunctionSeries(std::string name) = 0;

        static openmc::SpatialLegendreFilter* makeLegendreFilter(Real min,
                                                Real max, unsigned order,
                                                openmc::LegendreAxis axis);

        virtual void verifyOrders() = 0;
        std::vector<unsigned int> _orders;

        size_t _size;

        std::vector<FunctionSeries*> _functions;
        std::vector<std::vector<Real>> _coefficients;

        std::string _function_suffix;

};


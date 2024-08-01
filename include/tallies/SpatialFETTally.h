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
#include "CardinalEnums.h"

#include "openmc/tallies/filter_zernike.h"
#include "openmc/tallies/filter_sptl_legendre.h"

class SpatialFETTally : public TallyBase
{
public:
    static InputParameters validParams(); //done?

    SpatialFETTally(const InputParameters & parameters); //done?

    virtual void initializeTally() override; //done?

    virtual void resetTally() override; //copy pasted

    //TODO
    virtual Real storeResults(const std::vector<unsigned int> & var_numbers,
                            unsigned int local_score,
                            unsigned int global_score) override;
    
    //TODO
    virtual void storeExternalResults(const std::vector<unsigned int> & ext_var_numbers,
                                    unsigned int local_score,
                                    unsigned int global_score,
                                    const std::string & output_type) override;

protected:
    //done?
    template <typename T> 
    T* makeFETFilter(int32_t id);

    const unsigned int _order;

    const MooseEnum _func_type;
    const Real _radius;
    const std::vector<Real> _centroid;
    const openmc::LegendreAxis _leg_axis;

    int32_t _filter_index;
    int32_t _local_tally_index;

    
};

 

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

#include "FETallyBase.h"
#include "openmc/tallies/filter_sptl_legendre.h"

class SpatialLegendreTally : public FETallyBase
{
  public:

    static InputParameters validParams();

    SpatialLegendreTally(const InputParameters & parameters);

    virtual std::pair<unsigned int, std::vector<openmc::Filter *>> spatialFilter() override;

  protected:

    virtual int getNumBins() override {return (_orders[0] + 1) * (_orders[1] + 1) * (_orders[2] + 1);}

    virtual int getNumSptlFilters() override {return 3;}

    std::vector<unsigned> _orders;
    const Point & _min;
    const Point & _max;
};

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
#include "MutableCoefficientsInterface.h"

class Zernike
{
private:
    Real _x,_y;
    Real _order;
    Real
public:
    SpatialFETally(/* args */);
    ~SpatialFETally();
};


class SpatialFETally : public TallyBase
{
public:
    static InputParameters validParams();
    SpatialFETally(InputParameters & parameters);

    virtual void initializeTally() override; 

    virtual void resetTally() override;

    virtual Real storeResults(const std::vector<unsigned int> & var_numbers,
                            unsigned int local_score,
                            unsigned int global_score) override;

    virtual void storeExternalResults(const std::vector<unsigned int> & ext_var_numbers,
                                    unsigned int local_score,
                                    unsigned int global_score,
                                    const std::string & output_type) override;
protected:
    std::vector<Real> _current_coeffs;
    std::vector<

};



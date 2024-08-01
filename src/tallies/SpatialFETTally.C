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

//#ifdef ENABLE_OPENMC_COUPLING 
#include "SpatialFETTally.h"

registerMooseObject("CardinalApp",SpatialFETTally);

InputParameters
SpatialFETTally::validParams()
{
    auto params = TallyBase::validParams();
    params.addClassDescription(
    "This class implements spatial functional expansion tallies (FET).
    Importantly, this class does NOT expand angular components of multigroup cross sections.");
    params.addRequiredParam<unsigned int>("order","Order to expand tally eigenfunction to");
    params.addRequiredParam<MooseEnum>("function_type",getFETTypeEnum(),
    "Eigenfunction family to use for FET");
    params.addRequirdParam<Real>("radius",
    "Radius of the domain in which the FET is expanded over. For zernike this is simply
    the radius. For Legendre this is half of the domain, i.e. if the domain is (0,1) 
    the radius would be 0.5.")
    params.addRequiredParam<std::vector<Real>>("centroid",
    "Centroid of the expansion. For Zernike pass {x,y} and 
    for spatial legendre pass one value in a vector: {x} or {y} or {z}.")
    params.addParam<openmc::LegendreAxis>("legendre_axis",
    "Axis to expand the legendre FET along.")
    return params;
}

// 

SpatialFETTally::SpatialFETTally(const InputParameters & parameters)
    : TallyBase(parameters),
      _order(getParam<unsigned int>("order")),
      _func_type(getParam<MooseEnum>("function_type")),
      _radius(getParam<Real>("radius")),
      _centroids(getParam<std::vector<std::vector<Real>>>("centroids"))
{
    if (isParamValid("legendre_axis"))
    {
        _leg_axis(getParam<openmc::LegendreAxis>("legendre_axis"))
    }
}

void
SpatialFETTally::initializeTally()
{
    // Clear cached results.
    _local_sum_tally.clear();
    _local_sum_tally.resize(_tally_score.size(), 0.0);
    _local_mean_tally.clear();
    _local_mean_tally.resize(_tally_score.size(), 0.0);

    _current_tally.resize(_tally_score.size());
    _current_raw_tally.resize(_tally_score.size());
    _current_raw_tally_std_dev.resize(_tally_score.size());
    _previous_tally.resize(_tally_score.size());
    
    //initialize Filter vector and filter id
    std::vector<openmc::Filter *> filters{};
   
    int32_t _filter_index = openmc::model::tally_filters.size();
        
    switch (_func_type)
    {
    case tally::zernike:
        filters.pushback(makeFETFilter<openmc::ZernikeFilter>(_filter_index));
        break;
    case tally::zernike_radial:
        filters.pushback(makeFETFilter<openmc::ZernikeRadialFilter>(_filter_index));
    case tally::spatial_legendre:
        filters.pushback(makeFETFilter<openmc::SpatialLegendreFilter>(_filter_index))
    default:
        break;
    }

    _local_tally_index = openmc::model::tallies.size();
    _local_tally = openmc::Tally::create(_local_tally_index);
    _local_tally->set_scores(_tally_score);
    _local_tally->estimator_ = _estimator;
    _local_tally->set_filters(filters);
    applyTriggersToLocalTally(_local_tally);
    
}

void
SpatialFETTally::resetTally()
{
    for (int i = 0; i < _local_tally_index.size(); i++)
        {
            // Erase the tally.
            openmc::model::tallies.erase(openmc::model::tallies.begin() + _local_tally_index[i]);

            // Erase the filter(s).
            openmc::model::tally_filters.erase(openmc::model::tally_filters.begin() + _filter_index[i]);
        }
  
}

Real //what is this doing???
SpatialFETTally::storeResults(const std::vector<unsigned int> & var_numbers,
                            unsigned int local_score,
                            unsigned int global_score)
{

}

void //and this?
SpatialFETTally::storeExternalResults(const std::vector<unsigned int> & ext_var_numbers,
                                unsigned int local_score,
                                unsigned int global_score,
                                const std::string & output_type)
{
    
}
openmc::Filter*
template <typename T> // filter.cpp line 97 or filter.h line 165?
T* SpatialFETTally::makeFETFilter(int32_t id)
{   
    //
    T* _filter = openmc::Filter::create<T>(id);

    switch (_func_type)
    {
    case tally::zernike:
        //_filter = dynamic_cast<T *>(openmc::Filter::create("zernike"));
        _filter->set_x(_centroid[0]);
        _filter->set_y(_centroid[1]);
        _filter->set_r(_radius);
        break;
    case tally::zernike_radial:
        //_filter = dynamic_cast<T *>(openmc::Filter::create("zernikeradial"));
        _filter->set_x(_centroid[0]);
        _filter->set_y(_centroid[1]);
        _filter->set_r(_radius);
        break;
    case tally::spatial_legendre:
        //_filter = dynamic_cast<T *>(openmc::Filter::create("spatiallegendre"));
        _filter->set_minmax(_centroid[0]-_radius,_centroid[0]+_radius);
        break;
    default:
        break;
    }
    _filter->set_order(_order);
    return _filter;
}

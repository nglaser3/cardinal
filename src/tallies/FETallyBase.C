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

#ifdef ENABLE_OPENMC_COUPLING

#include "FETallyBase.h"

InputParameters
FETallyBase::validParams()
{
    auto params = TallyBase::validParams();
    params.addParam<std::string>("function_suffix","_function",
                            "The suffix to append to the score(s) as the name of the"
                            "function object(s) holding the functional expansion(s).");
    params.set<MultiMooseEnum>("output") = "UNRELAXED_TALLY";
    return params;
}

FETallyBase::FETallyBase(const InputParameters & parameters)
: TallyBase(parameters),
  _function_suffix(getParam<std::string>("function_suffix"))
{

    // Overriding auxvariable names, don't want to create any
    if (isParamValid("name"))
    {
      mooseWarning(this->_name + " does not have any ElementalAuxVariables " 
                   "associated with it! "+this->_name+" creates functions, " 
                   "the names of which are controllable by \"function_suffix\". "
                   "Clearing \"name\" parameter...");
    }
    _tally_name.clear();
    
    // OpenMC spatial FETs only support the collision estimator
    if (isParamValid("estimator"))
    {
        if (_estimator != openmc::TallyEstimator::COLLISION)
        paramError("estimator",
                    "Collision estimators are currently the only compatible "
                    "estimator type for Spatial Legendre expansion tallies!");
    }
    else
      _estimator = openmc::TallyEstimator::COLLISION;

    // initializing vector size for functions
    _functions.resize(_tally_score.size());
    for (int i = 0; i < _functions.size(); i++)
      _functions[i].resize(_ext_filters.size());
}

void
FETallyBase::resetTally()
{
  // Erase the tally.
  openmc::model::tallies.erase(openmc::model::tallies.begin() + _local_tally_index);

  for (int i = 0; i < getNumSptlFilters(); i++)
  {
    openmc::model::tally_filters.erase(openmc::model::tally_filters.begin() + _filter_index + i);
  };
}

void
FETallyBase::computeSumAndMean()
{ 
  auto num_bins = getNumBins();
  for (unsigned int score = 0; score < _tally_score.size(); ++score)
  {
    _local_sum_tally[score] = 0.0;
    _local_mean_tally[score] = 0.0;
    for (unsigned int ext = 0; ext < _num_ext_filter_bins; ++ext)
    {
      if (!_ext_bins_to_skip[ext])
      {
        auto coefficients = xt::view(_local_tally->results_,
                                 xt::range(ext * num_bins, ext * num_bins + num_bins),
                                 score,
                                 static_cast<int>(openmc::TallyResult::SUM));
        _local_sum_tally[score] += zerothMoment(coefficients);
        _local_mean_tally[score] += firstMoment(coefficients);
      };
    };
  };
}

void
FETallyBase::relaxAndNormalizeTally(unsigned int local_score,
                                      const Real & alpha,
                                            const Real & norm)
{
  // TODO
}

Real
FETallyBase::storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                 unsigned int local_score,
                                 unsigned int global_score,
                                 std::vector<xt::xtensor<double, 1>> tally_vals,
                                 bool norm_by_src_rate)
{
  Real total = 0.0;
  auto num_bins = getNumBins();
  for (unsigned int ext = 0; ext < _num_ext_filter_bins; ++ext)
    {
      auto xt_coeffs = xt::view(tally_vals[local_score],
                                xt::range(ext * num_bins, ext * num_bins + num_bins));
      total += zerothMoment(xt_coeffs);
      std::vector<Real> coeffs(xt_coeffs.begin(), xt_coeffs.end());
      _functions.at(local_score).at(ext)->setCoefficients(coeffs);
    };
    return total;
}
#endif

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
    //_tally_name.clear();
    
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
  for (unsigned int score = 0; score < _tally_score.size(); ++score)
  {
    std::cout<<"\n\n\n\n\033[92mCOMPUTESUMANDMEAN\033[0m\n\n\n\n"<<std::endl;
    auto xt_coeffs = xt::view(_local_tally->results_,
                                 xt::all(),
                                 score,
                                 static_cast<int>(openmc::TallyResult::SUM));
    std::vector<Real> coeffs(xt_coeffs.begin(), xt_coeffs.end());
    _function->setCoefficients(coeffs);
    auto [integral, volume] = computeIntegral(_function);
    _local_sum_tally[score] = integral;
    _local_mean_tally[score] = integral / volume;
  }
}

Real
FETallyBase::storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                 unsigned int local_score,
                                 unsigned int global_score,
                                 std::vector<xt::xtensor<double, 1>> tally_vals,
                                 bool norm_by_src_rate)
{
  std::cout<<"\n\n\n\n\033[92mSTORESRESULTSINNER\033[0m\n\n\n\n"<<std::endl;
  Real total = 0.0;
  auto xt_coeffs = xt::view(_local_tally->results_,
                            xt::all(),
                            local_score,
                            static_cast<int>(openmc::TallyResult::SUM));
  std::vector<Real> coeffs(xt_coeffs.begin(), xt_coeffs.end());
  _function->setCoefficients(coeffs);
  return total;
}

std::pair<Real, Real>
FETallyBase::computeIntegral(FunctionSeries* function)
{
  // Setup quadrature rule for integrating
  auto _fe = FEBase::build(3, FEType(1, FEFamily::LAGRANGE));
  auto _qr = QBase::build(QuadratureType::QGAUSS, 3, Order::TWENTIETH);
  _fe->attach_quadrature_rule(_qr.get());

  // grab the points and weights for integration
  const std::vector<Point> & points = _fe->get_xyz();
  const std::vector<Real> & weights = _fe->get_JxW();

  Real integral_value = 0.0;
  Real volume = 0;
  for(const auto * elem : _openmc_problem.getMooseMesh().getMesh().active_local_element_ptr_range())
  {
    /**
     * need to reinit the element, otherwise the points and weights are meaningless.
     * libMesh::FEBase::reinit(elem) operates on the points and weights vector, 
     * filling them with the correct points and weights for the current element
     */
    _fe->reinit(elem);
    for (unsigned int i = 0; i < points.size(); i++)
    {
      integral_value += weights[i] * function->evaluateValue(0.0, points[i]);
    }
    volume += elem->volume();
  }

  // bring everything back together after done from mpi (?) 
  const auto & comm = _openmc_problem.comm();
  comm.sum(integral_value);
  comm.sum(volume);

  return std::make_pair(integral_value, volume);
}
#endif

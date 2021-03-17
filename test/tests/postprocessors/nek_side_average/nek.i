[Problem]
  type = NekRSProblem
[]

[Mesh]
  type = NekRSMesh
  boundary = '1 2 3 4 5 6 7 8'
[]

[Executioner]
  type = Transient

  [TimeStepper]
    type = NekTimeStepper
  []
[]

[Outputs]
  [out]
    type = CSV
    hide = 'flux_integral'
    execute_on = 'final'
  []
[]

[Postprocessors]
  [temp_avg1]
    type = NekSideAverage
    field = temperature
    boundary = '1'
  []
  [temp_avg2]
    type = NekSideAverage
    field = temperature
    boundary = '2'
  []
  [temp_avg3]
    type = NekSideAverage
    field = temperature
    boundary = '3'
  []
  [temp_avg4]
    type = NekSideAverage
    field = temperature
    boundary = '4'
  []
  [temp_avg5]
    type = NekSideAverage
    field = temperature
    boundary = '5'
  []
  [temp_avg6]
    type = NekSideAverage
    field = temperature
    boundary = '6'
  []
  [temp_avg7]
    type = NekSideAverage
    field = temperature
    boundary = '7'
  []
  [temp_avg8]
    type = NekSideAverage
    field = temperature
    boundary = '8'
  []
  [pressure_avg1]
    type = NekSideAverage
    field = pressure
    boundary = '1'
  []
  [pressure_avg2]
    type = NekSideAverage
    field = pressure
    boundary = '2'
  []
  [pressure_avg3]
    type = NekSideAverage
    field = pressure
    boundary = '3'
  []
  [pressure_avg4]
    type = NekSideAverage
    field = pressure
    boundary = '4'
  []
  [pressure_avg5]
    type = NekSideAverage
    field = pressure
    boundary = '5'
  []
  [pressure_avg6]
    type = NekSideAverage
    field = pressure
    boundary = '6'
  []
  [pressure_avg7]
    type = NekSideAverage
    field = pressure
    boundary = '7'
  []
  [pressure_avg8]
    type = NekSideAverage
    field = pressure
    boundary = '8'
  []
  [flux_integral]
    type = Receiver
  []
[]

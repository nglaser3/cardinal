[Mesh]
  [sphere]
    type = FileMeshGenerator
    file = ../meshes/pincell.e
  []
[]

[Problem]
    type = OpenMCCellAverageProblem
    verbose = true
    power = 500.0
    initial_properties = xml
    cell_level = 0
    check_tally_sum = false

    [Tallies]
        [Zernike]
            type = ZernikeTally
            orders = '10 10'
            centroid = '0 0 5'
            radius = .4125
            half_range = 5.
            blocks = 1
        []
    []
[]

[Executioner]
    type = Transient
[]
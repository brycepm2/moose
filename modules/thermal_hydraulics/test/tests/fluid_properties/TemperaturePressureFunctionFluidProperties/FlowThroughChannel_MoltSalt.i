# Some inputs for the model
# Thermal Hydraulics model uses units kg, m, s, K
# ICs
T_in = 1100. # K
press = 1e5 # Pa

# Flow
m_dot_in = 1 #kg/s
r_in = 0.01 # this is id
A_in = '${fparse pi*r_in^2}'
D_h_in = '${fparse 2*r_in}'
v_in = '${fparse m_dot_in/(Bd*A_in)}'

# density
Ad = -0.552
Bd =  1903.7
beta_in = 0.547e-9

# specific heat
C0 = 989.6
C1 = 0.1046
C2 = 430

# thermal cond
K0 = 0.5047
K1 = 0.0001

# visco
M0 = 1.4965e-2
M1 = -2.91e-5
M2 = 1.784e-8

# energy/constants
T_0_in = 273.15      # K
P_0_in = 10e5
e_0_in = 0

# GLOBAL PARAMETERS
# used for params that are global to the sim
# NOTE: ~rdg_slope_reconstruction = full~ should be in almost all input files
#   it allows full slope reconstruction in the rDG scheme
# It can be convenient to put initial conditions here too
# Closures can also be added here
# Since helium is used everywhere also, that can be put in global params
[GlobalParams]
  initial_p = ${press}
  initial_vel = ${v_in}
  initial_T = ${T_in}
  gravity_vector = '0 0 0'

  rdg_slope_reconstruction = full
  scaling_factor_1phase = '1 1e-2 1e-4'
  closures = thm_closures
  fp = fp 
[]

[Functions]
  [rho]
    type = ParsedFunction
    expression = '(Ad*(x - T_0) + Bd)*(1 + beta*(y - P_0))'
    symbol_names ='Ad Bd T_0 beta P_0'
    symbol_values = '${Ad} ${Bd} ${T_0_in} ${beta_in} ${P_0_in}'
  []
  [k]
    type = ParsedFunction
    expression = 'K0 - K1*(x - T_0)'
    symbol_names ='K0 K1 T_0'
    symbol_values = '${K0} ${K1} ${T_0_in}'
  []
  [cp]
    type = ParsedFunction
    expression = 'C0 + C1*(x - T_0 - C2)'
    symbol_names ='C0 C1 C2 T_0'
    symbol_values ='${C0} ${C1} ${C2} ${T_0_in}'
  []
  [mu]
    type = ParsedFunction
    expression = '(M0 + M1*(x - T_0) + M2*(x - T_0)^2)'
    symbol_names ='M0 M1 M2 T_0'
    symbol_values ='${M0} ${M1} ${M2} ${T_0_in}'
  []
[]

# FLUID PROPS
[FluidProperties]
  [fp]
    type = TemperaturePressureFunctionFluidProperties
    rho = rho
    k = k
    mu = mu
    cp = cp
    T_ref = ${T_0_in}
    e_ref = ${e_0_in}
  []
[]

# Materials
[Materials]
  [Re_mat]
    type = ADReynoldsNumberMaterial
    Re = Re
    rho = rho
    vel = vel
    D_h = D_h
    mu = mu
    block = core_chan
  []
  [Pr_mat]
    type = ADPrandtlNumberMaterial
    Pr = Pr
    cp = cp
    mu = mu
    k = k
    block = core_chan
  []

[]

# CLOSURES
[Closures]
  [thm_closures]
    type = Closures1PhaseTHM
  []
[]

# COMPONENTS
# handled via Components blocks
# Components capture pipes, turbomachines, bcs, sources, etc
# They are lego pieces used to build model

# FlowChannel1Phase - a channel with a single phase
# position: location of channel start in 3D space
# orientation: direction vector of channel in 3D
# length: length of channel
# n_elems: number of elems to discretize
# A: CrossSec area
# D_h: hydraulic diameter
# Every flow channel defines a subdomain block using
# flow_channel_name:in and flow_channel_name:out
# NOTE: this does NOT set flow direction!!!

# BOUNDARY CONDITIONS
# each flow channel needs either a bc or a junction
# BCs are connected to flow channel via the input param

[Components]
  [inlet]
    type = InletMassFlowRateTemperature1Phase
    input = 'core_chan:in'
    m_dot = ${m_dot_in}
    T = ${T_in}
  []

  [core_chan]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '0 0 1'
    length = 1
    n_elems = 3
    A = ${A_in}
    D_h = ${D_h_in}
    f = 0.0
  []

  [outlet]
    type = Outlet1Phase
    input = 'core_chan:out'
    p = ${press}
  []
[]

# POSTPROCESSORS
# These compute single real values at locations
# here we want to watch pressure drop across channel
# SideAverageValue: compute on a side
[Postprocessors]
  [core_p_in]
    type = SideAverageValue
    boundary = core_chan:in
    variable = p
  []

  [core_p_out]
    type = SideAverageValue
    boundary = core_chan:out
    variable = p
  []

  # ParsedPostprocessor: define equation based on other pps
  [core_delta_p]
    type = ParsedPostprocessor
    pp_names = 'core_p_in core_p_out'
    expression = 'core_p_in - core_p_out'
  []

  [rho_out]
    type = SideAverageValue
    boundary = core_chan:out
    variable = rho
  []

  [Re_out]
    type = ADElementAverageMaterialProperty
    mat_prop = Re
    block = core_chan
  []

  [e_out]
    type = SideAverageValue
    boundary = core_chan:out
    variable = e 
  []
  [Pr_out]
    type = ADElementAverageMaterialProperty
    mat_prop = Pr 
    block = core_chan
  []

  [cv_out]
    type = ADElementAverageMaterialProperty
    mat_prop = cv 
    block = core_chan
  []
  [cp_out]
    type = ADElementAverageMaterialProperty
    mat_prop = cp 
    block = core_chan
  []
  [k_out]
    type = ADElementAverageMaterialProperty
    mat_prop = k 
    block = core_chan
  []
  [mu_out]
    type = ADElementAverageMaterialProperty
    mat_prop = mu 
    block = core_chan
  []

[]

# PRECONDITIONER
[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

# EXECUTIONER
# describes how to execute simulation
# type: Transient means over time steps
# NOTE: THM should alway use NEWTON with basic line search
# NOTE: also always reccoment given petsc options
[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = basic

  start_time = 0
  end_time = 1000
  dt = 500
  dtmin = 500

  petsc_options_iname = '-pc-type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6
  nl_max_its = 25
[]

# OUTPUTS
# Moose only prints to terminal by default
# need to setup output file if desired
[Outputs]
  exodus = true
  [console]
    type = Console
    max_rows = 1
    outlier_variable_norms = false
  []
  print_linear_residuals = false
  [csv]
    type = CSV
    file_base = FlowThroughChannel_MoltSalt_out 
    delimiter = ','
    show = 'core_delta_p'
  []
[]

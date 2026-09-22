from Run_ATES import run_DARTS
from Run_ATES import create_geomodel_confined_disc_opt
from Run_ATES import create_geomodel_uni_dxy
from plot_results import plot_timeseries, plot_cross_section

sim_name = "case_synthetic_no_obswell_gridext_v300_dens"

### Model dimensions
nx = 100 # correct would be 200
ny = 100 # correct would be 200
nz_aq = 20
nz_at_top = 2
nz_at_bot = 2

nx_large = 25 # correct would be 75
ny_large = 25 # correct would be 75

dx = 1
dy = 1
dz = 2
dz_at_top = 2
dz_at_bot = 2

dx_large = 4
dy_large = 4

start_z = 0

### Well locations (cell)
hwx = 50
hwy = 50

### Properties
permeability = 3.65 * 1e-11 # m2
permeability_mD = permeability * 1.01325 * 1e15 #1e12 m2 to D, 1e3
print("Horizontal permeability", permeability_mD, "mD")
anisotropy = 5
permeability_v = permeability_mD / anisotropy
print("Vertical permeability", permeability_v, "mD")

permeability_AT = 7.3 * 1e-14 # m2
permeability_mD_AT = permeability_AT * 1.01325 * 1e15
print("Horizontal permeability aquitard", permeability_mD_AT, "mD")
anisotropy_AT = 5
permeability_v_AT = permeability_mD_AT / anisotropy_AT
print("Vertical permeability aquitard", permeability_v_AT, "mD")

porosity = 0.3

solid_density = 2640 #kg/m3

solid_heat_cap = 710 #J/kg/K

solid_volumetric_hcap = (solid_density * solid_heat_cap)/ 1000 #kJ/m3
print("Volumetric heat capacity", solid_volumetric_hcap, "kJ/m3")

l_w = 0.58 #W/m/K
l_s = 2 #W/m/K
l_c = 1.7 #W/m/K

l_s_darts = l_s * ((24 * 3600)/1000) #kJ/K/day
l_c_darts = l_c * ((24 * 3600)/1000)
print("Thermal conductivity", l_s_darts, "kJ/K/day")
print("Thermal conductivity aquitard", l_c_darts, "kJ/K/day")

### Timestepping & Operational profile
# Block-function Operational profile
daysprofile = [90, 90, 90, 90]  # (days)
storage_periods = ['Charge', 'Rest', 'Discharge', 'Rest']

dt_max = 1 #day
dt_mult = 4

### Simulation settings
Tin = 20 # deg C
#Q_cell = 100000/(daysprofile[0]*nz_aq)  #100000/(daysprofile[0] * nz_aq) # m3/day #Note this is the flowrate per perforation
print("Injection time:", daysprofile[0])
Q_cell = 100000/(daysprofile[0])
print("Volumetric rate:", Q_cell)
well_diameter = 1 # m

#var_mod = create_geomodel_confined_disc_opt(nx, ny, nx_large, ny_large,
                                            #20, 11, 1,
                                            #dx, dy, dz, dx_large, dy_large,
                                            #0,
                                            #permeability_mD, permeability_v, porosity, volumetric_heat_cap, thermal_conductivity_DARTS,
                                            #permeability_mD_AT, permeability_v_AT, porosity, volumetric_heat_cap, thermal_conductivity_DARTS_AT,
                                            #int((nx/2) - 1 + nx_large), int((ny/2) - 1 + ny_large))

uni_mod = create_geomodel_uni_dxy(nx, ny, nz_at_top, nz_aq, nz_at_bot,
                                  dx, dy, dz, dz_at_top, dz_at_bot, start_z,
                                  permeability_mD, permeability_v, porosity, solid_volumetric_hcap, l_s_darts,
                                  permeability_mD_AT, permeability_v_AT, porosity, solid_volumetric_hcap, l_c_darts,
                                  hwx, hwy)


run_DARTS (sim_name,
           uni_mod,
           1,
           20, 5,
           Q_cell, daysprofile, storage_periods,
           set_transition_runtime = 1e-3,
           well_diameter = 1,
           n_points = 256) #ensure that n-points is correctly inferred in the other scripts

plot_timeseries(sim_name, "H_1")

time_step_list = [0, 1, 2, 3, 4]
plot_cross_section(sim_name,
                   time_step_list,
                   nx, ny, int(nz_at_top + nz_aq + nz_at_bot), hwx)

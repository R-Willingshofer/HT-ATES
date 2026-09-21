# Load model creation
from ATES_model import Model

import os
import xarray as xr
import pandas as pd
import numpy as np

def create_geomodel_uni_dxy(nx, ny,
                            nz_at_top, nz_aq, nz_at_bot,
                            dx, dy, dz, dz_at_top, dz_at_bot, start_z,
                            perm_h, perm_v, porosity, h_cap, t_cond,
                            perm_h_conf, perm_v_conf, porosity_conf, h_cap_conf, t_cond_conf,
                            hwx, hwy):

    "nz_conf includes the nz_actnum"

    #create dx arrays
    dx_array = np.tile(dx, nx)
    dy_array = np.tile(dy, ny)
    print(dx_array)

    x_centers = np.cumsum(dx_array) - 0.5 * dx
    y_centers = np.cumsum(dy_array) - 0.5 * dy

    dz_aq_array = np.tile(dz, nz_aq)
    dz_at_top_array = np.tile(dz_at_top, nz_at_top)
    dz_at_bot_array = np.tile(dz_at_bot, nz_at_bot)
    dz_array = np.concatenate([dz_at_top_array, dz_aq_array, dz_at_bot_array], axis = 0)
    print(dz_array)

    # Not entirely accurate
    z_centers = np.cumsum(dz_array) - 0.5 * dz

    # Property arrays aquifer
    permeability_xy_cells = np.full((nx, ny, nz_aq), perm_h)
    permeability_z_cells  = np.full((nx, ny, nz_aq), perm_v)
    porosity_cells        = np.full((nx, ny, nz_aq), porosity)
    heat_capacity_cells   = np.full((nx, ny, nz_aq), h_cap)
    thermal_conductivity_cells = np.full((nx, ny, nz_aq), t_cond)

    # Confining layers aquitard top
    at_top_permeability_xy_cells = np.full((nx, ny, nz_at_top), perm_h_conf)
    at_top_permeability_z_cells  = np.full((nx, ny, nz_at_top), perm_v_conf)
    at_top_porosity_cells        = np.full((nx, ny, nz_at_top), porosity_conf)
    at_top_heat_capacity_cells   = np.full((nx, ny, nz_at_top), h_cap_conf)
    at_top_thermal_conductivity_cells = np.full((nx, ny, nz_at_top), t_cond_conf)

    # Confining layers aquitard bottom
    at_bot_permeability_xy_cells = np.full((nx, ny, nz_at_bot), perm_h_conf)
    at_bot_permeability_z_cells  = np.full((nx, ny, nz_at_bot), perm_v_conf)
    at_bot_porosity_cells        = np.full((nx, ny, nz_at_bot), porosity_conf)
    at_bot_heat_capacity_cells   = np.full((nx, ny, nz_at_bot), h_cap_conf)
    at_bot_thermal_conductivity_cells = np.full((nx, ny, nz_at_bot), t_cond_conf)

    tot_permeability_xy = np.concatenate(
        [at_top_permeability_xy_cells, permeability_xy_cells, at_bot_permeability_xy_cells], axis = 2)
    tot_permeability_z = np.concatenate(
        [at_top_permeability_z_cells, permeability_z_cells, at_bot_permeability_z_cells], axis=2)
    tot_porosity = np.concatenate(
        [at_top_porosity_cells, porosity_cells, at_bot_porosity_cells], axis = 2)
    tot_heat_capacity = np.concatenate(
        [at_top_heat_capacity_cells, heat_capacity_cells, at_bot_heat_capacity_cells], axis = 2)
    tot_t_cond = np.concatenate(
        [at_top_thermal_conductivity_cells, thermal_conductivity_cells, at_bot_thermal_conductivity_cells], axis = 2)

    print("perm:", np.max(tot_permeability_xy), np.mean(tot_permeability_xy))
    print("vperm:", np.mean(tot_permeability_z), np.mean(tot_permeability_z))
    print("por:", np.mean(tot_porosity), np.mean(tot_porosity))
    print("heat:", np.mean(tot_heat_capacity), np.mean(tot_heat_capacity))
    print("tcond:", np.mean(tot_t_cond), np.mean(tot_t_cond))

    geomodel = xr.Dataset(
        data_vars={
            "permeability_xy": (("x", "y", "z"), tot_permeability_xy),
            "permeability_z":  (("x", "y", "z"), tot_permeability_z),
            "porosity":        (("x", "y", "z"), tot_porosity),
            "heat_capacity":   (("x", "y", "z"), tot_heat_capacity),
            "thermal_conductivity": (("x", "y", "z"), tot_t_cond)},

        coords={
            "x": x_centers,
            "y": y_centers,
            "z": z_centers},

        attrs={
            "dx" : dx_array,
            "dy" : dy_array,
            "dz" : dz_array,
            "hwx" : hwx,
            "hwy" : hwy,
            "nly_top" : nz_at_top,
            "nly_res" : nz_aq,
            "nly_bot" : nz_at_bot,
            "grid_type": "Uniform"})

    return geomodel

def run_DARTS (simulation_name,
               geomodel,
               nyears,
               Tin, TCutOff,
               volumetric_rate, operational_profile, storage_periods,
               set_transition_runtime = 1e-3,
               well_diameter = 1,
               n_points = 256):

    ### =========== Simulation settings ==============
    # ------------- Output directory ----------------
    output_directory = f"3D_{simulation_name}"  # Spatial output
    output_h5 = f"h5_{simulation_name}"  # DARTS output folder

    os.makedirs(output_h5, exist_ok=True)
    output_folder = output_h5

    output_well_data_excel = f"WD_{simulation_name}.xlsx"

    os.makedirs(output_directory, exist_ok=True)

    # ------------- Input geomodel -------------------
    perm_h = geomodel['permeability_xy'].values
    perm_v = geomodel['permeability_z'].values
    poro = geomodel['porosity'].values
    hcap = geomodel['heat_capacity'].values
    tcond = geomodel['thermal_conductivity'].values

    # Load XY-plane well indices
    hwx = geomodel.attrs['hwx']
    hwy = geomodel.attrs['hwy']

    dX_array = geomodel.attrs['dx']
    dY_array = geomodel.attrs['dy']
    dZ_array = geomodel.attrs['dz']

    nly_top = geomodel.attrs["nly_top"]
    nly_res = geomodel.attrs["nly_res"]
    nly_bot = geomodel.attrs["nly_bot"]

    print("nly:", nly_top, nly_res, nly_bot)

    # --------------- Initial conditions reservoir -------
    geothermal_grad = 0 #(K / km), the geothermal gradient for the initial condition of the reservoir

    # -------------- Operational parameters -----------------------
    # Well temperatures to Kelvin
    InjT = 273.15 + Tin # (K)
    #TCutOff = 273.15 + tcuo # (K), Serves as the injection temperature of the warm well

    depth_to_top = 0

    ### ============== Run Simulation =====================
    #Input model params here
    m = Model(dX_array, dY_array, dZ_array,
              nly_top, nly_res, nly_bot,
              perm_h, perm_v, poro, hcap, tcond,
              hwx, hwy, well_diameter,
              depth_to_top, geothermal_grad,
              n_points = n_points)
    m.init()
    m.set_output(output_folder = output_folder)

    iterr = 1
    for k in range(nyears):
        for i, runtime in enumerate(operational_profile):
            if storage_periods[i] == 'Charge':

                m.set_rate_hot(volumetric_rate, temp=InjT, func='inj')
                m.set_rate_cold(-1 * volumetric_rate, func='prod')
                m.set_rate_obs(0, func='prod')
                print('Operation: Charge')

            elif storage_periods[i] == 'Discharge':

                m.set_rate_hot(-1 * volumetric_rate, func='prod')
                m.set_rate_cold(volumetric_rate, temp=TCutOff, func='inj')
                m.set_rate_obs(0, func='prod')
                print('Operation: Discharge')

            elif storage_periods[i] == 'Rest':

                m.set_rate_hot(0, func='prod')
                m.set_rate_cold(0, func='prod')
                m.set_rate_obs(0, func='prod')
                print('Operation: Rest')

            m.run(runtime, restart_dt=set_transition_runtime)
            print("\nIterr :", iterr, "\tYear :", k, "\tRun Time :", runtime)
            print("\n")
            iterr += 1
    m.print_timers()
    m.print_stat()

    # Comment the following 2 lines to ensure no 3D output data is generated
    output_props = ['temperature', 'pressure']
    m.output.output_to_vtk(output_properties=output_props, output_directory = output_directory)


    # %%-----------------Write Results to Excel-----------------
    # output well information to Excel file
    td = pd.DataFrame.from_dict(m.physics.engine.time_data)
    #writer = pd.ExcelWriter(output_well_data_excel)
    #td.to_excel(writer, 'Sheet1')
    #writer.close()
    with pd.ExcelWriter(output_well_data_excel) as writer:
        td.to_excel(writer, sheet_name='Sheet1')

    # %%-----------------read H5-----------------
    #well_id, well_depth = write_well_perforation_id(m)
    #r = read_well_h5(well_block_id=well_id, well_block_depth=well_depth)
    #r.draw_combined_well_data()  # All wells on same subplots

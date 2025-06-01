from Basic_Ates_Model import Model
import pandas as pd
import xarray as xr
import os

#Code originally created by Taylan Akin
#Modified as part of the AES MSc thesis by Robin Willingshofer (2025)
#--------------Initialize arrays for batch simulations-------------------
geomod_names = ['PNA', 'RWK', 'BarSim']
geomod_ID = [1, 2, 3, 4]
Vin = [200000, 600000, 2000000]
Tinj = [90, 70]

#--------------Operational parameters--------------
Vin = Vin[0] #m3 per well per year
InjT = 273 + Tinj[0] # Injection temperature in Kelvin
TCutOff = 273 + 50

storage_periods=['Charge','Rest','Discharge','Rest']
daysprofile = [120, 60, 120, 65]
Q_hot = Vin/daysprofile[0]    #m3/day
Q_cold = Vin/daysprofile[0]    #m3/day
op_profile = [1, 0, 1, 0]

#-----------------Import geomodel-------------------
#geomod_name = f'{geomod_names}_{geomod_ID}_{Vin}_{InjT}'
#geomod_path = os.path.join("geomodels", 'PNA1_200k.nc')

# Load the dataset
geomodel = xr.load_dataset('PNA1_200k.nc')
XGR = geomodel['XGR'].values
YGR = geomodel['YGR'].values
ZGR = geomodel['ZGR'].values
porosity = geomodel['porosity'].values
horizontal_permeability = geomodel['permeability'].values
vertical_permeability = geomodel['kv2'].values
Cp = geomodel['Cv'].values
lam = geomodel['lam'].values

#----------------Output_dirs------------------------
output_directory = os.path.join(os.getcwd(), 'vtk_output_test')
os.makedirs(output_directory, exist_ok=True)
output_well_data_excel = os.path.join(output_directory, 'well_data_output_test.xlsx')

#-----------------Create DARTS ATES model--------------
m = Model(XGR, YGR, ZGR, porosity, horizontal_permeability, vertical_permeability, Cp, lam)
m.init()


#write initial conditions
m.output_to_vtk(ith_step=0, output_directory=output_directory)

set_transition_runtime = 1e-5
set_run_years = 3

flw_rates = [x * Q_hot for x in op_profile]
# coldrates = [x * Q_cold for x in op_profile]

iterr=1
for k in range(set_run_years):
    for i, runtime in enumerate(daysprofile):
        if storage_periods[i] == 'Charge':
            m.set_rate_hot(flw_rates[i], temp=InjT, func='inj')
            m.set_rate_cold(flw_rates[i], func='prod')
            print('Operation: Charge')

        elif storage_periods[i] == 'Discharge':
            m.set_rate_hot(flw_rates[i], func='prod')
            m.set_rate_cold(flw_rates[i], temp=TCutOff, func='inj')
            print('Operation: Discharge')

        elif storage_periods[i] == 'Rest':
            m.set_rate_hot(0, func='prod')
            m.set_rate_cold(0, func='prod')
            print('Operation: Rest')

        m.run(runtime, restart_dt=set_transition_runtime)
        m.output_to_vtk(ith_step=iterr, output_directory=output_directory)
        print("\nIterr :",iterr,"\tYear :",k, "\tRun Time :",runtime)
        print("\n")
        iterr+=1
m.print_stat()
#%%-----------------Write Results to Excel-----------------
# output well information to Excel file
td = pd.DataFrame.from_dict(m.physics.engine.time_data)
writer = pd.ExcelWriter(output_well_data_excel)
td.to_excel(writer, 'Sheet1')
writer.close()

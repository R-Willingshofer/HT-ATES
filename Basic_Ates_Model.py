from darts.physics.geothermal.physics import Geothermal
from darts.physics.geothermal.property_container import PropertyContainer
from darts.models.darts_model import DartsModel
from darts.engines import redirect_darts_output
from darts.reservoirs.cpg_reservoir import CPG_Reservoir
from darts.tools.gen_cpg_grid import gen_cpg_grid

redirect_darts_output('LogFile_Run_HT_ATES_DELFT_lam10.log')

from darts.engines import set_num_threads

set_num_threads(4)

import numpy as np

# %%
class Model(DartsModel):

    def __init__(self, XGR, YGR, ZGR, porosity, horizontal_permeability, vertical_permeability, Cp, lam, n_points=64):
        super().__init__()
        self.timer.node["initialization"].start()

        # === Grid dimensions ===
        self.nx = len(XGR) - 1
        self.ny = len(YGR) - 1
        self.nz = len(ZGR) - 1

        dx_list = np.diff(XGR)
        dy_list = np.diff(YGR)
        dz_list = np.diff(ZGR)

        self.n_ly_cap = 5  # layer number in cap rock
        self.n_ly_res = 120  # layer number in reservoir rock
        self.n_ly_bttm = 5  # layer number in bottom formation

        # === Flatten properties ===
        perm_x = horizontal_permeability.flatten()
        perm_y = horizontal_permeability.flatten()
        perm_z = vertical_permeability.flatten()
        poro = porosity.flatten()

        self.depth_to_top = 115  # [m]

        # === Discretize grid ===
        arrays = gen_cpg_grid(
            nx=self.nx, ny=self.ny, nz=self.nz,
            dx=dx_list, dy=dy_list, dz=dz_list,
            start_z=self.depth_to_top,
            permx=perm_x, permy=perm_y, permz=perm_z, poro=poro
        )

        self.reservoir = CPG_Reservoir(self.timer, arrays)
        self.reservoir.discretize()

        # === Thermal properties ===
        self.reservoir.hcap[:] = Cp.flatten()
        self.reservoir.conduction[:] = lam.flatten()

        #self.reservoir.set_boundary_volume(yz_minus=1e15, yz_plus=1e15, xz_minus=1e15, xz_plus=1e15)
        self.reservoir.set_boundary_volume(
            xy_minus=1e15,
            xy_plus=1e15,
            yz_minus=1e15,
            yz_plus=1e15,
            xz_minus=1e15,
            xz_plus=1e15
        )

        # === Physics ===
        property_container =PropertyContainer()
        self.physics = Geothermal(self.timer, n_points, 0.1, 150, 500, 15000, cache=False)
        self.physics.add_property_region(property_container)
        self.physics.init_physics()

        # === Time-stepping ===
        self.params.first_ts = 1e-5
        self.params.mult_ts = 2
        self.params.max_ts = 30

        # === Solver tolerances ===
        self.params.tolerance_newton = 1e-4

        self.timer.node["initialization"].stop()


    def set_wells(self):
        top_ind = self.n_ly_cap + 1
        btm_ind = self.n_ly_cap + self.n_ly_res

        self.reservoir.add_well("H1")
        self.reservoir.add_well("L1")

        for i in range(top_ind, btm_ind + 1):
            print("ind: ", i)
            # Change the index of the well location to the actual thing, also find a more elegant way of providing this (in the gridding function ideally)
            self.reservoir.add_perforation("H1", cell_index=(int(19), int(self.ny / 2), i), verbose=True)
            self.reservoir.add_perforation("L1", cell_index=(int(51), int(self.ny / 2), i), verbose=True)


    def set_initial_conditions(self):
        self.physics.set_nonuniform_initial_conditions(self.reservoir.mesh,
                                                   pressure_grad=100,  # bar/km
                                                   temperature_grad=18,  # C/km
                                                   T_at_ref_depth=(273.15 + 10))  # Temp at Surface 10 Celcius


    def set_boundary_conditions(self):
        # activate wells with rate control for inejctor and producer
        for i, w in enumerate(self.reservoir.wells):
        # if 'INJ' in w.name:
            if 'H' in w.name:
                w.control = self.physics.new_rate_water_prod(0)
            else:
                w.control = self.physics.new_rate_water_prod(0)

    # def set_rate_hot(self, rate, welln=0, temp=300, func='inj'):


    def set_rate_hot(self, rate, temp=300, func='inj'):
    # w = self.reservoir.wells[welln]
        for w in self.reservoir.wells:
            if 'H' in w.name:
                if func == 'inj':
                    w.control = self.physics.new_rate_water_inj(rate, temp)
                # w.constraint = self.physics.new_bhp_water_inj(self.midrespress + self.bhp_limit, temp)
                if func == 'prod':
                    w.control = self.physics.new_rate_water_prod(rate)
                # w.constraint = self.physics.new_bhp_prod(self.midrespress - self.bhp_limit)

    # def set_rate_cold(self, rate, welln=1, temp=300, func='inj'):


    def set_rate_cold(self, rate, temp=300, func='inj'):
    # w = self.reservoir.wells[welln]
        for w in self.reservoir.wells:
            if 'L' in w.name:
                if func == 'inj':
                    w.control = self.physics.new_rate_water_inj(rate, temp)
                # w.constraint = self.physics.new_bhp_water_inj(self.midrespress + self.bhp_limit, temp)
                if func == 'prod':
                    w.control = self.physics.new_rate_water_prod(rate)
                # w.constraint = self.physics.new_bhp_prod(self.midrespress - self.bhp_limit)
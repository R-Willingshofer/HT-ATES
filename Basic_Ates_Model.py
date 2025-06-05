from darts.discretizer import value_vector
from darts.physics.geothermal.physics import Geothermal
from darts.physics.geothermal.property_container import PropertyContainer
from darts.models.darts_model import DartsModel
from darts.engines import redirect_darts_output, well_control_iface, sim_params
from darts.reservoirs.struct_reservoir import StructReservoir

redirect_darts_output('LogFile_Run_HT_ATES_DELFT.log')

from darts.engines import set_num_threads

set_num_threads(4)

import numpy as np


# %%
class Model(DartsModel):

    def __init__(self, XGR, YGR, ZGR, porosity, horizontal_permeability, vertical_permeability, Cp, lam, hwx, hwy, cwx, cwy, n_points=64):
        super().__init__()
        self.timer.node["initialization"].start()

        #Set Well indices (test)
        self.hwx = hwx
        self.hwy = hwy
        self.cwx = cwx
        self.cwy = cwy

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

        self.reservoir = StructReservoir(self.timer, nx=self.nx, ny=self.ny, nz=self.nz,
                                         dx=dx_list, dy=dy_list, dz=dz_list,
                                         permx=perm_x, permy=perm_y, permz=perm_z,
                                         poro=poro, start_z=self.depth_to_top)
        # === Thermal properties ===
        self.reservoir.hcap[:] = Cp.flatten()
        self.reservoir.rcond[:] = lam.flatten()

        #
        self.reservoir.boundary_volumes['yz_minus'] = 1e20
        self.reservoir.boundary_volumes['yz_plus'] = 1e20
        self.reservoir.boundary_volumes['xz_minus'] = 1e20
        self.reservoir.boundary_volumes['xz_plus'] = 1e20

        # ------------create pre-defined physics for geothermal------------
        property_container = PropertyContainer()
        self.physics = Geothermal(self.timer, n_points, 0.1, 150, 500, 7500, cache=False)
        self.physics.add_property_region(property_container)
        self.physics.init_physics()

        self.set_sim_params(first_ts=1e-3, mult_ts=8, max_ts=30, runtime=3650, tol_newton=1e-4, tol_linear=1e-8,
                            it_newton=20, it_linear=40, newton_type=sim_params.newton_global_chop,
                            newton_params=value_vector([1]))
        self.timer.node["initialization"].stop()

    def set_wells(self):
        top_ind = self.n_ly_cap + 1
        btm_ind = self.n_ly_cap + self.n_ly_res

        self.reservoir.add_well("H1")
        self.reservoir.add_well("L1")

        for i in range(top_ind, btm_ind + 1):
            self.reservoir.add_perforation("H1",
                                           cell_index=(self.hwx, self.hwy, i),
                                           verbose=True,
                                           multi_segment=False,
                                           well_indexD=0)
                                           # well_index=50,
                                           # well_indexD=0)
            self.reservoir.add_perforation("L1",
                                           cell_index=(self.cwx, self.cwy, i),
                                           verbose=True,
                                           multi_segment=False,
                                           well_indexD=0)
                                           #  well_index=50,
                                           # well_indexD=0)


    def set_initial_conditions(self):
        input_depth = [0., np.amax(self.reservoir.mesh.depth)]
        input_distribution = {'pressure': [1., 1. + input_depth[1] * 100. / 1000],
                              'temperature': [283.15, 283.15 + input_depth[1] * 18. / 1000]}
        return self.physics.set_initial_conditions_from_depth_table(self.reservoir.mesh,
                                                                    input_distribution=input_distribution,
                                                                    input_depth=input_depth)
    def set_well_controls(self, h_func=None, l_func=None):
        for i, w in enumerate(self.reservoir.wells):
            if 'H' in w.name:
                self.physics.set_well_controls(wctrl=w.control, control_type=well_control_iface.VOLUMETRIC_RATE,
                                               is_inj=True, target=0., phase_name='water', inj_composition=[],
                                               inj_temp=273.15 + 90)
            else:
                self.physics.set_well_controls(wctrl=w.control, control_type=well_control_iface.VOLUMETRIC_RATE,
                                               is_inj=False, target=0., phase_name='water')

    def set_rate_hot(self, rate, temp=300, func='inj'):
        for w in self.reservoir.wells:
            if 'H' in w.name:
                if func == 'inj':
                    self.physics.set_well_controls(wctrl=w.control,
                                                   control_type=well_control_iface.VOLUMETRIC_RATE,
                                                   is_inj=True, target=rate, phase_name='water', inj_composition=[],
                                                   inj_temp=temp)
                if func == 'prod':
                    self.physics.set_well_controls(wctrl=w.control,
                                                   control_type=well_control_iface.VOLUMETRIC_RATE,
                                                   is_inj=False, target=rate, phase_name='water')

    def set_rate_cold(self, rate, temp=300, func='inj'):
        # w = self.reservoir.wells[welln]
        for w in self.reservoir.wells:
            if 'L' in w.name:
                if func == 'inj':
                    self.physics.set_well_controls(wctrl=w.control,
                                                   control_type=well_control_iface.VOLUMETRIC_RATE,
                                                   is_inj=True, target=rate, phase_name='water', inj_composition=[],
                                                   inj_temp=temp)

                if func == 'prod':
                    self.physics.set_well_controls(wctrl=w.control,
                                                   control_type=well_control_iface.VOLUMETRIC_RATE,
                                                   is_inj=False, target=rate, phase_name='water')

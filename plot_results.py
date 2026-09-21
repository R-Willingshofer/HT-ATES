import os
import matplotlib.pyplot as plt
import pandas as pd
import pyvista as pv
import numpy as np

def plot_timeseries (sim_name, well_name):
    df = pd.read_excel(f'WD_{sim_name}.xlsx')
    time = df['time']
    rate = df[f'{well_name} : L rate (m3/day)']
    temp = df[f'{well_name} : temperature (K)']
    bhpr = df[f'{well_name} : BHP (bar)']

    fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, figsize = (6,8))
    ax[0].plot(time, rate)
    ax[0].set_ylabel('(m3/day)')
    ax[0].set_title('Flow rate')
    ax[0].grid()

    ax[1].plot(time, temp - 273.15)
    ax[1].set_ylabel('(deg C)')
    ax[1].set_title('Temperature')
    ax[1].grid()

    ax[2].plot(time, bhpr)
    ax[2].set_xlabel('Time (days)')
    ax[2].set_title('BHP')
    ax[2].set_ylabel('(bar)')
    ax[2].grid()

    plt.show()

def plot_cross_section (sim_name, time_step_list,
                        nx, ny, nz, hwx):

    folder = f'3D_{sim_name}'

    # Read all sections first
    sections = []

    for time_step in time_step_list:

        file_name = f'solution_ts{time_step}.vts'
        path = os.path.join(folder, file_name)
        grid = pv.read(path)
        temp = (np.asarray(grid.cell_data["temperature[K]"]).reshape((nz, ny, nx)).transpose(2, 1, 0))

        # K -> °C
        temp -= 273.15

        # YZ section through well
        sections.append(temp[hwx, :, :].T)

    # Common temperature limits
    vmin = min(section.min() for section in sections)
    vmax = max(section.max() for section in sections)

    fig, ax = plt.subplots(
        nrows=len(time_step_list),
        ncols=1,
        sharex=True,
        figsize=(6, 2.5 * len(time_step_list)),
        squeeze=False)

    ax = ax[:, 0]

    for i, (time_step, section_array) in enumerate(zip(time_step_list, sections)):

        im = ax[i].imshow(section_array, aspect='auto', origin='upper',vmin=vmin,  vmax=vmax, cmap = 'coolwarm')

        ax[i].set_ylabel('Z (cell index)')
        ax[i].set_title(f'Time step {time_step}')

    ax[-1].set_xlabel('Y (cell index)')

    fig.colorbar(im, ax=ax, label='Temperature (°C)', pad=0.02)

    plt.tight_layout()
    plt.show()

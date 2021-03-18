
import numpy as np
from matplotlib import pyplot as plt

def plot_e(sim_dir, simulation_dictionary, species, step_ii,flim, simulation_has_run = True, do_save = False):
    sd = simulation_dictionary
    output_dir = 'simulation_output/' + species + '/'
    sim_dir_str = ''
    if sim_dir is not None:
        sim_dir_str = sim_dir
        if sim_dir[-1] != '/':
            sim_dir_str += '/'
        output_dir = sim_dir_str + output_dir

    if not simulation_has_run:
        print('unable to plot; simulation has not run or had errors')
        return
    
    simtime = step_ii * sd["dt"]
    
    # output_dir = sim_dir + 'simulation_output/'
    xs = np.fromfile(output_dir + f'xs/xs_{step_ii}')
    es = np.fromfile(output_dir + f'es/es_{step_ii}')
    
    fig, ax1 = plt.subplots()
    ax1.set_title(f't={simtime:.03f}')
    ax1.plot(xs,es,'.')
    ax1.set_xlabel('x')

    ax1.set_xlim(sd['xmin'], sd['xmax'])
    # ax1.set_ylim(sd['pmin'],sd['pmax'])

    ax1.grid()
    plt.tight_layout()
    if do_save:
        plt.savefig(sim_dir_str + f'e_field_time_{step_ii}.png')
    plt.close()
# end plot_e

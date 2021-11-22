"""
Collection of functions for creating and running FARRSIGHT simulations

Contains
---
SimType : Enumerated type
make_dirs [DEPRECATED]:
generate_standard_names_dirs :
dict_to_deck :
deck_to_dict :
update_dictionary:
run_sim :
plot_phase_space :
phase_movie :
logf_movie :
panel_height_movie : 
sim_diagnostics_sample :

Dependencies
---
standard Python distribution (i.e. through Anaconda)

ffmpeg for movie writing

"""

import json # dump, load
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import os # path.exists, makedirs, getcwd
import shutil # rmtree, copy2
import subprocess # run, PIPE
import sys # path.append
import time # time

can_do_movie = True
try: 
    from matplotlib.animation import FFMpegWriter
    import matplotlib.animation as manimation
except:
    print('Unable to load ffmpeg.  Movie writer not accessible')
    can_do_movie = False


import matplotlib
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Polygon
plt.rcParams.update({'font.size': 18})

FIG_DPI = 500
MOV_DPI = 300
LOW_DPI = 200
panels_cmap = 'Greys'

from enum import IntEnum
class BoundaryConditions(IntEnum):
    PERIODIC = 0
    OPEN = 1
class SimType(IntEnum):
    WEAK_LD = 1
    STRONG_LD = 2
    STRONG_TWO_STREAM = 3
    COLDER_TWO_STREAM = 4
    FRIEDMAN_BEAM = 5
class Quadrature(IntEnum):
    TRAP = 0
    SIMPSONS = 1

sim_type_to_flim = {SimType.WEAK_LD : (0, 0.44),
                SimType.STRONG_LD : (0, 0.47),
                SimType.STRONG_TWO_STREAM : (0, 0.45),
                SimType.COLDER_TWO_STREAM : (0, 0.8)}

def make_dirs(project_name, sim_group, sim_name, 
                          root_dir = None):
    """
    [DEPRECATED] 
    Make directories to hold simulation data and return simulation directory
    
    If any of the following directories don't exist, make them:
    1. root_dir
    2. root_dir/project_name
    3. root_dir/project_name/sim_group
    
    If sim_dir = root_dir/project_name/sim_group/sim_name exists,
    then it possibly contains old data.  To eliminate confusion,
    this is erased and recreated.
    
    Inside sim_dir, create
    1. es
    2. fs
    3. xs
    4. vs
    5. q_ws
    6. panels
    directories
    
    Parameters
    ---
    project_name
    sim_group
    sim_name
    root_dir : optional, where all this is stored.  If not provided, data is created in working directory
    
    Returns
    ---
    sim_dir : string, <root_dir>/<project_name>/<sim_group>/<sim_name>/
    directories_found : boolean
    """
    directories_found = True
    if root_dir is not None:
        if root_dir[-1] != '/':
            root_dir += '/'
        simulations_dir = root_dir + 'simulations/' 
    else:
        simulations_dir = 'simulations/'
    proj_dir = simulations_dir + project_name + '/'
    sim_group_dir = proj_dir + sim_group + '/'
    sim_dir = sim_group_dir + sim_name + '/'

    if root_dir is not None:
        if not os.path.exists(root_dir):
            directories_found = False
            os.makedirs(root_dir)
    if not os.path.exists(simulations_dir):
        directories_found = False
        os.makedirs(simulations_dir)
    if not os.path.exists(proj_dir):
        directories_found = False
        os.makedirs(proj_dir)
    if not os.path.exists(sim_group_dir):
        directories_found = False
        os.makedirs(sim_group_dir)
    if not os.path.exists(sim_dir):
        directories_found = False
        os.makedirs(sim_dir)

    return sim_dir, directories_found

def generate_standard_names_dirs(simulation_dictionary, root_dir=None):

    directories_found = True
    if root_dir is not None:
        if root_dir[-1] != '/':
            root_dir += '/'
        simulations_dir = root_dir + 'simulations/' 
    else:
        simulations_dir = 'simulations/'

    sd = simulation_dictionary
    bcs_string = BoundaryConditions(0).name
    if 'bcs' in sd:
        bcs_string = BoundaryConditions(sd['bcs']).name
    bcs_string += '_bcs'
        
    tc_string = ''
    # if 'use_treecode' not in sd
    if 'use_treecode' in sd and sd['use_treecode']:
        if 0 <= sd['beta'] <= 1:
            tc_string = '_beta_%.2f'%sd['beta']# f'_beta_{sim_dict['beta']:.2f}'
        else:
            tc_string = '_mac_%.2f_degree_%d_msource_%d_maxtarget_%d'%(sd['mac'],sd['degree'], sd['max_source'], sd['max_target'])
    else:
        tc_string = '_dsum'

    if 'adaptively_refine' in sd and sd['adaptively_refine']:
        amr_string = 'amr_max_height_%i'%sd['max_height']
        amr_string += '_epsilons'
        for eps in sd['amr_epsilons']:
            amr_string += '_%.04f'%(eps)
    else:
        amr_string = 'no_amr'
#         sim_group = ''
#         sim_name = ''
    tf = sd['dt'] * sd['num_steps']
    physical_parameters = SimType(sd['sim_type']).name + '_' + bcs_string + '_vth_%.3f_vstr_%.3f_amp_%.3f_normal_k_%.3f_tf_%.1f'%(sd['vth'], sd['vstr'], sd['amp'], sd['normalized_wavenumber'], tf)
    if 'remesh_period' in sd:
        remesh_period = sd['remesh_period']
    else:
        remesh_period = 1
    if 'quadrature' in sd:
        quad_str = Quadrature(sd['quadrature']).name
    else:
        quad_str = 'TRAP'
    if 'v_height' in sd:
        v_height = sd['v_height']
    else:
        v_height = 0
    numerical_parameters = '%s_quadrature_height0_%i_v_height_%i_vm_%.1f_g_eps_%.5f_dt_%.4f_remesh_period_%i_diag_freq_%i'%(quad_str,sd['initial_height'], v_height, sd['vmax'], sd['greens_epsilon'], sd['dt'], remesh_period,sd['diag_period'])
    amr_treecode_paramters = amr_string + tc_string
#         sim_name = f'height0_{self.initial_height}_vm_{self.vmax:.1f}_g_eps_{self.greens_epsilon:.3f}_dt_{self.dt:.3f}_tf_{self.tf:.1f}_diag_freq_{self.diag_freq}' + tc_string
    sim_dir = simulations_dir + sd['project_name'] + '/' + physical_parameters + '/'
    sim_dir += numerical_parameters + '/' + amr_treecode_paramters + '/'
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)
    # sim_dir, directories_found = make_dirs(sd['project_name'], sim_group, sim_name,root_dir=root_dir)
    return sim_dir, directories_found
# end generate_standard_names_dirs

def dict_to_deck(sim_dictionary,deck_dir=None, deck_name='deck.json'):
    """writes contents of a simulation dictionary to an input deck file

    Parameters
    ----------
    sim_dictionary : dictionary
        contains simulation data
    deck_dir : str, optional
        location of input deck, by default None
    deck_name : str, optional
        name of input deck, by default 'deck.json'
    """
    if deck_dir is None:
        input_deck = deck_name
    else:
        if deck_dir[-1] != '/':
            deck_dir += '/'
        input_deck = deck_dir + deck_name
    with open(input_deck, 'w') as outfile:
        json.dump(sim_dictionary, outfile)

def deck_to_dict(deck_dir = None, deck_name = None):
    """
    open input deck .json file and load data into a dictionary
    
    Returns
    ---
    simulation_dictionary : dictionary, contains data from input deck
    """
    if deck_dir is not None:
        if deck_dir[-1] != '/':
            deck_dir += '/'
    else:
        deck_dir = ''
    if deck_name is not None:
        input_deck = deck_dir + deck_name
    else:
        input_deck = deck_dir + 'deck.json'
    with open(input_deck, 'r') as infile:
        simulation_dictionary = json.load(infile)
    return simulation_dictionary
# end deck_to_dict

def update_dictionary(deck_dir = None, deck_name = None, 
                from_existing = True, **sim_vars):
    simulation_dictionary = {}
    if from_existing:
        try:
            simulation_dictionary = deck_to_dict(deck_dir, deck_name)
        except FileNotFoundError:
            print('unable to find input deck ', deck_name, ', looking in ', deck_dir)
            print('Creating default dictionary')
            from_existing = False
    if not from_existing:
        simulation_dictionary = \
        {'project_name':'template',\
        'xmin' : 0.0, 'xmax' : 4*np.pi,\
        'vmin' : -6.0, 'vmax' : 6.0,\
        'bcs' : BoundaryConditions.PERIODIC.value,\
        'sim_type' : 2,\
        'normalized_wavenumber':1.0,\
        'amp':0.5, 'vth' : 1.0, 'vstr':0.0,\
        'initial_height' : 5,\
        'v_height' : 0,\
        'max_height' : 5,\
        'greens_epsilon' : 0.2,\
        'quadrature' : 0,\
        'use_treecode' : 0,\
        'beta' : -1.0,\
        'mac' : 0.8,\
        'degree' : 4,\
        'max_source' : 200,\
        'max_target' : 200,\
        'num_steps' : 60,\
        'remesh_period' : 1,\
        'diag_period' :1,\
        'dt' : 0.5,\
        'adaptively_refine' : 0,\
        'amr_epsilons':[]}
    simulation_dictionary.update(sim_vars)
    dict_to_deck(simulation_dictionary, deck_dir=deck_dir, deck_name=deck_name)
    return simulation_dictionary
# end of update_dictionary

def run_sim(sim_dir=None, deck_dir = None, deck_name = None, force_deck_cwd = False, use_gpu = False):
    sim_dir_str = ''
    if sim_dir is not None:
        sim_dir_str = sim_dir
        if sim_dir[-1] != '/':
            sim_dir_str += '/'
    # if deck_dir = None:
    #     deck_dir = sim_dir
    # elif deck_dir[-1] != '/':
    #     deck_dir += '/'
    # if deck_name is None:
    #     input_deck = deck_dir + 'deck.json'
    # else:
    #     input_deck = deck_dir + deck_name

    output_dir = sim_dir_str + 'simulation_output/'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    if not os.path.exists(output_dir + 'panels/'):
        os.makedirs(output_dir + 'panels/')
    if not os.path.exists(output_dir + 'xs/'):
        os.makedirs(output_dir + 'xs/')
    if not os.path.exists(output_dir + 'vs/'):
        os.makedirs(output_dir + 'vs/')
    if not os.path.exists(output_dir + 'es/'):
        os.makedirs(output_dir + 'es/')
    if not os.path.exists(output_dir + 'fs/'):
        os.makedirs(output_dir + 'fs/')
    if not os.path.exists(output_dir + 'qws/'):
        os.makedirs(output_dir + 'qws/')

    if use_gpu:
        farrsight_exe = 'farrsight_gpu'
    else:
        farrsight_exe = 'farrsight_cpu'
    farrsight_args = [farrsight_exe]


    deck_dir_str = ''
    input_deck = ''
    if deck_dir is not None:
        deck_dir_str = deck_dir
        if deck_dir[-1] != '/':
            deck_dir_str += '/'
    elif not force_deck_cwd:
        deck_dir_str = sim_dir_str

    if deck_name is not None:
        input_deck = deck_dir_str + deck_name
    else:
        input_deck = deck_dir_str + 'deck.json'
    # print('input deck', input_deck)
    # print('sim_dir', sim_dir)
    if sim_dir_str != deck_dir_str:
        shutil.copy2(input_deck, sim_dir_str)

    if sim_dir is not None:
        farrsight_args.append(sim_dir_str)
        if  deck_dir is not None and not force_deck_cwd:
            farrsight_args.append(deck_dir_str)
        if deck_name is not None:
            farrsight_args.append(deck_name)
    farrsight_args = [str(arg) for arg in farrsight_args]

    print('farrsight_args', farrsight_args)

    print('running sim...')
    t1 = time.time()
    with open(sim_dir_str + 'sim_out.txt','w') as log:
        with open(sim_dir_str + 'sim_err.txt','w') as err_log:
            proc = subprocess.run(farrsight_args,stdout=log,stderr=err_log)
    t2 = time.time()

    if proc.returncode == 0:
        print('done!')
        simulation_has_run = True
    else:
        simulation_has_run = False
        print('simulation stopped with error code %i'%proc.returncode)

    with open(sim_dir_str + 'sim_out.txt','a') as log:
        log.write(f'compute time {t2-t1:.03f}s')
    print(f'simulation compute time {t2-t1:.03f}s')
    return simulation_has_run
# end run_sim

def run_from_dict_standard_tree(simulation_dictionary, root_dir=None, use_gpu=False):

    sim_dir, directories_found = generate_standard_names_dirs(simulation_dictionary,root_dir)
    dict_to_deck(simulation_dictionary, deck_dir=sim_dir)
    simulation_has_run = run_sim(sim_dir, use_gpu=use_gpu)
    return sim_dir, directories_found, simulation_has_run


def plot_phase_space(sim_dir, simulation_dictionary, step_ii,flim, simulation_has_run = True, do_save = False):
    sd = simulation_dictionary
    output_dir = 'simulation_output/'
    sim_dir_str = ''
    if sim_dir is not None:
        sim_dir_str = sim_dir
        if sim_dir[-1] != '/':
            sim_dir_str += '/'
        output_dir = sim_dir_str + 'simulation_output/'

    if not simulation_has_run:
        print('unable to plot; simulation has not run or had errors')
        return
    
    simtime = step_ii * sd["dt"]
    
    # output_dir = sim_dir + 'simulation_output/'
    xs = np.fromfile(output_dir + f'xs/xs_{step_ii}')
    vs = np.fromfile(output_dir + f'vs/vs_{step_ii}')
    fs = np.fromfile(output_dir + f'fs/fs_{step_ii}')
    es = np.fromfile(output_dir + f'es/es_{step_ii}')
    panels = np.fromfile(output_dir + f'panels/leaf_point_inds_{step_ii}',dtype='int32')
    num_panels = int(panels.size/9)
    panels = np.reshape(panels, (num_panels,9))
    panels_fs = np.zeros(4*num_panels)
    
    
    fig,ax = plt.subplots(figsize=(8,6))
    ax0 = ax

    patches = []
    for ii, panel in enumerate(panels):
        panel_xs = xs[panel]
        panel_vs = vs[panel]
        panel_fs = fs[panel]
#             weights = np.array([1,4,1,4,16,4,1,4,1])
#             panels_fs[ii] = 1./36. * np.dot(panel_fs,weights)

        p0 = [0,1,4,3]
        panels_fs[4*ii] = .25*sum(panel_fs[p0])
        rect_pts = np.vstack([panel_xs[p0],panel_vs[p0]]).T
        patches.append(Polygon(rect_pts))

        p1 = [1,2,5,4]
        panels_fs[4*ii+1] = .25*sum(panel_fs[p1])
        rect_pts = np.vstack([panel_xs[p1],panel_vs[p1]]).T
        patches.append(Polygon(rect_pts))

        p2 = [3,4,7,6]
        panels_fs[4*ii+2] = .25*sum(panel_fs[p2])
        rect_pts = np.vstack([panel_xs[p2],panel_vs[p2]]).T
        patches.append(Polygon(rect_pts))

        p3 = [4,5,8,7]
        panels_fs[4*ii+3] = .25*sum(panel_fs[p3])
        rect_pts = np.vstack([panel_xs[p3],panel_vs[p3]]).T
        patches.append(Polygon(rect_pts))

    p = PatchCollection(patches, cmap=matplotlib.cm.jet)
#     if do_show_panels:
#         p.set_ec('black')
    p.set_array(panels_fs)
    p.set_clim(flim)
    ax0.add_collection(p)
    cb = fig.colorbar(p, ax=ax0)
    # cb.set_label('f')

    ax0.set_xlim(sd['xmin'], sd['xmax'])
    ax0.set_ylim(sd['vmin'],sd['vmax'])
    ax0.set_xlabel('x')
    ax0.set_ylabel('v')
    # ax0.set_title(f't={simtime:.03f}')
    plt.tight_layout()
    
    plt.savefig(sim_dir_str + f'phase_space_{simtime:.2f}.pdf')
    plt.savefig(sim_dir_str + f'phase_space_{simtime:.2f}.png',dpi=FIG_DPI)
    plt.close()
# end plot_phase_space


def plot_logf(sim_dir, simulation_dictionary, step_ii,flim=(-20,0), simulation_has_run = True, do_save = False):
    sd = simulation_dictionary
    output_dir = 'simulation_output/'
    sim_dir_str = ''
    if sim_dir is not None:
        sim_dir_str = sim_dir
        if sim_dir[-1] != '/':
            sim_dir_str += '/'
        output_dir = sim_dir_str + 'simulation_output/'

    if not simulation_has_run:
        print('unable to plot; simulation has not run or had errors')
        return
    
    simtime = step_ii * sd["dt"]
    
    # output_dir = sim_dir + 'simulation_output/'
    xs = np.fromfile(output_dir + f'xs/xs_{step_ii}')
    vs = np.fromfile(output_dir + f'vs/vs_{step_ii}')
    fs = np.fromfile(output_dir + f'fs/fs_{step_ii}')
    es = np.fromfile(output_dir + f'es/es_{step_ii}')
    panels = np.fromfile(output_dir + f'panels/leaf_point_inds_{step_ii}',dtype='int32')
    num_panels = int(panels.size/9)
    panels = np.reshape(panels, (num_panels,9))
    panels_fs = np.zeros(4*num_panels)
    
    
    fig,ax = plt.subplots(figsize=(8,6))
    ax0 = ax

    patches = []
    for ii, panel in enumerate(panels):
        panel_xs = xs[panel]
        panel_vs = vs[panel]
        panel_fs = fs[panel]
#             weights = np.array([1,4,1,4,16,4,1,4,1])
#             panels_fs[ii] = 1./36. * np.dot(panel_fs,weights)

        p0 = [0,1,4,3]
        panels_fs[4*ii] = .25*sum(panel_fs[p0])
        rect_pts = np.vstack([panel_xs[p0],panel_vs[p0]]).T
        patches.append(Polygon(rect_pts))

        p1 = [1,2,5,4]
        panels_fs[4*ii+1] = .25*sum(panel_fs[p1])
        rect_pts = np.vstack([panel_xs[p1],panel_vs[p1]]).T
        patches.append(Polygon(rect_pts))

        p2 = [3,4,7,6]
        panels_fs[4*ii+2] = .25*sum(panel_fs[p2])
        rect_pts = np.vstack([panel_xs[p2],panel_vs[p2]]).T
        patches.append(Polygon(rect_pts))

        p3 = [4,5,8,7]
        panels_fs[4*ii+3] = .25*sum(panel_fs[p3])
        rect_pts = np.vstack([panel_xs[p3],panel_vs[p3]]).T
        patches.append(Polygon(rect_pts))

    p = PatchCollection(patches, cmap=matplotlib.cm.jet)
    p.set_array(np.log10(panels_fs))
    p.set_clim(flim)
    ax0.add_collection(p)
    cb = fig.colorbar(p, ax=ax0)
    # cb.set_label('f')

    ax0.set_xlim(sd['xmin'], sd['xmax'])
    ax0.set_ylim(sd['vmin'],sd['vmax'])
    ax0.set_xlabel('x')
    ax0.set_ylabel('v')
    # ax0.set_title(f't={simtime:.03f}')
    plt.tight_layout()
    
    plt.savefig(sim_dir_str + f'logf_{simtime:.2f}.pdf')
    plt.savefig(sim_dir_str + f'logf_{simtime:.2f}.png',dpi=FIG_DPI)
    plt.close()
# end plot_logf

def plot_height(sim_dir, simulation_dictionary, iter_num, height_range = [7,11],simulation_has_run=True):

    sd = simulation_dictionary
    Lx = sd['xmax'] - sd['xmin']
    output_dir = 'simulation_output/'
    sim_dir_str = ''
    if sim_dir is not None:
        sim_dir_str = sim_dir
        if sim_dir[-1] != '/':
            sim_dir_str += '/'
        output_dir = sim_dir_str + 'simulation_output/'
    print('starting panel height plot')
    t1 =time.time()
    if not simulation_has_run:
        print('unable to plot; simulation has not run or had errors')
        return

    fig, ax = plt.subplots(figsize=(8,6))

            
    ncolors = height_range[1] - height_range[0] + 1
    mymap = cm.get_cmap(panels_cmap,ncolors)
    
    ax.set_xlim([sd['xmin'], sd['xmax']])
    ax.set_ylim([sd['vmin'],sd['vmax']])
    ax.set_xlabel('x')
    ax.set_ylabel('v')


    panels = np.fromfile(output_dir + f'panels/leaf_point_inds_{iter_num}',dtype='int32')
    num_panels = int(panels.size/9)
    panels = np.reshape(panels, (num_panels,9))
    panels_fs = np.zeros(num_panels)
    xs = np.fromfile(output_dir + f'xs/xs_{iter_num}')
    vs = np.fromfile(output_dir + f'vs/vs_{iter_num}')

    patches = []
    pvert = [0,2,8,6]
    for ii, panel in enumerate(panels):
        panel_xs = xs[panel]
        panel_vs = vs[panel]
        dx = panel_xs[8] - panel_xs[0]
        
        panels_fs[ii] = np.log2(Lx / dx)
        rect_pts = np.vstack([panel_xs[pvert],panel_vs[pvert]]).T
        patches.append(Polygon(rect_pts))

    p = PatchCollection(patches, cmap=mymap)
    p.set_array(panels_fs)
    p.set_clim(height_range[0] - 0.5, height_range[1] + 0.5)
    ax.add_collection(p)
    plt_ticks = np.arange(height_range[0],height_range[1]+1)
    cb = fig.colorbar(p, ax=ax, ticks=plt_ticks)
    # cb.set_label('panel height')
    # ax.set_title(f't={iter_num*sd["dt"]:.3f}')

    plt.tight_layout()
    plt.savefig(sim_dir_str + f'panel_height_image_t_{iter_num*sd["dt"]}.png',dpi=FIG_DPI)
    plt.savefig(sim_dir_str + f'panel_height_image_t_{iter_num*sd["dt"]}.pdf')
    print('panel height plot done!')
    t2 = time.time()
    print(f'panels plot took {t2-t1:.3f}s')
    plt.close()
# end plot height

def plot_e(sim_dir, simulation_dictionary, step_ii,flim, simulation_has_run = True, do_save = False):
    sd = simulation_dictionary
    output_dir = 'simulation_output/'
    sim_dir_str = ''
    if sim_dir is not None:
        sim_dir_str = sim_dir
        if sim_dir[-1] != '/':
            sim_dir_str += '/'
        output_dir = sim_dir_str + 'simulation_output/'

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
    ax1.set_ylim(sd['vmin'],sd['vmax'])

    ax1.grid()
    plt.tight_layout()
    if do_save:
        plt.savefig(sim_dir_str + f'e_field_time_{step_ii}.png',dpi=FIG_DPI)
        plt.savefig(sim_dir_str + f'e_field_time_{step_ii}.pdf')
    plt.close()
# end plot_e

    # %%time
# do_show_panels = False
def phase_movie(sim_dir, simulation_dictionary,do_show_panels,flim=(0,.3), simulation_has_run=True, can_do_movie=True):
    sd = simulation_dictionary
    output_dir = 'simulation_output/'
    sim_dir_str = ''
    if sim_dir is not None:
        sim_dir_str = sim_dir
        if sim_dir[-1] != '/':
            sim_dir_str += '/'
        output_dir = sim_dir_str + 'simulation_output/'
    print('starting phase space movie')
    t1 =time.time()
    if not simulation_has_run:
        print('unable to plot; simulation has not run or had errors')
        return
    if not can_do_movie:
        print('unable to load ffmpeg writer, movie capability inaccessible')
        return
    panel_string = ''

    if do_show_panels:
        panel_string = 'show_panels_'

    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='phase space', artist='Matplotlib',
                    comment='')
    writer = FFMpegWriter(fps=5, metadata=metadata)

    fig, ax = plt.subplots(figsize=(8,6))

    with writer.saving(fig, sim_dir_str + panel_string + 'phase_space'+ ".mp4", dpi=FIG_DPI):


#     xs = np.fromfile(output_dir + f'xs/xs_{step_ii}')
#     vs = np.fromfile(output_dir + f'vs/vs_{step_ii}')
#     fs = np.fromfile(output_dir + f'fs/fs_{step_ii}')
#     es = np.fromfile(output_dir + f'es/es_{step_ii}')
#     panels = np.fromfile(output_dir + f'panels/leaf_point_inds_{step_ii}',dtype='int32')
#     num_panels = int(panels.size/9)
#     panels = np.reshape(panels, (num_panels,9))
#     panels_fs = np.zeros(4*num_panels)
#         flim = (0,.3)

        panels = np.fromfile(output_dir + 'panels/leaf_point_inds_0',dtype='int32')
        num_panels = int(panels.size/9)
        panels = np.reshape(panels, (num_panels,9))
        panels_fs = np.zeros(4*num_panels)
        xs = np.fromfile(output_dir + 'xs/xs_0')
        vs = np.fromfile(output_dir + 'vs/vs_0')
        fs = np.fromfile(output_dir + 'fs/fs_0')

        patches = []
        for ii, panel in enumerate(panels):
            panel_xs = xs[panel]
            panel_vs = vs[panel]
            panel_fs = fs[panel]
#             weights = np.array([1,4,1,4,16,4,1,4,1])
#             panels_fs[ii] = 1./36. * np.dot(panel_fs,weights)
            
            p0 = [0,1,4,3]
            panels_fs[4*ii] = .25*sum(panel_fs[p0])
            rect_pts = np.vstack([panel_xs[p0],panel_vs[p0]]).T
            patches.append(Polygon(rect_pts))
            
            p1 = [1,2,5,4]
            panels_fs[4*ii+1] = .25*sum(panel_fs[p1])
            rect_pts = np.vstack([panel_xs[p1],panel_vs[p1]]).T
            patches.append(Polygon(rect_pts))
            
            p2 = [3,4,7,6]
            panels_fs[4*ii+2] = .25*sum(panel_fs[p2])
            rect_pts = np.vstack([panel_xs[p2],panel_vs[p2]]).T
            patches.append(Polygon(rect_pts))
            
            p3 = [4,5,8,7]
            panels_fs[4*ii+3] = .25*sum(panel_fs[p3])
            rect_pts = np.vstack([panel_xs[p3],panel_vs[p3]]).T
            patches.append(Polygon(rect_pts))

        p = PatchCollection(patches, cmap=matplotlib.cm.jet)
        if do_show_panels:
            p.set_ec('black')
        p.set_array(panels_fs)
        p.set_clim(flim)
        ax.add_collection(p)
        cb = fig.colorbar(p, ax=ax)
        cb.set_label('f')

        ax.set_xlim([sd['xmin'], sd['xmax']])
        ax.set_ylim([sd['vmin'],sd['vmax']])
        ax.set_xlabel('x')
        ax.set_ylabel('v')
        ax.set_title(f't=0.000')

        fig.canvas.draw()
        writer.grab_frame()


        print_update_frequency = int(np.ceil((sd['num_steps']+1)/6/sd['diag_period']))
        print_update_counter = 0

        for iter_num in range(sd['diag_period'], sd['num_steps'] + 1, sd['diag_period']):
            cb.remove()
            ax.collections.pop()

            panels = np.fromfile(output_dir + f'panels/leaf_point_inds_{iter_num}',dtype='int32')
            num_panels = int(panels.size/9)
            panels = np.reshape(panels, (num_panels,9))
            panels_fs = np.zeros(num_panels*4)
            xs = np.fromfile(output_dir + f'xs/xs_{iter_num}')
            vs = np.fromfile(output_dir + f'vs/vs_{iter_num}')
            fs = np.fromfile(output_dir + f'fs/fs_{iter_num}')

            patches = []
            for ii, panel in enumerate(panels):
                panel_xs = xs[panel]
                panel_vs = vs[panel]
                panel_fs = fs[panel]
                
                p0 = [0,1,4,3]
                panels_fs[4*ii] = .25*sum(panel_fs[p0])
                rect_pts = np.vstack([panel_xs[p0],panel_vs[p0]]).T
                patches.append(Polygon(rect_pts))

                p1 = [1,2,5,4]
                panels_fs[4*ii+1] = .25*sum(panel_fs[p1])
                rect_pts = np.vstack([panel_xs[p1],panel_vs[p1]]).T
                patches.append(Polygon(rect_pts))

                p2 = [3,4,7,6]
                panels_fs[4*ii+2] = .25*sum(panel_fs[p2])
                rect_pts = np.vstack([panel_xs[p2],panel_vs[p2]]).T
                patches.append(Polygon(rect_pts))

                p3 = [4,5,8,7]
                panels_fs[4*ii+3] = .25*sum(panel_fs[p3])
                rect_pts = np.vstack([panel_xs[p3],panel_vs[p3]]).T
                patches.append(Polygon(rect_pts))
                

            p = PatchCollection(patches, cmap=matplotlib.cm.jet)
            if do_show_panels:
                p.set_ec('black')
            p.set_array(panels_fs)
            p.set_clim(flim)
            ax.add_collection(p)
            cb = fig.colorbar(p, ax=ax)
            cb.set_label('f')
            ax.set_title(f't={iter_num*sd["dt"]:.3f}')


            fig.canvas.draw()
            writer.grab_frame()


            if print_update_counter == print_update_frequency:
                print(f'Movie is about {iter_num/(sd["num_steps"]+1)*100 :0.0f}% complete')
                print_update_counter = 0
#                 plt.savefig(mya.sim_dir + f'phase_space_image_t_{iter_num*dt}.svg')
                plt.savefig(sim_dir_str + panel_string + f'phase_space_image_t_{iter_num*sd["dt"]}.png',dpi=FIG_DPI)
                plt.savefig(sim_dir_str + panel_string + f'phase_space_image_t_{iter_num*sd["dt"]}.pdf')
            print_update_counter += 1

    plt.savefig(sim_dir_str + panel_string + f'phase_space_image_t_{iter_num*sd["dt"]}.png',dpi=FIG_DPI)
    plt.savefig(sim_dir_str + panel_string + f'phase_space_image_t_{iter_num*sd["dt"]}.pdf')
    print('phase space movie done!')
    t2 = time.time()
    print(f'phase space movie took {t2-t1:.3f}s')
    plt.close()
# end phase space movie

def phase_movie_standard_tree(simulation_dictionary,root_dir=None,do_show_panels=False,flim=(0,.3), simulation_has_run=True, can_do_movie=True):

    sim_dir, directories_found = generate_standard_names_dirs(simulation_dictionary,root_dir)
    phase_movie(sim_dir, simulation_dictionary, do_show_panels=do_show_panels, flim=flim, simulation_has_run=simulation_has_run, can_do_movie=can_do_movie)
# logf movie
# %%time
def logf_movie(sim_dir, simulation_dictionary, simulation_has_run = True, can_do_movie = True, flim=(-12,0)):
    sd = simulation_dictionary
    output_dir = 'simulation_output/'
    sim_dir_str = ''
    if sim_dir is not None:
        sim_dir_str = sim_dir
        if sim_dir[-1] != '/':
            sim_dir_str += '/'
        output_dir = sim_dir + 'simulation_output/'

    print('starting logf movie')
    t1 = time.time()
    if not simulation_has_run:
        print('unable to plot; simulation has not run or had errors')
        return
    if not can_do_movie:
        print('unable to load ffmpeg writer, movie capability inaccessible')
        return
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='logf movie', artist='Matplotlib',
                    comment='')
    writer = FFMpegWriter(fps=5, metadata=metadata)

    fig, ax = plt.subplots(figsize=(8,6))

    with writer.saving(fig, sim_dir_str+ 'logf_movie.mp4', dpi=FIG_DPI):


    #     flim = (-8,0)

        panels = np.fromfile(output_dir + 'panels/leaf_point_inds_0',dtype='int32')
        num_panels = int(panels.size/9)
        panels = np.reshape(panels, (num_panels,9))
        panels_fs = np.zeros(4*num_panels)
        xs = np.fromfile(output_dir + 'xs/xs_0')
        vs = np.fromfile(output_dir + 'vs/vs_0')
        fs = np.fromfile(output_dir + 'fs/fs_0')

        patches = []
        for ii, panel in enumerate(panels):
            panel_xs = xs[panel]
            panel_vs = vs[panel]
            panel_fs = fs[panel]

            p0 = [0,1,4,3]
            panels_fs[4*ii] = .25*sum(panel_fs[p0])
            rect_pts = np.vstack([panel_xs[p0],panel_vs[p0]]).T
            patches.append(Polygon(rect_pts))

            p1 = [1,2,5,4]
            panels_fs[4*ii+1] = .25*sum(panel_fs[p1])
            rect_pts = np.vstack([panel_xs[p1],panel_vs[p1]]).T
            patches.append(Polygon(rect_pts))

            p2 = [3,4,7,6]
            panels_fs[4*ii+2] = .25*sum(panel_fs[p2])
            rect_pts = np.vstack([panel_xs[p2],panel_vs[p2]]).T
            patches.append(Polygon(rect_pts))

            p3 = [4,5,8,7]
            panels_fs[4*ii+3] = .25*sum(panel_fs[p3])
            rect_pts = np.vstack([panel_xs[p3],panel_vs[p3]]).T
            patches.append(Polygon(rect_pts))

        p = PatchCollection(patches, cmap=matplotlib.cm.jet)

        p.set_array(np.log10(panels_fs))
        p.set_clim(flim)
        ax.add_collection(p)
        cb = fig.colorbar(p, ax=ax)
        cb.set_label('log f')

        ax.set_xlim(sd["xmin"], sd["xmax"])
        ax.set_ylim(sd["vmin"], sd["vmax"])
        ax.set_xlabel('x')
        ax.set_ylabel('v')
        ax.set_title(f't=0.000')

        fig.canvas.draw()
        writer.grab_frame()


        print_update_frequency = int(np.ceil((sd["num_steps"]+1)/6/sd["diag_period"]))
        print_update_counter = 0

        for iter_num in range(sd["diag_period"], sd["num_steps"] + 1, sd["diag_period"]):
            cb.remove()
            ax.collections.pop()

            panels = np.fromfile(output_dir + f'panels/leaf_point_inds_{iter_num}',dtype='int32')
            num_panels = int(panels.size/9)
            panels = np.reshape(panels, (num_panels,9))
            panels_fs = np.zeros(4*num_panels)
            xs = np.fromfile(output_dir + f'xs/xs_{iter_num}')
            vs = np.fromfile(output_dir + f'vs/vs_{iter_num}')
            fs = np.fromfile(output_dir + f'fs/fs_{iter_num}')

            patches = []
            for ii, panel in enumerate(panels):
                panel_xs = xs[panel]
                panel_vs = vs[panel]
                panel_fs = fs[panel]

                p0 = [0,1,4,3]
                panels_fs[4*ii] = .25*sum(panel_fs[p0])
                rect_pts = np.vstack([panel_xs[p0],panel_vs[p0]]).T
                patches.append(Polygon(rect_pts))

                p1 = [1,2,5,4]
                panels_fs[4*ii+1] = .25*sum(panel_fs[p1])
                rect_pts = np.vstack([panel_xs[p1],panel_vs[p1]]).T
                patches.append(Polygon(rect_pts))

                p2 = [3,4,7,6]
                panels_fs[4*ii+2] = .25*sum(panel_fs[p2])
                rect_pts = np.vstack([panel_xs[p2],panel_vs[p2]]).T
                patches.append(Polygon(rect_pts))

                p3 = [4,5,8,7]
                panels_fs[4*ii+3] = .25*sum(panel_fs[p3])
                rect_pts = np.vstack([panel_xs[p3],panel_vs[p3]]).T
                patches.append(Polygon(rect_pts))

    #             panel_vs = vs[panel]
    #             panel_fs = fs[panel]
    #             panels_fs[ii] = np.log10(.25*sum(panel_fs))
    #             rect_pts = np.vstack([panel_xs,panel_vs]).T
    #             patches.append(Polygon(rect_pts))

            p = PatchCollection(patches, cmap=matplotlib.cm.jet)

            p.set_array(np.log10(panels_fs))
            p.set_clim(flim)
            ax.add_collection(p)
            cb = fig.colorbar(p, ax=ax)
            cb.set_label('log f')
            ax.set_title(f't={iter_num*sd["dt"]:.3f}')


            fig.canvas.draw()
            writer.grab_frame()


            if print_update_counter == print_update_frequency:
                print(f'Movie is about {iter_num/(sd["num_steps"]+1)*100 :0.0f}% complete')
                print_update_counter = 0
    #             plt.savefig(mya.sim_dir + f'logf_phase_space_image_t_{iter_num*dt}.svg')
                plt.savefig(sim_dir_str + f'logf_phase_space_image_t_{iter_num*sd["dt"]}.png',dpi=FIG_DPI)
                plt.savefig(sim_dir_str + f'logf_phase_space_image_t_{iter_num*sd["dt"]}.pdf')
            print_update_counter += 1

        plt.savefig(sim_dir_str + f'logf_phase_space_final_image_t_{iter_num*sd["dt"]}.pdf')
        plt.savefig(sim_dir_str + f'logf_phase_space_final_image_t_{iter_num*sd["dt"]}.png',dpi=FIG_DPI)
    print('done with logf movie!')
    t2 = time.time()
    print(f'logf movie took {t2 - t1:.3f}s')
    plt.close()
# end logf movie

def logf_movie_standard_tree(simulation_dictionary,root_dir=None,simulation_has_run = True, can_do_movie = True, flim=(-8,0)):

    sim_dir, directories_found = generate_standard_names_dirs(simulation_dictionary,root_dir)
    logf_movie(sim_dir, simulation_dictionary)

def panel_height_movie(sim_dir, simulation_dictionary, height_range = [7,11],simulation_has_run=True, can_do_movie=True):
    sd = simulation_dictionary
    Lx = sd['xmax'] - sd['xmin']
    output_dir = 'simulation_output/'
    sim_dir_str = ''
    if sim_dir is not None:
        sim_dir_str = sim_dir
        if sim_dir[-1] != '/':
            sim_dir_str += '/'
        output_dir = sim_dir_str + 'simulation_output/'
    print('starting panel height movie')
    t1 =time.time()
    if not simulation_has_run:
        print('unable to plot; simulation has not run or had errors')
        return
    if not can_do_movie:
        print('unable to load ffmpeg writer, movie capability inaccessible')
        return

    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='panel heights', artist='Matplotlib',
                    comment='')
    writer = FFMpegWriter(fps=5, metadata=metadata)

    fig, ax = plt.subplots(figsize=(8,6))

    with writer.saving(fig, sim_dir_str + 'panel_heights'+ ".mp4", dpi=FIG_DPI):


        panels = np.fromfile(output_dir + 'panels/leaf_point_inds_0',dtype='int32')
        num_panels = int(panels.size/9)
        panels = np.reshape(panels, (num_panels,9))
        panels_fs = np.zeros(num_panels)
        xs = np.fromfile(output_dir + 'xs/xs_0')
        vs = np.fromfile(output_dir + 'vs/vs_0')

        patches = []
        pvert = [0,2,8,6]
        for ii, panel in enumerate(panels):
            panel_xs = xs[panel]
            panel_vs = vs[panel]
            dx = panel_xs[8] - panel_xs[0]
            
            panels_fs[ii] = np.log2(Lx / dx)
            rect_pts = np.vstack([panel_xs[pvert],panel_vs[pvert]]).T
            patches.append(Polygon(rect_pts))
            
        ncolors = height_range[1] - height_range[0] + 1
        # mymap = cm.get_cmap('gist_rainbow_r',ncolors)
        mymap = cm.get_cmap(panels_cmap,ncolors)
        p = PatchCollection(patches, mymap)
        p.set_array(panels_fs)
        p.set_clim(height_range[0] - 0.5, height_range[1] + 0.5)
        ax.add_collection(p)
        # cb = fig.colorbar(p, ax=ax)
        plt_ticks = np.arange(height_range[0],height_range[1]+1)
        cb = fig.colorbar(p, ax=ax, ticks=plt_ticks)
        cb.set_label('panel x-level')

        ax.set_xlim([sd['xmin'], sd['xmax']])
        ax.set_ylim([sd['vmin'],sd['vmax']])
        ax.set_xlabel('x')
        ax.set_ylabel('v')
        ax.set_title(f't=0.000')

        fig.canvas.draw()
        writer.grab_frame()


        print_update_frequency = int(np.ceil((sd['num_steps']+1)/6/sd['diag_period']))
        print_update_counter = 0

        for iter_num in range(sd['diag_period'], sd['num_steps'] + 1, sd['diag_period']):
            cb.remove()
            ax.collections.pop()

            panels = np.fromfile(output_dir + f'panels/leaf_point_inds_{iter_num}',dtype='int32')
            num_panels = int(panels.size/9)
            panels = np.reshape(panels, (num_panels,9))
            panels_fs = np.zeros(num_panels)
            xs = np.fromfile(output_dir + f'xs/xs_{iter_num}')
            vs = np.fromfile(output_dir + f'vs/vs_{iter_num}')

            patches = []
            for ii, panel in enumerate(panels):
                panel_xs = xs[panel]
                panel_vs = vs[panel]
                dx = panel_xs[8] - panel_xs[0]
                
                panels_fs[ii] = np.log2(Lx / dx)
                rect_pts = np.vstack([panel_xs[pvert],panel_vs[pvert]]).T
                patches.append(Polygon(rect_pts))

            p = PatchCollection(patches, cmap=mymap)
            p.set_array(panels_fs)
            p.set_clim(height_range[0] - 0.5, height_range[1] + 0.5)
            ax.add_collection(p)
            cb = fig.colorbar(p, ax=ax, ticks=plt_ticks)
            cb.set_label('panel x-level')
            ax.set_title(f't={iter_num*sd["dt"]:.3f}')


            fig.canvas.draw()
            writer.grab_frame()


            if print_update_counter == print_update_frequency:
                print(f'Movie is about {iter_num/(sd["num_steps"]+1)*100 :0.0f}% complete')
                print_update_counter = 0
#                 plt.savefig(mya.sim_dir + f'phase_space_image_t_{iter_num*dt}.svg')
                plt.savefig(sim_dir_str + f'panel_height_image_t_{iter_num*sd["dt"]}.png',dpi=FIG_DPI)
                plt.savefig(sim_dir_str + f'panel_height_image_t_{iter_num*sd["dt"]}.pdf')
            print_update_counter += 1

    plt.savefig(sim_dir_str + f'panel_height_image_t_{iter_num*sd["dt"]}.pdf')
    plt.savefig(sim_dir_str + f'panel_height_image_t_{iter_num*sd["dt"]}.png',dpi=FIG_DPI)
    print('panel height movie done!')
    t2 = time.time()
    print(f'panels movie took {t2-t1:.3f}s')
    plt.close()
# end panels movie

def panel_height_movie_standard_tree(simulation_dictionary,root_dir=None,height_range=(7,11), simulation_has_run=True, can_do_movie=True):

    sim_dir, directories_found = generate_standard_names_dirs(simulation_dictionary,root_dir)
    panel_height_movie(sim_dir, simulation_dictionary, height_range=height_range, simulation_has_run=simulation_has_run, can_do_movie=can_do_movie)


def sim_diagnostics_sample(simulation_dictionary, sim_dir = None, test_times=[45.0]):
    """Generate and plot diagnostics 

    This is a convenience function.  Diagnostics will need to be developed on a project-by-project basis.


    Parameters
    ----------
    sim_dir : string, optional
        where to find simulation output, by default None
    """
    t1 = time.time()
    sim_dir_str = ''
    if sim_dir is not None:
        sim_dir_str = sim_dir
        if sim_dir[-1] != '/':
            sim_dir_str += '/'
    sd = simulation_dictionary
    num_steps = sd['num_steps']
    diag_freq = sd['diag_period']
    dt = sd['dt']
    try:
        q = sd['q']
        m = sd['m']
        qm = sd['qm']
    except KeyError:
        print('q, m, qm not found in simulation dictionary')
        print ('using default values of -1, 1, -1, respectively')
        q = -1
        m = 1
        qm = q / m
    Lx = sd['xmax'] - sd['xmin']
    kx = 2*np.pi / Lx
    simtype = sd['sim_type']
    xmin = sd['xmin']
    xmax = sd['xmax']
    vmin = sd['vmin']
    vmax = sd['vmax']
    # tf = dt * num_steps

    output_dir = sim_dir_str + 'simulation_output/'
    # file_list = os.listdir('simulation_output/fs')
    # iterations = []
    # for file_name in file_list:
    #     iterations.append(int(file_name[3:]))
    # iterations = np.sort(iterations)
    # diag_times = dt * iterations

    diag_times = np.arange(0,num_steps+1,diag_freq) * dt
    diag_times.tofile(output_dir + 'diag_times')
    tf = diag_times[-1]
    iterations = (diag_times / dt).astype('int')
    iterations.tofile(output_dir + 'iterations')
    total_charge = np.zeros_like(diag_times)
    total_momentum = np.zeros_like(diag_times)
    total_kinetic = np.zeros_like(diag_times)
    total_potential = np.zeros_like(diag_times)
    total_entropy = np.zeros_like(diag_times)
    total_negative_area = np.zeros_like(diag_times)
    num_negative = np.zeros_like(diag_times)
    num_points = np.zeros_like(diag_times)
    frac_negative = np.zeros_like(diag_times)
    frac_negative_mass = np.zeros_like(diag_times) 

    min_f = np.zeros_like(diag_times)
    max_f = np.zeros_like(diag_times)

    l1f = np.zeros_like(diag_times)
    l2f = np.zeros_like(diag_times)

    # Lv = FS.vmin - FS.vmax

    for ii, t in enumerate(diag_times):
        iter_num = ii * diag_freq
        simtime = iter_num * dt

        xs = np.fromfile(output_dir + f'xs/xs_{iter_num}')
        vs = np.fromfile(output_dir + f'vs/vs_{iter_num}')
        fs = np.fromfile(output_dir + f'fs/fs_{iter_num}')
        qws = np.fromfile(output_dir + f'qws/qws_{iter_num}')
        es = np.fromfile(output_dir + f'es/es_{iter_num}')
        panel_point_inds = np.fromfile(output_dir + f'panels/leaf_point_inds_{iter_num}',dtype='int32')
        pinds = np.reshape(panel_point_inds, (int(panel_point_inds.size/9),9))
        
        num_points[ii] = xs.size

        total_charge[ii] = np.sum(qws)
        total_momentum[ii] = 1/qm * np.dot(vs, qws)
        total_kinetic[ii] = 1/qm * np.dot(vs **2, qws)
        
        # potential
        x_sort, sort_inds = np.unique(xs, return_index=True)
        E_sort = es[sort_inds]
        mod_x_inds = np.nonzero((x_sort >= 0) * (x_sort < Lx))
        modx = x_sort[mod_x_inds]
        modE = E_sort[mod_x_inds]
        dxs = .5 * np.hstack([modx[1]+Lx-modx[-1], modx[2:] - modx[:-2], modx[0]+Lx - modx[-2]])
        if np.nonzero(dxs<=0)[0].size > 0:
            raise ValueError
        total_potential[ii] = np.dot(dxs, modE**2)
        # entropy
        where_negative = np.nonzero(fs < 0)
        total_negative_area[ii] = 1/q * np.sum(qws[where_negative])
        num_negative[ii] = where_negative[0].size
        frac_negative[ii] = num_negative[ii] / fs.size
        frac_negative_mass[ii] = abs(np.sum(qws[where_negative]) / np.sum(abs(qws)) ) 
        where_positive = np.nonzero(fs > 0)
        total_entropy[ii] = 1/q * np.dot(fs[where_positive] * np.log10(fs[where_positive]), qws[where_positive])
        
        min_f[ii] = np.min(fs)
        max_f[ii] = np.max(fs)
        l1f[ii] = 1/abs(q) * np.sum(abs(qws))
        l2f[ii] = 1/q * np.dot(qws,fs)

    # write diagnostic arrays to file
    num_points.tofile(output_dir + 'num_points')
    total_charge.tofile(output_dir + 'total_charge')
    total_momentum.tofile(output_dir + 'total_momentum')
    total_kinetic.tofile(output_dir + 'total_kinetic')
    total_potential.tofile(output_dir + 'total_potential')
    frac_negative.tofile(output_dir + 'frac_negative')
    frac_negative_mass.tofile(output_dir + 'frac_negative_mass')
    total_entropy.tofile(output_dir + 'total_entropy')
    min_f.tofile(output_dir + 'min_f')
    max_f.tofile(output_dir + 'max_f')
    l1f.tofile(output_dir + 'l1f')
    l2f.tofile(output_dir + 'l2f')

    # plot diagnostics
    plt.rcParams['font.size']=18
    diag_fig_size=(8,6)
    # L2 E diagnostic
    # plt.figure()
    plt.figure(figsize=diag_fig_size)
    # plt.title(r'$||E||_2$')
    l2e = np.sqrt(total_potential)
    plt.semilogy(diag_times, l2e)

    #weak LD:
    # if type_str == 'weak_LD':
    if simtype == 1:
        g_t = .1533
        peak_ind = 0
        plt.semilogy(diag_times, l2e[peak_ind]*np.exp(-g_t*(diag_times - diag_times[peak_ind])), '-.')

    #strong LD:
    # if type_str == 'strong_LD':
    if simtype == 2:
    # if FS.simtype is SimType.STRONG_LD:
        g_strong1 = .2920
        g_strong2 = .0815
        damping1_ind = int(0.5/dt * 5 / diag_freq)
        damping1_times = diag_times[:int(0.5/dt * 25 / diag_freq)]
        grow1_start_ind = int(0.5/dt * 28 / diag_freq)
        grow1_times = diag_times[grow1_start_ind:int(0.5/dt * 90 / diag_freq)]
        grow1_ind = int(0.5/dt * 19 / diag_freq)
        grow1_ind_e = grow1_start_ind + grow1_ind
        plt.semilogy(damping1_times,l2e[damping1_ind]*np.exp(- g_strong1 * (damping1_times - damping1_times[damping1_ind])),'-.',label=r'$\sim e^{- %.05f t}$'%g_strong1)
        plt.semilogy(grow1_times,l2e[grow1_ind_e]*np.exp(g_strong2*(grow1_times-grow1_times[grow1_ind])),'-.k')

    plt.grid()
    #weak LD
    # if type_str == 'weak_LD':
    if simtype == 1:
    # if FS.simtype is SimType.WEAK_LD:
    #     pass
        tf = sd['num_steps'] * sd['dt']
        if tf > 85:
            plt.ylim([1e-16,1e-1])
        else:
            plt.ylim([1e-8,1e-1])
    #strong LD
    # if type_str == 'strong_LD' or type_str == 'strong_2str':
    if simtype in [2,3,4]:
        plt.ylim([1e-3, 1e1])
    plt.xlim(0,tf)
    plt.xlabel('t')
    plt.ylabel(r'$||E||_2$')
    plt.tight_layout()
    plt.savefig(sim_dir_str + 'l2E.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'l2E.pdf')
    plt.close()
    # end l2e diagnostic
    #---------------------------------

    # plt.figure()
    plt.figure(figsize=diag_fig_size)
    # plt.title(r'maximum of f')
    plt.xlabel('t')
    plt.plot(diag_times, max_f )
    plt.grid()
    plt.xlim(0,tf)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.savefig(sim_dir_str + 'maxf_conservation.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'maxf_conservation.pdf')
    plt.close()
    #---------------------------------
    plt.figure(figsize=diag_fig_size)
    # plt.figure()
    # plt.title(r'relative variation in maximum of f')
    plt.xlabel('t')
    plt.plot(diag_times, (max_f - max_f[0])/max_f[0])
    plt.grid()
    plt.xlim(0,tf)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tight_layout()
    plt.savefig(sim_dir_str + 'relative_maxf_conservation.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'relative_maxf_conservation.pdf')
    plt.close()
    #---------------------------------
    plt.figure(figsize=diag_fig_size)
    # plt.figure()
    # plt.title(r'minimum of f')
    plt.xlabel('t')
    plt.plot(diag_times, min_f)
    plt.grid()
    plt.xlim(0,tf)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tight_layout()
    plt.savefig(sim_dir_str + 'minf_conservation.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'minf_conservation.pdf')
    plt.close()
    #---------------------------------
    plt.figure(figsize=diag_fig_size)
    # plt.figure()
    # plt.title(r'relative variation in minimum of f')
    plt.xlabel('t')
    plt.plot(diag_times, (min_f - min_f[0])/max_f[0])
    plt.grid()
    plt.xlim(0,tf)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tight_layout()
    plt.savefig(sim_dir_str + 'relative_minf_conservation.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'relative_minf_conservation.pdf')
    plt.close()
    #---------------------------------
    plt.figure(figsize=diag_fig_size)
    # plt.figure()
    # plt.title(r'fraction of negative f values')
    plt.xlabel('t')
    plt.xlim(0,tf)
    plt.plot(diag_times, frac_negative)
    plt.grid()
    plt.tight_layout()
    plt.savefig(sim_dir_str + 'frac_negatives.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'frac_negatives.pdf')
    plt.close()
    #---------------------------------
    plt.figure(figsize=diag_fig_size)
    # plt.figure()
    # plt.title(r'fraction of negative mass')
    plt.xlabel('t')
    plt.xlim(0,tf)
    plt.plot(diag_times, frac_negative_mass)
    plt.grid()
    plt.tight_layout()
    plt.savefig(sim_dir_str + 'frac_negative_mass.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'frac_negative_mass.pdf')
    plt.close()
    #---------------------------------

    plt.figure(figsize=diag_fig_size)
    # plt.figure()
    plt.plot(diag_times,num_points,'.',label='number of points')
    n_diags = sd['num_steps']+1
    initial_nx = 2**(sd['initial_height']+1) + 1
    max_nx = 2**(sd['max_height']+1) + 1
    plt.plot(diag_times, initial_nx**2 * np.ones_like(num_points),label=r'$%i^2$ points'%initial_nx)
    plt.plot(diag_times, max_nx**2 * np.ones_like(num_points),label=r'$%i^2$ points'%max_nx)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.grid()
    plt.xlabel('t')
    plt.xlim(0,tf)
    # plt.title('number of points in simulation')
    plt.legend()

    plt.tight_layout()
    plt.savefig(sim_dir_str + 'num_points.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'num_points.pdf')
    plt.close()
    #---------------------------------

    plt.figure(figsize=diag_fig_size)
    # plt.figure()
    # plt.title(r'relative variation in total charge')
    plt.xlabel('t')
    plt.xlim(0,tf)
    plt.plot(diag_times, (total_charge - total_charge[0])/total_charge[0])
    plt.grid()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tight_layout()
    plt.savefig(sim_dir_str + 'relative_charge_conservation.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'relative_charge_conservation.pdf')
    plt.close()
    #---------------------------------
    plt.figure(figsize=diag_fig_size)
    # plt.figure()
    # plt.title(r'variation in total momentum')
    plt.xlabel('t')
    plt.xlim(0,tf)
    plt.plot(diag_times, total_momentum - total_momentum[0])
    plt.grid()
    plt.tight_layout()
    plt.savefig(sim_dir_str + 'momentum_conservation.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'momentum_conservation.pdf')
    plt.close()
    #---------------------------------
    plt.figure(figsize=diag_fig_size)
    # plt.figure()
    # plt.title(r'relative variation in $L_1(f^n)$')
    plt.xlabel('t')
    plt.xlim(0,tf)
    plt.plot(diag_times, (l1f - l1f[0])/l1f[0])
    plt.grid()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tight_layout()
    plt.savefig(sim_dir_str + 'relative_l1f_conservation.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'relative_l1f_conservation.pdf')
    plt.close()
    #---------------------------------
    plt.figure(figsize=diag_fig_size)
    # plt.figure()
    # plt.title(r'relative variation in $L_2(f^n)$')
    plt.xlabel('t')
    plt.xlim(0,tf)
    plt.plot(diag_times, (l2f - l2f[0])/l2f[0])
    plt.grid()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tight_layout()
    plt.savefig(sim_dir_str + 'relative_l2f_conservation.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'relative_l2f_conservation.pdf')
    plt.close()
    #---------------------------------
    plt.figure(figsize=diag_fig_size)
    # plt.figure()
    # plt.title('relative energy variation')
    plt.xlabel('t')
    plt.xlim(0,tf)
    total_energy = total_potential + total_kinetic
    plt.plot(diag_times, (total_energy - total_energy[0])/total_energy[0])
    plt.grid()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tight_layout()
    plt.savefig(sim_dir_str + 'relative_energy_conservation.png',dpi=LOW_DPI)
    plt.savefig(sim_dir_str + 'relative_energy_conservation.pdf')
    plt.close()
    # done plotting
    t2 = time.time()
    print(f'Done plotting diagnostics')
    print(f'Diagnostic collection and plot time {t2-t1:.3f}s')
    #--------------- f(v) plots --------
    nx0 = 2**(sd["initial_height"]+1) + 1
    dx0 = (sd["xmax"] - sd["xmin"]) / nx0
    for test_time in test_times:
        if simtype in [SimType.WEAK_LD, SimType.STRONG_LD]:# create new diagnostic : f(x=a,t=b,v)
            try:
                step_ii = iterations[np.nonzero(diag_times <= test_time)[0][-1]]
                xs = np.fromfile(output_dir + f'xs/xs_{step_ii}')
                vs = np.fromfile(output_dir + f'vs/vs_{step_ii}')
                fs = np.fromfile(output_dir + f'fs/fs_{step_ii}')
                inds = np.lexsort([vs,xs])
                xtest = 4 * np.pi - 0.005
                inds_xtest = np.nonzero(abs(xs - xtest) <= dx0/4.0)
                vs_xtest = vs[inds_xtest]
                fs_xtest = fs[inds_xtest]
                sort_inds = np.argsort(vs_xtest)
                # plt.figure()
                plt.figure(figsize=diag_fig_size)
                ax = plt.gca()
                simtime = diag_times[step_ii]
                # plt.title(f'f(t={simtime:.2f}, x={xtest:.2f}, v)')
                plt.plot(vs_xtest[sort_inds],fs_xtest[sort_inds],'r')
                plt.xlabel('v')
                plt.grid()
                plt.xlim(4.8,6.3)
                plt.ylim(-2e-6,6e-6)
                ax.ticklabel_format(axis='y',style='sci',scilimits=(0,1))
                plt.savefig(sim_dir_str + f'fv_t_{simtime:.0f}.png',dpi=LOW_DPI)
                plt.savefig(sim_dir_str + f'fv_t_{simtime:.0f}.pdf')
                plt.close()
            except IndexError:
                print(f'time {test_time:.1f} is less than all simulation times')
        if simtype in [SimType.STRONG_TWO_STREAM]:
            # test_iter = int(test_time / dt)
            try:
                test_iter = iterations[np.nonzero(diag_times <= test_time)[0][-1]]

                xs = np.fromfile(output_dir + f'xs/xs_{test_iter}')
                vs = np.fromfile(output_dir + f'vs/vs_{test_iter}')
                fs = np.fromfile(output_dir + f'fs/fs_{test_iter}')
                inds = np.lexsort([vs,xs])
                xtest = 2 * np.pi
                inds_xtest = np.nonzero(abs(xs - xtest) <= dx0 / 4.0)
                vs_xtest = vs[inds_xtest]
                fs_xtest = fs[inds_xtest]
                sort_inds = np.argsort(vs_xtest)

                # plt.figure()
                plt.figure(figsize=diag_fig_size)
                simtime = diag_times[test_iter]
                # plt.title(f'f(t={simtime:.2f}, x={xtest:.2f}, v)')
                plt.plot(vs_xtest[sort_inds],fs_xtest[sort_inds],'r')
                plt.xlabel('v')
                plt.ylabel('f')
                plt.grid()
                plt.ylim(-0.05,0.45)
                plt.xlim(-6,6)
                plt.savefig(sim_dir_str + f'fv_t_{simtime:.1f}.png',dpi=LOW_DPI)
                plt.savefig(sim_dir_str + f'fv_t_{simtime:.1f}.pdf')
                plt.close()
            except IndexError:
                print(f'time {test_time:.1f} is less than all simulation times')


        


# end sim_diagnostics_sample
    
def diagnostics_sample_standard_tree(simulation_dictionary,root_dir=None):

    sim_dir, directories_found = generate_standard_names_dirs(simulation_dictionary,root_dir)
    sim_diagnostics_sample(simulation_dictionary, sim_dir=sim_dir)

    

# main
if __name__ == '__main__':
    """
    Running as a script
    Use cases
    * (not implemented) run and/or plot simulation in directory as specified by command-line arguments
    * (not implemented) construct simulation directory from file and copy deck from current directory to constructed directory
    * (not implemented) construct deck and directory from command-line arguments
    """
    import argparse

    parser = argparse.ArgumentParser(description='Python command-line interface for FARRSIGHT')
    parser.add_argument('--root_dir', '-rd', help='Root directory of FARRSIGHT simulations, default to current working directory (cwd)', default=None)
    parser.add_argument('--sim_dir', '-sd', help='Location within root_dir of this FARRSIGHT simulation, default to cwd', default=None)
    parser.add_argument('--deck_dir', '-dd', help='Location of FARRSIGHT input deck, default to sim_dir', default=None)
    parser.add_argument('--standard_tree', '-std', action='store_true', help='Get simulation directory from deck in deck directory (default cwd) + root directory')
    parser.add_argument('--deck_dir_cwd', '-cwd', action='store_true', help='force use of deck found in current working directory')
    parser.add_argument('--deck_name', '-dn', help='Name of input_deck.  Default is "deck.json"', default='deck.json')
    parser.add_argument('--make_deck', '-make', help='')
    parser.add_argument('--run','-r',action='store_true', help='Run the simulation using deck in deck directory, store output in sim_dir')
    parser.add_argument('--gpu',dest='use_gpu', action='store_true', help='use gpu')
    parser.add_argument('--phase_movie','-pha',action='store_true', help='make phase space movie from data in sim_dir')
    parser.add_argument('--show_panels','-p',action='store_true', help='show panels in phase space movie')
    parser.add_argument('--logf_plot','-logp',action='store_true',help='generate log f plots at times in plot_times argument')
    parser.add_argument('--logf_movie','-log',action='store_true', help='make log f movie')
    parser.add_argument('--panels_movie','-panel',action='store_true', help='make panel height movie')
    parser.add_argument('--plot_diagnostics','-pd',action='store_true',help='plot diagnostics')
    parser.add_argument('--plot_times','-pt', type=float, nargs='*', help='times at which to generate phase space plots')
    parser.add_argument('--dict_args','-args', help='optional edits to deck, of the form $ python FARRSIGHT_Interface.py ...args... name value type.' +\
                        ' For example,  dt 0.5 float',nargs='*')
    args = parser.parse_args()

    # convert other arguments to a dictionary

    dict_args = {}
    if args.dict_args is not None:
        num_dict_entries = len(args.dict_args)//3
        wrong_num_entries = len(args.dict_args) % 3
        if 'amr_epsilons' in args.dict_args:
            dict_args['amr_epsilons'] = []
        for ii in range(num_dict_entries):
            if args.dict_args[3*ii] == 'amr_epsilons':
                dict_args['amr_epsilons'].append(float(args.dict_args[3*ii+1]))
            elif (ii != num_dict_entries or wrong_num_entries == 0):
                if (args.dict_args[3*ii+2] == 'int'):
                    dict_args[args.dict_args[3*ii]] = int(args.dict_args[3*ii+1])
                elif (args.dict_args[3*ii+2] == 'float'):
                    dict_args[args.dict_args[3*ii]] = float(args.dict_args[3*ii+1])
                else:
                    dict_args[args.dict_args[3*ii]] = args.dict_args[3*ii+1]
        # if 'amr_epsilons' in dict_args:
        #     dict_args['amr_epsilons'] = [dict_args['amr_epsilons']]
    

    root_dir_str = ''
    if args.root_dir is not None:
        root_dir_str = args.root_dir
        if root_dir_str[-1] != '/':
            root_dir_str += '/'
        print('root directory ', root_dir_str)
    else:
        print('no root dir provided, working in ', os.getcwd() )
    simulation_dictionary = None
    sim_dir_str = ''
    dictionaries_found = False
    if args.standard_tree:
        deck_dir = None
        deck_str = 'deck.json'
        if args.deck_dir is None or args.deck_dir_cwd:
            deck_str = ''
            deck_dir = None
        else:
            deck_str = args.deck_dir
            if args.deck_dir[-1] != '/':
                deck_str += '/'
            deck_dir = deck_dir_str
        if args.deck_name is None:
            deck_str += 'deck.json'
        else:
            deck_str += args.deck_name
        # get sim_dir from deck
        # simulation_dictionary = deck_to_dict(deck_dir = deck_dir, deck_name=args.deck_name)
        simulation_dictionary = update_dictionary(deck_dir=deck_dir, deck_name = args.deck_name, **dict_args)
        sim_dir_str, directories_found = generate_standard_names_dirs(simulation_dictionary, root_dir=args.root_dir)
        sim_dir = sim_dir_str
        shutil.copy2(deck_str, sim_dir_str)
        # end standard tree from deck option
    else:
        sim_dir = None
        if args.sim_dir is not None:
            sim_dir_str = root_dir_str + args.sim_dir
            if sim_dir_str[-1] != '/':
                sim_dir_str += '/'
            sim_dir = sim_dir_str
            if not os.path.exists(sim_dir_str):
                os.makedirs(sim_dir_str)
        else:
            if args.root_dir is not None:
                sim_dir_str = root_dir_str
                sim_dir = sim_dir_str
                if not os.path.exists(sim_dir_str):
                    os.makedirs(sim_dir_str)

        deck_dir = None
        if not args.deck_dir_cwd:
            if args.deck_dir is None:
                deck_dir = sim_dir
            else:
                deck_dir = args.deck_dir
                if args.deck_dir[-1] != '/':
                    deck_dir += '/'
        # simulation_dictionary = deck_to_dict(deck_dir = deck_dir, deck_name=args.deck_name)
        simulation_dictionary = update_dictionary(deck_dir=deck_dir, deck_name = args.deck_name, **dict_args)
        
        dictionaries_found = True
        # end not-from-deck option
    sd = simulation_dictionary
    if sim_dir is not None:
        print('sim directory', sim_dir_str)

    if args.run:
        run_sim(sim_dir=sim_dir, deck_dir=deck_dir, deck_name=args.deck_name, use_gpu=args.use_gpu)

    if args.plot_diagnostics:
        iterations = np.arange(0,sd['num_steps']+1, sd['diag_period'])
        final_iter = iterations[-1]
        sim_diagnostics_sample(simulation_dictionary, sim_dir=sim_dir)
        if SimType(sd['sim_type']) is not SimType.FRIEDMAN_BEAM:
            # plot_phase_space(sim_dir, simulation_dictionary, int(45.0 / simulation_dictionary['dt']/simulation_dictionary['diag_period']), sim_type_to_flim[SimType(simulation_dictionary['sim_type'])])
            plot_phase_space(sim_dir, simulation_dictionary, final_iter, sim_type_to_flim[SimType(simulation_dictionary['sim_type'])])
        if args.show_panels or args.panels_movie:
            plot_height(sim_dir, sd, final_iter)

    if args.plot_times is not None:
        t1 = time.time()
        iterations = np.arange(0,sd['num_steps']+1,sd['diag_period'])
        diag_ts = iterations * sd['dt']
        for t0 in args.plot_times:
            print('generating time',t0,'phase space')
            t0_ind = np.nonzero(diag_ts <= t0)[0][-1]
            t0_iter = int(iterations[t0_ind])
            t0_sim = diag_ts[t0_ind]
            if SimType(sd['sim_type']) is SimType.COLDER_TWO_STREAM:
                flim = (0,0.5/np.sqrt(2*np.pi)/sd['vth'])
            else:
                flim = sim_type_to_flim[SimType(sd['sim_type'])]
            plot_phase_space(sim_dir,sd,t0_iter,flim)
            if args.logf_movie or args.logf_plot:
                plot_logf(sim_dir, sd, t0_iter)
                
            if sd['adaptively_refine'] or args.show_panels or args.panels_movie:
                hlim2 = sd['max_height']
                if 'v_height' in sd:
                    hlim2 -= sd['v_height']
                plot_height(sim_dir, sd, t0_iter, height_range=[sd['initial_height'],hlim2])
        t2 = time.time()
        print(f'spent {t2-t1:.1f}s making phase space (and log or height) images')

    if args.phase_movie:
        simtype = simulation_dictionary['sim_type']
        if simtype == 1:
            flim = (0,0.44)
        elif simtype == 2:
            flim = (0,0.47)
        elif simtype == 3:
            flim = (0,0.45)
        elif simtype == 4:
            flim = (0,0.8)
        do_show_panels = False
        phase_movie(sim_dir, simulation_dictionary, do_show_panels, flim=flim, can_do_movie=can_do_movie)
        if args.show_panels:
            do_show_panels = True
            phase_movie(sim_dir, simulation_dictionary, do_show_panels, flim=flim, can_do_movie=can_do_movie)

    if args.panels_movie and sd['adaptively_refine']==1:
       hlim2 = sd['max_height']
       if 'v_height' in sd:
           hlim2 -= sd['v_height']
       panel_height_movie(sim_dir, simulation_dictionary,\
           height_range=[simulation_dictionary['initial_height'],hlim2],\
           can_do_movie=can_do_movie)

    if args.logf_movie:
        logf_movie(sim_dir, simulation_dictionary, can_do_movie=can_do_movie)
    

    print('done with this run')

# %%

"""
Collection of functions for creating and running FARRSIGHT simulations

Contains
---
SimType : Enumerated type
make_dirs [DEPRECATED]:
generate_standard_names_dirs :
dict_to_deck :
deck_to_dict :
run_sim :
plot_phase_space :
phase_movie :
logf_movie :
sim_diagnostics_sample :

Dependencies
---
standard Python distribution (i.e. through Anaconda)

ffmpeg for movie writing

"""

import json # dump, load
import numpy as np
from matplotlib import pyplot as plt
import os # path.exists, makedirs
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
plt.rcParams.update({'font.size': 14})

from enum import IntEnum
class SimType(IntEnum):
    WEAK_LD = 1
    STRONG_LD = 2
    STRONG_TWO_STREAM = 3
    COLDER_TWO_STREAM = 4

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
        amr_string = 'max_height_%i'%sd['max_height']
        amr_string += '_amr_epsilons'
        for eps in sd['amr_epsilons']:
            amr_string += '_%.03f'%(eps)
    else:
        amr_string = 'no_amr'
#         sim_group = ''
#         sim_name = ''
    tf = sd['dt'] * sd['num_steps']
    physical_parameters = 'vth_%.3f_vstr_%.3f_amp_%.3f_normal_k_%.3f_tf_%.1f'%(sd['vth'], sd['vstr'], sd['amp'], sd['normalized_wavenumber'], tf)
    numerical_parameters = 'height0_%i_vm_%.1f_g_eps_%.3f_dt_%.3f_diag_freq_%i'%(sd['initial_height'], sd['vmax'], sd['greens_epsilon'], sd['dt'], sd['diag_period'])
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

def make_deck(deck_dir = None, deck_name = None, 
                from_existing = False, **sim_vars):
    simulation_dictionary = {}
    if from_existing:
        simulation_dictionary = deck_to_dict(deck_dir, deck_name)
    else:
        simulation_dictionary = \
        {'project_name':'template',\
        'xmin' : 0.0, 'xmax' : 4*np.pi,\
        'vmin' : -6.0, 'vmax' : 6.0,\
        'sim_type' : 2,\
        'normalized_wavenumber':1.0,\
        'amp':0.5, 'vth' : 1.0, 'vstr':0.0,\
        'initial_height' : 5,\
        'max_height' : 5,\
        'greens_epsilon' : 0.2,\
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
# end of make_deck

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
    
    
    fig,ax = plt.subplots()
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
    cb.set_label('f')

    ax0.set_xlim(sd['xmin'], sd['xmax'])
    ax0.set_ylim(sd['vmin'],sd['vmax'])
    ax0.set_ylabel('v')
    ax0.set_title(f'phase space with {num_panels} panels, t={simtime:.03f}')
    
plt.close()
# end plot_phase_space

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
    ax1.set_title(f'electric field, t={simtime:.03f}')
    ax1.plot(xs,es,'.')
    ax1.set_xlabel('x')

    ax1.set_xlim(sd['xmin'], sd['xmax'])
    ax1.set_ylim(sd['vmin'],sd['vmax'])

    ax1.grid()
    plt.tight_layout()
    if do_save:
        plt.savefig(sim_dir_str + f'e_field_time_{step_ii}.png')
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

    fig, ax = plt.subplots(figsize=(12,10))

    with writer.saving(fig, sim_dir_str + panel_string + 'phase_space'+ ".mp4", dpi=100):


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
        ax.set_title(f'phase space with {num_panels} panels, t=0.000')

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
            ax.set_title(f'phase space with {num_panels} panels, t={iter_num*sd["dt"]:.3f}')


            fig.canvas.draw()
            writer.grab_frame()


            if print_update_counter == print_update_frequency:
                print(f'Movie is about {iter_num/(sd["num_steps"]+1)*100 :0.0f}% complete')
                print_update_counter = 0
#                 plt.savefig(mya.sim_dir + f'phase_space_image_t_{iter_num*dt}.svg')
                plt.savefig(sim_dir_str + panel_string + f'phase_space_image_t_{iter_num*sd["dt"]}.png')
            print_update_counter += 1

    plt.savefig(sim_dir_str + panel_string + f'phase_space_image_t_{iter_num*sd["dt"]}.png')
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
def logf_movie(sim_dir, simulation_dictionary, simulation_has_run = True, can_do_movie = True, flim=(-8,0)):
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

    fig, ax = plt.subplots()

    with writer.saving(fig, sim_dir_str+ 'logf_movie.mp4', dpi=100):



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
        cb.set_label('f')

        ax.set_xlim(sd["xmin"], sd["xmax"])
        ax.set_ylim(sd["vmin"], sd["vmax"])
        ax.set_xlabel('x')
        ax.set_ylabel('v')
        ax.set_title(f'log(f), t=0.000')

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
            cb.set_label('f')
            ax.set_title(f'log(f), t={iter_num*sd["dt"]:.3f}')


            fig.canvas.draw()
            writer.grab_frame()


            if print_update_counter == print_update_frequency:
                print(f'Movie is about {iter_num/(sd["num_steps"]+1)*100 :0.0f}% complete')
                print_update_counter = 0
    #             plt.savefig(mya.sim_dir + f'logf_phase_space_image_t_{iter_num*dt}.svg')
                plt.savefig(sim_dir_str + f'logf_phase_space_image_t_{iter_num*sd["dt"]}.png')
            print_update_counter += 1

        plt.savefig(sim_dir_str + f'logf_phase_space_final_image_t_{iter_num*sd["dt"]}.png')
    print('done with logf movie!')
    t2 = time.time()
    print(f'logf movie took {t2 - t1:.3f}s')
    plt.close()
# end logf movie

def logf_movie_standard_tree(simulation_dictionary,root_dir=None,simulation_has_run = True, can_do_movie = True, flim=(-8,0)):

    sim_dir, directories_found = generate_standard_names_dirs(simulation_dictionary,root_dir)
    logf_movie(sim_dir, simulation_dictionary)

def sim_diagnostics_sample(simulation_dictionary, sim_dir = None):
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

    output_dir = sim_dir_str + 'simulation_output/'
    diag_times = np.arange(0,num_steps+1,diag_freq) * dt
    total_charge = np.zeros_like(diag_times)
    total_momentum = np.zeros_like(diag_times)
    total_kinetic = np.zeros_like(diag_times)
    total_potential = np.zeros_like(diag_times)
    total_entropy = np.zeros_like(diag_times)
    total_negative_area = np.zeros_like(diag_times)
    num_negative = np.zeros_like(diag_times)
    num_points = np.zeros_like(diag_times)
    frac_negative = np.zeros_like(diag_times) 

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
        where_positive = np.nonzero(fs > 0)
        total_entropy[ii] = 1/q * np.dot(fs[where_positive] * np.log10(fs[where_positive]), qws[where_positive])
        
        l1f[ii] = 1/q * np.sum(abs(qws))
        l2f[ii] = 1/q * np.dot(qws,fs)


    # plot diagnostics
    # L2 E diagnostic
    plt.figure()
    plt.title('l2 E')
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
        damping1_ind = int(0.5/dt * 5)
        damping1_times = diag_times[:int(0.5/dt * 25)]
        grow1_start_ind = int(0.5/dt * 28)
        grow1_times = diag_times[grow1_start_ind:int(0.5/dt * 90)]
        grow1_ind = int(0.5/dt * 19)
        grow1_ind_e = grow1_start_ind + grow1_ind
        plt.semilogy(damping1_times,l2e[damping1_ind]*np.exp(- g_strong1 * (damping1_times - damping1_times[damping1_ind])),'-.',label=r'$\sim e^{- %.05f t}$'%g_strong1)
        plt.semilogy(grow1_times,l2e[grow1_ind_e]*np.exp(g_strong2*(grow1_times-grow1_times[grow1_ind])),'-.k')

    plt.grid()
    #weak LD
    # if type_str == 'weak_LD':
    if simtype == 1:
    # if FS.simtype is SimType.WEAK_LD:
    #     pass
        plt.ylim([1e-8,1e-1])
    #strong LD
    # if type_str == 'strong_LD' or type_str == 'strong_2str':
    if simtype in [2,3,4]:
        plt.ylim([1e-3, 1e1])
    plt.savefig(sim_dir_str + 'l2E.png')
    plt.close()
    # end l2e diagnostic
    #---------------------------------

    plt.figure()
    plt.title(r'fraction of negative f values')
    plt.xlabel('t')
    plt.plot(diag_times, frac_negative)
    plt.grid()
    plt.savefig(sim_dir_str + 'frac_negatives.png')
    plt.close()
    #---------------------------------

    plt.figure()
    plt.plot(diag_times,num_points,'.',label='number of points')
    n_diags = sd['num_steps']+1
    initial_nx = 2**(sd['initial_height']+1) + 1
    max_nx = 2**(sd['max_height']+1) + 1
    plt.plot(diag_times, initial_nx**2 * np.ones_like(num_points),label=r'$%i^2$ points'%initial_nx)
    plt.plot(diag_times, max_nx**2 * np.ones_like(num_points),label=r'$%i^2$ points'%max_nx)
    plt.grid()
    plt.xlabel('t')
    plt.title('number of points in simulation')
    plt.legend()

    plt.savefig(sim_dir + 'num_points.png')
    plt.close()
    #---------------------------------

    plt.figure()
    plt.title(r'percent variation in total charge, $\sum_i q w^n_i - \sum_i q w^0_i$')
    plt.xlabel('t')
    plt.plot(diag_times, 100 * (total_charge - total_charge[0])/total_charge[0])
    plt.grid()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.savefig(sim_dir_str + 'percent_charge_conservation.png')
    plt.close()
    #---------------------------------
    plt.figure()
    plt.title(r'variation in total momentum, $\sum_i m v^n_i w^n_i - \sum_i m v^0_i w^0_i$')
    plt.xlabel('t')
    plt.plot(diag_times, total_momentum - total_momentum[0])
    plt.grid()
    plt.savefig(sim_dir_str + 'momentum_conservation.png')
    plt.close()
    #---------------------------------
    plt.figure()
    plt.title(r'percent variation in $L_2(f^n)$, $\sum_i (f^n_i)^2 w^n_i - \sum_i (f^n_i)^2 w^0_i$')
    plt.xlabel('t')
    plt.plot(diag_times, 100*(l2f - l2f[0])/l2f[0])
    plt.grid()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.savefig(sim_dir_str + 'percent_l2f_conservation.png')
    plt.close()
    #---------------------------------
    plt.figure()
    plt.title('percent energy variation')
    plt.xlabel('t')
    total_energy = total_potential + total_kinetic
    plt.plot(diag_times, 100 * (total_energy - total_energy[0])/total_energy[0])
    plt.grid()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.savefig(sim_dir_str + 'relative_energy_conservation.png')
    plt.close()
    # done plotting
    t2 = time.time()
    print(f'Done plotting diagnostics')
    print(f'Diagnostic collection and plot time {t2-t1:.3f}s')
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
    parser.add_argument('--run',action='store_true', help='Run the simulation using deck in deck directory, store output in sim_dir')
    parser.add_argument('--gpu',dest='use_gpu', action='store_true', help='boolean switch for using gpu')
    parser.add_argument('--phase_movie','-pha',action='store_true', help='use this flag to make phase space movie from data in sim_dir')
    parser.add_argument('--logf_movie','-log',action='store_true', help='use this flag to make log f movie')
    parser.add_argument('--plot_diagnostics','-pd',action='store_true',help='use this flag to plot diagnostics')
    args = parser.parse_args()

    root_dir_str = ''
    if args.root_dir is not None:
        root_dir_str = args.root_dir
        if root_dir_str[-1] != '/':
            root_dir_str += '/'
    simulation_dictionary = None
    sim_dir_str = ''
    dictionaries_found = False
    if args.standard_tree:
        deck_dir = None
        deck_dir_str = 'deck.json'
        if args.deck_dir is None or args.deck_dir_cwd:
            deck_dir_str = ''
            deck_dir = None
        else:
            deck_dir_str = args.deck_dir
            if args.deck_dir[-1] != '/':
                deck_dir_str += '/'
            deck_dir = deck_dir_str
        if args.deck_name is None:
            deck_dir_str += 'deck.json'
        else:
            deck_dir_str += args.deck_name
        # get sim_dir from deck
        simulation_dictionary = deck_to_dict(deck_dir = deck_dir, deck_name=args.deck_name)
        sim_dir_str, directories_found = generate_standard_names_dirs(simulation_dictionary, root_dir=args.root_dir)
        sim_dir = sim_dir_str
        shutil.copy2(deck_dir_str, sim_dir_str)
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
        simulation_dictionary = deck_to_dict(deck_dir = deck_dir, deck_name=args.deck_name)
        dictionaries_found = True
        # end not-from-deck option
    if args.run:
        run_sim(sim_dir=sim_dir, deck_dir=deck_dir, deck_name=args.deck_name, use_gpu=args.use_gpu)
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

    if args.logf_movie:
        logf_movie(sim_dir, simulation_dictionary, can_do_movie=can_do_movie)
    
    if args.plot_diagnostics:
        sim_diagnostics_sample(simulation_dictionary, sim_dir=sim_dir)


    # parser.add_argument('project_name')
    # parser.add_argument('xmin', type=float)
    # parser.add_argument('xmax', type=float)
    # parser.add_argument('vmin', type=float)
    # parser.add_argument('vmax', type=float)
    # parser.add_argument('simtype', type=int, choices=range(1,5), help='1: weak LD, 2: strong LD, 3: strong two-stream, 4: ''colder'' two-stream')
    # parser.add_argument('vth',type=float)
    # parser.add_argument('vstr', type=float)
    # parser.add_argument('normalized_wavenumber', metavar='kn', type=float)
    # parser.add_argument('amp', type=float)
    # parser.add_argument('initial_height', metavar='height', type=int)
    # parser.add_argument('greens_epsilon', metavar='eps', type=float)
    # parser.add_argument('num_steps', metavar='nt', type=int)
    # parser.add_argument('diag_freq', type=int)
    # parser.add_argument('dt', type=float)
    # parser.add_argument('--gpu',dest='use_gpu', action='store_true', help='boolean switch for using gpu')
    # parser.add_argument('--treecode','-tc',dest='use_treecode', action='store_true', help='boolean switch for using treecode')
    # parser.add_argument('--beta','-tc_b',type=float,default=-1.0, help='tunable treecode accuracy parameter. Barytree defaults to using this if 0<= beta <= 1.  Beta->0 is less accurate, Beta->1 is more accurate')
    # parser.add_argument('--mac','-tc_mac',type=float,default=-1.0,help='treecode multipole acceptance criterion. 0 <= mac <= 1')
    # parser.add_argument('--degree','-tc_d',type=int,default=-1,help='degree of treecode interpolation, must be nonnegative integer')
    # parser.add_argument('--maxSourceLeafSize','-s_l',type=int,default=200,help='maximum number of particles per source leaf')
    # parser.add_argument('--maxTargetLeafSize','-t_l',type=int,default=200,help='maximum number of particles per target leaf')
    # parser.add_argument('--root_dir',help='where to store simulation data')
    # args = parser.parse_args()

    # remesh_period = 1
    # FS = FarrsightInterface(args.project_name, args.xmin, args.xmax, args.vmin, args.vmax,\
    #                         SimType(args.simtype), 
    #                         args.normalized_wavenumber, args.amp, args.vth, args.vstr,\
    #                         args.initial_height,args.greens_epsilon,\
    #                         args.use_treecode, args.beta, args.mac,\
    #                         args.degree, args.maxSourceLeafSize, args.maxTargetLeafSize,\
    #                         args.num_steps, remesh_period, args.diag_freq, args.dt,\
    #                         root_dir = args.root_dir, can_do_movie=can_do_movie)

    # FS.run_farrsight(use_gpu = args.use_gpu)




# %%

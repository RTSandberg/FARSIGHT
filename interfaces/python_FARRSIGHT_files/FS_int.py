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

import numpy as np
import os # path.exists, makedirs, getcwd
import shutil # copy2

# import sys # path.append
import time # time

can_do_movie = True
try: 
    from matplotlib.animation import FFMpegWriter
    import matplotlib.animation as manimation
except:
    print('Unable to load ffmpeg.  Movie writer not accessible')
    can_do_movie = False


# import matplotlib
# from matplotlib.collections import PatchCollection
# from matplotlib.patches import Rectangle, Polygon
# plt.rcParams.update({'font.size': 18})
import sys
sys.path.append('/Users/ryansand/Documents/plasma_codes/FARRSIGHT/biquadratic_rk4_source/interfaces/python_FARRSIGHT_files')
# sys.path.append('/Users/ryansand/Documents/PlasmaProjects/Lagrangian_particle_method/heat_lamps_projects/amr/08_biquadratic_work/05_biquadratic_cpp/biquadratic_rk4_source')

#------ custom imports ---
# import python_FARRSIGHT_files
from . import FARRSIGHT_types as FST
# from python_FARRSIGHT_files import FARRSIGHT_types as FST
from python_FARRSIGHT_files import diagnostic_plot as dp
from python_FARRSIGHT_files import deck_utilities as deck
from python_FARRSIGHT_files import logf
from python_FARRSIGHT_files import make_dirs
from python_FARRSIGHT_files import panel_height
from python_FARRSIGHT_files import phase_space_vis as phase
from python_FARRSIGHT_files import plot_e
from python_FARRSIGHT_files import run
# --- end imports ---


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
        simulation_dictionary = deck.update_dictionary(deck_dir=deck_dir, deck_name = args.deck_name, **dict_args)
        sim_dir_str, directories_found = deck.generate_standard_names_dirs(simulation_dictionary, root_dir=args.root_dir)
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
        simulation_dictionary = deck.update_dictionary(deck_dir=deck_dir, deck_name = args.deck_name, **dict_args)
        
        dictionaries_found = True
        # end not-from-deck option
    sd = simulation_dictionary
    if sim_dir is not None:
        print('sim directory', sim_dir_str)

    if args.run:
        run.run_sim(sim_dir=sim_dir, deck_dir=deck_dir, deck_name=args.deck_name, use_gpu=args.use_gpu)

    if args.plot_diagnostics:
        iterations = np.arange(0,sd['num_steps']+1, sd['diag_period'])
        final_iter = iterations[-1]
        for species in sd['species_list']:
            sp_name = species['name']
            dp.sim_diagnostics_sample(simulation_dictionary, sp_name, sim_dir=sim_dir)
            if FST.SimType(sd['sim_type']) is not FST.SimType.FRIEDMAN_BEAM:
                # phase.plot_phase_space(sim_dir, simulation_dictionary, int(45.0 / simulation_dictionary['dt']/simulation_dictionary['diag_period']), FST.sim_type_to_flim[FST.SimType(simulation_dictionary['sim_type'])])
                phase.plot_phase_space(sim_dir, simulation_dictionary, sp_name, final_iter, FST.sim_type_to_flim[FST.SimType(simulation_dictionary['sim_type'])])
            if args.show_panels or args.panels_movie:
                panel_height.plot_height(sim_dir, sd, sp_name, final_iter)

    if args.plot_times is not None:
        t1 = time.time()
        iterations = np.arange(0,sd['num_steps']+1,sd['diag_period'])
        diag_ts = iterations * sd['dt']
        for species in sd['species_list']:
            for t0 in args.plot_times:
                print('generating time',t0,'phase space')
                t0_ind = np.nonzero(diag_ts <= t0)[0][-1]
                t0_iter = int(iterations[t0_ind])
                t0_sim = diag_ts[t0_ind]
                if FST.ICsType(species['ics_type']) is FST.ICs.COLDER_TWO_STREAM:
                    flim = (0,0.5/np.sqrt(2*np.pi)/species['pth'])
                else:
                    flim = FST.ics_type_to_flim[FST.ICsType(species['ics_type'])]
                phase.plot_phase_space(sim_dir,sd,sp_name, t0_iter,flim)
                if args.logf_movie or args.logf_plot:
                    logf.plot_logf(sim_dir, sd, sp_name, t0_iter)
                    
                if species['adaptively_refine'] or args.show_panels or args.panels_movie:
                    if 'p_height' in species:
                        hmax = species['max_height'] - species['p_height']
                    else:
                        hmax = species['max_height']
                    panel_height.plot_height(sim_dir, sd, sp_name, t0_iter, height_range=[species['initial_height'],hmax])
        t2 = time.time()
        print(f'spent {t2-t1:.1f}s making phase space (and log or height) images')

    if args.phase_movie:
        for species in sd['species_list']:
            sp_name = species['name']
            icstype = species['ics_type']
            if icstype == 1:
                flim = (0,0.44)
            elif icstype == 2:
                flim = (0,0.47)
            elif icstype == 3:
                flim = (0,0.45)
            elif icstype == 4:
                flim = (0,0.8)
            do_show_panels = False
            phase.phase_movie(sim_dir, simulation_dictionary, sp_name, do_show_panels, flim=flim, can_do_movie=can_do_movie)
            if args.show_panels:
                do_show_panels = True
                phase.phase_movie(sim_dir, simulation_dictionary, sp_name, do_show_panels, flim=flim, can_do_movie=can_do_movie)

    for species in sd['species_list']:
        sp_name = species['name']
        if args.panels_movie:
            if 'p_height' in species:
                hmax = species['max_height'] - species['p_height']
            else:
                hmax = species['max_height']
            panel_height.panel_height_movie(sim_dir, simulation_dictionary, sp_name,\
                height_range=[species['initial_height'],hmax],\
                can_do_movie=can_do_movie)

        if args.logf_movie:
            logf.logf_movie(sim_dir, simulation_dictionary, sp_name, can_do_movie=can_do_movie)
    

    print('done with this run')

# %%

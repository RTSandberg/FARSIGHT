
import os # makedirs, path.exists
import shutil # rmtree
import subprocess # run, PIPE
import time

import deck_utilities
import make_dirs

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


    sd = deck_utilities.deck_to_dict(deck_dir = deck_dir_str, deck_name = deck_name)
    
    output_dir = sim_dir_str + 'simulation_output/'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    for species in sd['species_list']:
        species_dir = output_dir + species['name'] + '/'
        if not os.path.exists(species_dir):
            os.makedirs(species_dir)
        if not os.path.exists(species_dir + 'panels/'):
            os.makedirs(species_dir + 'panels/')
        if not os.path.exists(species_dir + 'xs/'):
            os.makedirs(species_dir + 'xs/')
        if not os.path.exists(species_dir + 'vs/'):
            os.makedirs(species_dir + 'vs/')
        if not os.path.exists(species_dir + 'es/'):
            os.makedirs(species_dir + 'es/')
        if not os.path.exists(species_dir + 'fs/'):
            os.makedirs(species_dir + 'fs/')
        if not os.path.exists(species_dir + 'qws/'):
            os.makedirs(species_dir + 'qws/')

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

    sim_dir, directories_found = make_dirs.generate_standard_names_dirs(simulation_dictionary,root_dir)
    deck_utilities.dict_to_deck(simulation_dictionary, deck_dir=sim_dir)
    simulation_has_run = run_sim(sim_dir, use_gpu=use_gpu)
    return sim_dir, directories_found, simulation_has_run


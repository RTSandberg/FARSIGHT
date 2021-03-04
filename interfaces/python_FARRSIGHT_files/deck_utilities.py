import json # dump, load

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
        {'sim_name':'template',\
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


import os

import FARRSIGHT_types as FST

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
    bcs_string = FST.BoundaryConditions(0).name
    if 'bcs' in sd:
        bcs_string = FST.BoundaryConditions(sd['bcs']).name
    bcs_string += '_bcs'
        
    tc_string = ''
    # if 'use_treecode' not in sd
    if 'use_treecode' in sd and sd['use_treecode']:
        if 0 <= sd['beta'] <= 1:
            tc_string = 'tc_beta_%.2f'%sd['beta']# f'_beta_{sim_dict['beta']:.2f}'
        else:
            tc_string = 'tc_mac_%.2f_degree_%d_msource_%d_maxtarget_%d'%(sd['mac'],sd['degree'], sd['max_source'], sd['max_target'])
    else:
        tc_string = 'dsum'


    tf = sd['dt'] * sd['num_steps']
    physical_parameters =  bcs_string + '_nspecies_%i_tf_%.3f'%(len(sd['species_list']),tf)
    if 'remesh_period' in sd:
        remesh_period = sd['remesh_period']
    else:
        remesh_period = 1
    if 'quadrature' in sd:
        quad_str = FST.Quadrature(sd['quadrature']).name
    else:
        quad_str = 'TRAP'
    numerical_parameters = '%s_quadrature_vm_%.1f_g_eps_%.5f_dt_%.4f_remesh_period_%i_diag_freq_%i'%(quad_str,sd['vmax'], sd['greens_epsilon'], sd['dt'], remesh_period,sd['diag_period'])
#         sim_name = f'height0_{self.initial_height}_vm_{self.vmax:.1f}_g_eps_{self.greens_epsilon:.3f}_dt_{self.dt:.3f}_tf_{self.tf:.1f}_diag_freq_{self.diag_freq}' + tc_string
    sim_dir = simulations_dir + sd['sim_name'] + '/' + physical_parameters + '/'
    sim_dir += numerical_parameters + '/' + tc_string + '/'
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)
    # sim_dir, directories_found = make_dirs(sd['project_name'], sim_group, sim_name,root_dir=root_dir)
    return sim_dir, directories_found
# end generate_standard_names_dirs



import time
import numpy as np
from matplotlib import pyplot as plt

# import matplotlib
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Polygon
plt.rcParams.update({'font.size': 18})

import FARRSIGHT_types as FST

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
    plt.savefig(sim_dir_str + 'l2E.png')
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
    plt.savefig(sim_dir_str + 'maxf_conservation.png')
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
    plt.savefig(sim_dir_str + 'relative_maxf_conservation.png')
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
    plt.savefig(sim_dir_str + 'minf_conservation.png')
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
    plt.savefig(sim_dir_str + 'relative_minf_conservation.png')
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
    plt.savefig(sim_dir_str + 'frac_negatives.png')
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
    plt.savefig(sim_dir_str + 'frac_negative_mass.png')
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
    plt.savefig(sim_dir_str + 'num_points.png')
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
    plt.savefig(sim_dir_str + 'relative_charge_conservation.png')
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
    plt.savefig(sim_dir_str + 'momentum_conservation.png')
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
    plt.savefig(sim_dir_str + 'relative_l1f_conservation.png')
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
    plt.savefig(sim_dir_str + 'relative_l2f_conservation.png')
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
    plt.savefig(sim_dir_str + 'relative_energy_conservation.png')
    plt.close()
    # done plotting
    t2 = time.time()
    print(f'Done plotting diagnostics')
    print(f'Diagnostic collection and plot time {t2-t1:.3f}s')
    #--------------- f(v) plots --------
    nx0 = 2**(sd["initial_height"]+1) + 1
    dx0 = (sd["xmax"] - sd["xmin"]) / nx0
    for test_time in test_times:
        if simtype in [FST.SimType.WEAK_LD, FST.SimType.STRONG_LD]:# create new diagnostic : f(x=a,t=b,v)
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
                plt.savefig(sim_dir_str + f'fv_t_{simtime:.0f}.png')
                plt.close()
            except IndexError:
                print(f'time {test_time:.1f} is less than all simulation times')
        if simtype in [FST.SimType.STRONG_TWO_STREAM]:
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
                plt.savefig(sim_dir_str + f'fv_t_{simtime:.1f}.png')
                plt.close()
            except IndexError:
                print(f'time {test_time:.1f} is less than all simulation times')


        


# end sim_diagnostics_sample
    
def diagnostics_sample_standard_tree(simulation_dictionary,root_dir=None):

    sim_dir, directories_found = generate_standard_names_dirs(simulation_dictionary,root_dir)
    sim_diagnostics_sample(simulation_dictionary, sim_dir=sim_dir)

  
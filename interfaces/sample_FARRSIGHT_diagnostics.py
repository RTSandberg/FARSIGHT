

# if __name__ == '__main__':
import numpy as np
from matplotlib import pyplot as plt
# import os # path.exists, makedirs
# import subprocess # run, PIPE
import time # time
# import shutil # rmtree
import sys # sys.path
import matplotlib
can_do_movie = True
try: 
    from matplotlib.animation import FFMpegWriter
    import matplotlib.animation as manimation
except:
    print('Unable to load ffmpeg.  Movie writer not accessible')
    can_do_movie = False
# important: will need to import FARRSIGHT_Interface somehow!
# sys.path.append(<path to FARRSIGHT_Interface>)
import FARRSIGHT_Interface as FS_int
plt.rcParams.update({'font.size': 14})

# create FARRSIGHT interface

import argparse

parser = argparse.ArgumentParser(description='Python command-line interface for FARRSIGHT')
parser.add_argument('project_name')
parser.add_argument('xmin', type=float)
parser.add_argument('xmax', type=float)
parser.add_argument('vmin', type=float)
parser.add_argument('vmax', type=float)
parser.add_argument('simtype', type=int, choices=range(1,5), help='1: weak LD, 2: strong LD, 3: strong two-stream, 4: ''colder'' two-stream')
parser.add_argument('vth',type=float)
parser.add_argument('vstr', type=float)
parser.add_argument('normalized_wavenumber', metavar='kn', type=float)
parser.add_argument('amp', type=float)
parser.add_argument('initial_height', metavar='height', type=int)
parser.add_argument('greens_epsilon', metavar='eps', type=float)
parser.add_argument('num_steps', metavar='nt', type=int)
parser.add_argument('diag_freq', type=int)
parser.add_argument('dt', type=float)
# parser.add_argument('--gpu',dest='use_gpu', action='store_true', help='boolean switch for using gpu')
parser.add_argument('--treecode','-tc',dest='use_treecode', action='store_true', help='boolean switch for using treecode')
parser.add_argument('--beta','-tc_b',type=float,default=-1.0, help='tunable treecode accuracy parameter. Barytree defaults to using this if 0<= beta <= 1.  Beta->0 is less accurate, Beta->1 is more accurate')
parser.add_argument('--mac','-tc_mac',type=float,default=-1.0,help='treecode multipole acceptance criterion. 0 <= mac <= 1')
parser.add_argument('--degree','-tc_d',type=int,default=-1,help='degree of treecode interpolation, must be nonnegative integer')
parser.add_argument('--maxSourceLeafSize','-s_l',type=int,default=200,help='maximum number of particles per source leaf')
parser.add_argument('--maxTargetLeafSize','-t_l',type=int,default=200,help='maximum number of particles per target leaf')
parser.add_argument('--root_dir',help='where to store simulation data')
parser.add_argument('--phase_movie','-pm',dest='do_phase_movie',action='store_true',help='make phase space movie')
parser.add_argument('--logf_movie','-logm',dest='do_logf_movie',action='store_true',help='make log  movie')
args = parser.parse_args()

remesh_period = 1
FS = FS_int.FarrsightInterface(args.project_name, args.xmin, args.xmax, args.vmin, args.vmax,\
                        FS_int.SimType(args.simtype), args.normalized_wavenumber, args.amp, 
                        args.vth, args.vstr,\
                        args.initial_height, args.greens_epsilon, args.use_treecode, args.beta, args.mac,\
                        args.degree, args.maxSourceLeafSize, args.maxTargetLeafSize,\
                        args.num_steps, remesh_period,\
                        args.diag_freq, args.dt, root_dir = args.root_dir, can_do_movie=can_do_movie)


if args.do_phase_movie:
    FS.make_movie(do_show_panels=False)
if args.do_logf_movie:
    FS.logf_movie()

# get diagnostics
print('collecting diagnostics')
t1 = time.time()
output_dir = FS.sim_dir + 'simulation_output/'
diag_times = np.arange(0,FS.num_steps+1,FS.diag_freq) * FS.dt
total_charge = np.zeros_like(diag_times)
total_momentum = np.zeros_like(diag_times)
total_kinetic = np.zeros_like(diag_times)
total_potential = np.zeros_like(diag_times)
total_entropy = np.zeros_like(diag_times)
total_negative_area = np.zeros_like(diag_times)
num_negative = np.zeros_like(diag_times)

l1f = np.zeros_like(diag_times)
l2f = np.zeros_like(diag_times)

Lv = FS.vmin - FS.vmax

for ii, t in enumerate(diag_times):
    iter_num = ii * FS.diag_freq
    simtime = iter_num * FS.dt

    xs = np.fromfile(output_dir + f'xs/xs_{iter_num}')
    vs = np.fromfile(output_dir + f'vs/vs_{iter_num}')
    fs = np.fromfile(output_dir + f'fs/fs_{iter_num}')
    qws = np.fromfile(output_dir + f'qws/qws_{iter_num}')
    es = np.fromfile(output_dir + f'es/es_{iter_num}')
    panel_point_inds = np.fromfile(output_dir + f'panels/leaf_point_inds_{iter_num}',dtype='int32')
    pinds = np.reshape(panel_point_inds, (int(panel_point_inds.size/9),9))
    
    total_charge[ii] = np.sum(qws)
    total_momentum[ii] = 1/FS.qm * np.dot(vs, qws)
    total_kinetic[ii] = 1/FS.qm * np.dot(vs **2, qws)
    
    # potential
    x_sort, sort_inds = np.unique(xs, return_index=True)
    E_sort = es[sort_inds]
    mod_x_inds = np.nonzero((x_sort >= 0) * (x_sort < FS.Lx))
    modx = x_sort[mod_x_inds]
    modE = E_sort[mod_x_inds]
    dxs = .5 * np.hstack([modx[1]+FS.Lx-modx[-1], modx[2:] - modx[:-2], modx[0]+FS.Lx - modx[-2]])
    if np.nonzero(dxs<=0)[0].size > 0:
        raise ValueError
    total_potential[ii] = np.dot(dxs, modE**2)
    # entropy
    where_negative = np.nonzero(fs < 0)
    total_negative_area[ii] = 1/FS.q * np.sum(qws[where_negative])
    num_negative[ii] = where_negative[0].size
    where_positive = np.nonzero(fs > 0)
    total_entropy[ii] = 1/FS.q * np.dot(fs[where_positive] * np.log10(fs[where_positive]), qws[where_positive])
    
    l1f[ii] = 1/FS.q * np.sum(abs(qws))
    l2f[ii] = 1/FS.q * np.dot(qws,fs)

# done collecting diagnostics
t2 = time.time()
print('Done collecting diagnostics')
print(f'Took {t2-t1:.3f}s to collect')

print('plotting diagnostics')

plt.figure()
plt.title('l2 E')
l2e = np.sqrt(total_potential)
plt.semilogy(diag_times, l2e)

#weak LD:
# if type_str == 'weak_LD':
if FS.simtype is FS_int.SimType.WEAK_LD:
    g_t = .1533
    peak_ind = 0
    plt.semilogy(diag_times, l2e[peak_ind]*np.exp(-g_t*(diag_times - diag_times[peak_ind])), '-.')

#strong LD:
# if type_str == 'strong_LD':
# if simtype == 2:
if FS.simtype is FS_int.SimType.STRONG_LD:
    g_strong1 = .2920
    g_strong2 = .0815
    damping1_ind = 5
    damping1_times = diag_times[:25]
    grow1_start_ind = 28
    grow1_times = diag_times[grow1_start_ind:90]
    grow1_ind = 19
    grow1_ind_e = grow1_start_ind + grow1_ind
    plt.semilogy(damping1_times,l2e[damping1_ind]*np.exp(- g_strong1 * (damping1_times - damping1_times[damping1_ind])),'-.',label=r'$\sim e^{- %.05f t}$'%g_strong1)
    plt.semilogy(grow1_times,l2e[grow1_ind_e]*np.exp(g_strong2*(grow1_times-grow1_times[grow1_ind])),'-.k')

plt.grid()
#weak LD
# if type_str == 'weak_LD':
# if simtype == 1:
if FS.simtype is FS_int.SimType.WEAK_LD:
#     pass
    plt.ylim([1e-8,1e-1])
#strong LD
# if type_str == 'strong_LD' or type_str == 'strong_2str':
if FS.simtype in [2,3,4]:
    plt.ylim([1e-3, 1e1])
plt.savefig(FS.sim_dir + 'l2E.png') 
plt.close()

plt.figure()
plt.title(r'number negative f values')
plt.xlabel('t')
plt.plot(diag_times, num_negative)
plt.grid()
plt.savefig(FS.sim_dir + 'num_negatives.png')
plt.close()

plt.figure()
plt.title(r'percent variation in total charge, $\sum_i q w^n_i - \sum_i q w^0_i$')
plt.xlabel('t')
plt.plot(diag_times, 100 * (total_charge - total_charge[0])/total_charge[0])
plt.grid()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig(FS.sim_dir + 'percent_charge_conservation.png')
plt.close()

plt.figure()
plt.title(r'variation in total momentum, $\sum_i m v^n_i w^n_i - \sum_i m v^0_i w^0_i$')
plt.xlabel('t')
plt.plot(diag_times, total_momentum - total_momentum[0])
plt.grid()
plt.savefig(FS.sim_dir + 'momentum_conservation.png')
plt.close()

plt.figure()
plt.title(r'percent variation in $L_2(f^n)$, $\sum_i (f^n_i)^2 w^n_i - \sum_i (f^n_i)^2 w^0_i$')
plt.xlabel('t')
plt.plot(diag_times, 100*(l2f - l2f[0])/l2f[0])
plt.grid()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig(FS.sim_dir + 'percent_l2f_conservation.png')
plt.close()

plt.figure()
plt.title('percent energy variation')
plt.xlabel('t')
total_energy = total_potential + total_kinetic
plt.plot(diag_times, 100 * (total_energy - total_energy[0])/total_energy[0])
plt.grid()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig(FS.sim_dir + 'relative_energy_conservation.png')
plt.close()
print('done plotting diagnostics')

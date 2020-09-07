import numpy as np
from matplotlib import pyplot as plt
import os # path.exists, makedirs
import subprocess # run, PIPE
import time # time
import shutil # rmtree
import sys

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

def make_dirs(project_name, sim_group, sim_name, 
                          root_dir = None):
    """
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

    output_dir = sim_dir + 'simulation_output/'
    if not os.path.exists(output_dir):
        directories_found = False
        os.makedirs(output_dir)
    if not os.path.exists(output_dir + 'panels/'):
        directories_found = False
        os.makedirs(output_dir + 'panels/')
    if not os.path.exists(output_dir + 'xs/'):
        directories_found = False
        os.makedirs(output_dir + 'xs/')
    if not os.path.exists(output_dir + 'vs/'):
        directories_found = False
        os.makedirs(output_dir + 'vs/')
    if not os.path.exists(output_dir + 'es/'):
        directories_found = False
        os.makedirs(output_dir + 'es/')
    if not os.path.exists(output_dir + 'fs/'):
        directories_found = False
        os.makedirs(output_dir + 'fs/')
    if not os.path.exists(output_dir + 'qws/'):
        directories_found = False
        os.makedirs(output_dir + 'qws/')
    return sim_dir, directories_found

class FarrsightInterface:
    """object for interacting with FARRSIGHT
    
    Members
    ---
    
    """
    def __init__(self, project_name,\
                xmin, xmax,\
                vmin, vmax, simtype,\
                normalized_wavenumber, amp, vth, vstr,\
                initial_height, greens_epsilon,\
                use_treecode, beta,\
                mac, degree, max_source, max_target,\
                num_steps, remesh_period, diag_freq, dt,\
                root_dir = None, can_do_movie = False):
        self.project_name = project_name
        self.xmin = xmin
        self.xmax = xmax
        self.Lx = xmax - xmin
        self.vmin = vmin
        self.vmax = vmax
        self.Lv = vmax - vmin
        self.simtype = simtype
        self.q = -1.0
        self.m = 1.0
        self.qm = self.q / self.m
        self.normalized_wavenumber = normalized_wavenumber
        self.amp = amp
        self.vth = vth
        self.vstr = vstr
        self.initial_height = initial_height
        self.greens_epsilon = greens_epsilon
        if use_treecode:
            self.use_treecode = 1
        else:
            self.use_treecode = 0
        self.beta = beta
        self.mac = mac
        self.degree = degree
        self.max_source = max_source
        self.max_target = max_target
        self.num_steps = num_steps
        self.remesh_period = remesh_period
        self.diag_freq = diag_freq
        self.dt = dt
        
        self.tf = self.num_steps * self.dt
        
        self.can_do_movie = can_do_movie

        self.simulation_has_run = False
        
        self.generate_standard_names_dirs()
    
    def __repr__(self):
        return f'{self.__class__.__name__}(...)'
    
    def __str__(self):
        return f'A {self.__class__.__name__} for a '+\
            f'level {self.initial_height} sim'
    
    def run_farrsight(self, use_gpu = False):
        # clear output dirs if they exist
        # make dirs if not made
        
        output_dir = self.sim_dir + 'simulation_output/'
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
        farrsight_args = [farrsight_exe, self.sim_dir,\
                          self.xmin, self.xmax,\
                          self.vmin, self.vmax,\
                          self.simtype,\
                          self.normalized_wavenumber,\
                          self.amp, self.vth, self.vstr,\
                          self.initial_height, \
                          self.greens_epsilon,\
                          self.use_treecode, self.beta,\
                          self.mac, self.degree, \
                          self.max_source, self.max_target,\
                          self.num_steps, self.remesh_period,\
                          self.diag_freq, self.dt]
        farrsight_args = [str(arg) for arg in farrsight_args]
        
        print('running sim...')
        t1 = time.time()
        with open(self.sim_dir + 'sim_out.txt','w') as log:
            with open(self.sim_dir + 'sim_err.txt','w') as err_log:
                proc = subprocess.run(farrsight_args,stdout=log,stderr=err_log)
        t2 = time.time()

        if proc.returncode == 0:
            print('done!')
            self.simulation_has_run = True
        else:
            self.simulation_has_run = False
            print('simulation stopped with error code %i'%proc.returncode)

        with open(self.sim_dir + 'sim.txt','a') as log:
            log.write(f'compute time {t2-t1:.03f}s')
        print(f'simulation compute time {t2-t1:.03f}s')
    # end run_farrsight
    
    def generate_standard_names_dirs(self):
        tc_string = ''
        if self.use_treecode:
            if 0 <= self.beta <= 1:
                tc_string = f'_beta_{self.beta:.2f}'
            else:
                tc_string = f'_mac_{self.mac:.2f}_degree_{self.degree}_msource_{self.max_source}_maxtarget_{self.max_target}'
        else:
            tc_string = '_dsum'

        sim_group = f'vth_{self.vth:.3f}_vstr_{self.vstr:.3f}_amp_{self.amp:.3f}_normal_k_{self.normalized_wavenumber:.3f}'
        sim_name = f'height0_{self.initial_height}_vm_{self.vmax:.1f}_g_eps_{self.greens_epsilon:.3f}_dt_{self.dt:.3f}_tf_{self.tf:.1f}_diag_freq_{self.diag_freq}' + tc_string

        self.sim_dir, directories_found = make_dirs(self.project_name, sim_group, sim_name)
        self.simulation_has_run = directories_found
    # end generate_standard_names_dirs

    def plot_things(self,step_ii,flim, do_save = False):
        
        if not self.simulation_has_run:
            print('unable to plot; simulation has not run or had errors')
            return
        
        simtime = step_ii * self.dt
        
        output_dir = self.sim_dir + 'simulation_output/'
        xs = np.fromfile(output_dir + f'xs/xs_{step_ii}')
        vs = np.fromfile(output_dir + f'vs/vs_{step_ii}')
        fs = np.fromfile(output_dir + f'fs/fs_{step_ii}')
        es = np.fromfile(output_dir + f'es/es_{step_ii}')
        panels = np.fromfile(output_dir + f'panels/leaf_point_inds_{step_ii}',dtype='int32')
        num_panels = int(panels.size/9)
        panels = np.reshape(panels, (num_panels,9))
        panels_fs = np.zeros(4*num_panels)
        
        
        fig,ax = plt.subplots(2,1)
        ax0 = ax[0]
        ax1 = ax[1]

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

        ax0.set_xlim(self.xmin,self.xmax)
        ax0.set_ylim(self.vmin,self.vmax)
        ax0.set_ylabel('v')
        ax0.set_title(f'phase space with {num_panels} panels, t={simtime:.03f}')
        
    #     plt.figure()
        ax1.set_title(f'electric field with {num_panels} panels, t={simtime:.03f}')
        ax1.plot(xs,es,'.')
        ax1.set_xlabel('x')
        plt.grid()
        plt.tight_layout()
        if do_save:
            plt.savefig(self.sim_dir + f'phase_space_time_{step_ii}.png')
    plt.close()
    # end plot_things

    # %%time
# do_show_panels = False
def make_movie(self,do_show_panels,flim=(0,.3)):

    print('starting phase space movie')
    t1 =time.time()
    if not self.simulation_has_run:
        print('unable to plot; simulation has not run or had errors')
        return
    if not self.can_do_movie:
        print('unable to load ffmpeg writer, movie capability inaccessible')
        return
    panel_string = ''

    if do_show_panels:
        panel_string = 'show_panels_'

    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='cpp test', artist='Matplotlib',
                    comment='')
    writer = FFMpegWriter(fps=5, metadata=metadata)

    fig, ax = plt.subplots(figsize=(12,10))

    with writer.saving(fig, self.sim_dir + panel_string + 'phase_space'+ ".mp4", dpi=100):


#     xs = np.fromfile(output_dir + f'xs/xs_{step_ii}')
#     vs = np.fromfile(output_dir + f'vs/vs_{step_ii}')
#     fs = np.fromfile(output_dir + f'fs/fs_{step_ii}')
#     es = np.fromfile(output_dir + f'es/es_{step_ii}')
#     panels = np.fromfile(output_dir + f'panels/leaf_point_inds_{step_ii}',dtype='int32')
#     num_panels = int(panels.size/9)
#     panels = np.reshape(panels, (num_panels,9))
#     panels_fs = np.zeros(4*num_panels)
#         flim = (0,.3)

        panels = np.fromfile(self.sim_dir + 'simulation_output/panels/leaf_point_inds_0',dtype='int32')
        num_panels = int(panels.size/9)
        panels = np.reshape(panels, (num_panels,9))
        panels_fs = np.zeros(4*num_panels)
        xs = np.fromfile(self.sim_dir + 'simulation_output/xs/xs_0')
        vs = np.fromfile(self.sim_dir + 'simulation_output/vs/vs_0')
        fs = np.fromfile(self.sim_dir + 'simulation_output/fs/fs_0')

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

        ax.set_xlim([self.xmin, self.xmax])
        ax.set_ylim([self.vmin,self.vmax])
        ax.set_xlabel('x')
        ax.set_ylabel('v')
        ax.set_title(f'phase space with {num_panels} panels, t=0.000')

        fig.canvas.draw()
        writer.grab_frame()


        print_update_frequency = int(np.ceil((self.num_steps+1)/6/self.diag_freq))
        print_update_counter = 0

        for iter_num in range(self.diag_freq, self.num_steps + 1, self.diag_freq):
            cb.remove()
            ax.collections.pop()

            panels = np.fromfile(self.sim_dir + f'simulation_output/panels/leaf_point_inds_{iter_num}',dtype='int32')
            num_panels = int(panels.size/9)
            panels = np.reshape(panels, (num_panels,9))
            panels_fs = np.zeros(num_panels*4)
            xs = np.fromfile(self.sim_dir + f'simulation_output/xs/xs_{iter_num}')
            vs = np.fromfile(self.sim_dir + f'simulation_output/vs/vs_{iter_num}')
            fs = np.fromfile(self.sim_dir + f'simulation_output/fs/fs_{iter_num}')

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
            ax.set_title(f'phase space with {num_panels} panels, t={iter_num*self.dt:.3f}')


            fig.canvas.draw()
            writer.grab_frame()


            if print_update_counter == print_update_frequency:
                print(f'Movie is about {iter_num/(self.num_steps+1)*100 :0.0f}% complete')
                print_update_counter = 0
#                 plt.savefig(mya.sim_dir + f'phase_space_image_t_{iter_num*dt}.svg')
                plt.savefig(self.sim_dir + f'phase_space_image_t_{iter_num*self.dt}.png')
            print_update_counter += 1
    print('phase space movie done!')
    t2 = time.time()
    print(f'phase space movie took {t2-t1:.3f}s')
    plt.close()
FarrsightInterface.make_movie = make_movie

# logf movie
# %%time
def logf_movie(self,flim=(-8,0)):
    print('starting logf movie')
    t1 = time.time()
    if not self.simulation_has_run:
        print('unable to plot; simulation has not run or had errors')
        return
    if not self.can_do_movie:
        print('unable to load ffmpeg writer, movie capability inaccessible')
        return
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='logf movie', artist='Matplotlib',
                    comment='')
    writer = FFMpegWriter(fps=5, metadata=metadata)

    fig, ax = plt.subplots()

    with writer.saving(fig, self.sim_dir+ 'logf_movie.mp4', dpi=100):



    #     flim = (-8,0)

        panels = np.fromfile(self.sim_dir + 'simulation_output/panels/leaf_point_inds_0',dtype='int32')
        num_panels = int(panels.size/9)
        panels = np.reshape(panels, (num_panels,9))
        panels_fs = np.zeros(4*num_panels)
        xs = np.fromfile(self.sim_dir + 'simulation_output/xs/xs_0')
        vs = np.fromfile(self.sim_dir + 'simulation_output/vs/vs_0')
        fs = np.fromfile(self.sim_dir + 'simulation_output/fs/fs_0')

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

        ax.set_xlim(self.xmin, self.xmax)
        ax.set_ylim(self.vmin, self.vmax)
        ax.set_xlabel('x')
        ax.set_ylabel('v')
        ax.set_title(f'log(f), t=0.000')

        fig.canvas.draw()
        writer.grab_frame()


        print_update_frequency = int(np.ceil((self.num_steps+1)/6/self.diag_freq))
        print_update_counter = 0

        for iter_num in range(self.diag_freq, self.num_steps + 1, self.diag_freq):
            cb.remove()
            ax.collections.pop()

            panels = np.fromfile(self.sim_dir + f'simulation_output/panels/leaf_point_inds_{iter_num}',dtype='int32')
            num_panels = int(panels.size/9)
            panels = np.reshape(panels, (num_panels,9))
            panels_fs = np.zeros(4*num_panels)
            xs = np.fromfile(self.sim_dir + f'simulation_output/xs/xs_{iter_num}')
            vs = np.fromfile(self.sim_dir + f'simulation_output/vs/vs_{iter_num}')
            fs = np.fromfile(self.sim_dir + f'simulation_output/fs/fs_{iter_num}')

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
            ax.set_title(f'log(f), t={iter_num*self.dt:.3f}')


            fig.canvas.draw()
            writer.grab_frame()


            if print_update_counter == print_update_frequency:
                print(f'Movie is about {iter_num/(self.num_steps+1)*100 :0.0f}% complete')
                print_update_counter = 0
    #             plt.savefig(mya.sim_dir + f'logf_phase_space_image_t_{iter_num*dt}.svg')
                plt.savefig(self.sim_dir + f'logf_phase_space_image_t_{iter_num*self.dt}.png')
            print_update_counter += 1

        plt.savefig(self.sim_dir + f'logf_phase_space_final_image_t_{iter_num*self.dt}.png')
    print('done with logf movie!')
    t2 = time.time()
    print(f'logf movie took {t2 - t1:.3f}s')
    plt.close()
FarrsightInterface.logf_movie = logf_movie


# main
if __name__ == '__main__':

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
    parser.add_argument('--gpu',dest='use_gpu', action='store_true', help='boolean switch for using gpu')
    parser.add_argument('--treecode','-tc',dest='use_treecode', action='store_true', help='boolean switch for using treecode')
    parser.add_argument('--beta','-tc_b',type=float,default=-1.0, help='tunable treecode accuracy parameter. Barytree defaults to using this if 0<= beta <= 1.  Beta->0 is less accurate, Beta->1 is more accurate')
    parser.add_argument('--mac','-tc_mac',type=float,default=-1.0,help='treecode multipole acceptance criterion. 0 <= mac <= 1')
    parser.add_argument('--degree','-tc_d',type=int,default=-1,help='degree of treecode interpolation, must be nonnegative integer')
    parser.add_argument('--maxSourceLeafSize','-s_l',type=int,default=200,help='maximum number of particles per source leaf')
    parser.add_argument('--maxTargetLeafSize','-t_l',type=int,default=200,help='maximum number of particles per target leaf')
    parser.add_argument('--root_dir',help='where to store simulation data')
    args = parser.parse_args()

    remesh_period = 1
    FS = FarrsightInterface(args.project_name, args.xmin, args.xmax, args.vmin, args.vmax,\
                            SimType(args.simtype), 
                            args.normalized_wavenumber, args.amp, args.vth, args.vstr,\
                            args.initial_height,args.greens_epsilon,\
                            args.use_treecode, args.beta, args.mac,\
                            args.degree, args.maxSourceLeafSize, args.maxTargetLeafSize,\
                            args.num_steps, remesh_period, args.diag_freq, args.dt,\
                            root_dir = args.root_dir, can_do_movie=can_do_movie)

    FS.run_farrsight(use_gpu = args.use_gpu)





import time
import numpy as np
from matplotlib import pyplot as plt
# import matplotlib
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Polygon

can_do_movie = True
try: 
    from matplotlib.animation import FFMpegWriter
    import matplotlib.animation as manimation
except:
    print('Unable to load ffmpeg.  Movie writer not accessible')
    can_do_movie = False

import make_dirs

plt.rcParams.update({'font.size': 18})

def plot_logf(sim_dir, simulation_dictionary, species, step_ii,flim=(-20,0), simulation_has_run = True, do_save = False):
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
    ps = np.fromfile(output_dir + f'ps/ps_{step_ii}')
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
        panel_ps = ps[panel]
        panel_fs = fs[panel]
#             weights = np.array([1,4,1,4,16,4,1,4,1])
#             panels_fs[ii] = 1./36. * np.dot(panel_fs,weights)

        p0 = [0,1,4,3]
        panels_fs[4*ii] = .25*sum(panel_fs[p0])
        rect_pts = np.vstack([panel_xs[p0],panel_ps[p0]]).T
        patches.append(Polygon(rect_pts))

        p1 = [1,2,5,4]
        panels_fs[4*ii+1] = .25*sum(panel_fs[p1])
        rect_pts = np.vstack([panel_xs[p1],panel_ps[p1]]).T
        patches.append(Polygon(rect_pts))

        p2 = [3,4,7,6]
        panels_fs[4*ii+2] = .25*sum(panel_fs[p2])
        rect_pts = np.vstack([panel_xs[p2],panel_ps[p2]]).T
        patches.append(Polygon(rect_pts))

        p3 = [4,5,8,7]
        panels_fs[4*ii+3] = .25*sum(panel_fs[p3])
        rect_pts = np.vstack([panel_xs[p3],panel_ps[p3]]).T
        patches.append(Polygon(rect_pts))

    p = PatchCollection(patches, cmap=plt.cm.jet)
    p.set_array(np.log10(panels_fs))
    p.set_clim(flim)
    ax0.add_collection(p)
    cb = fig.colorbar(p, ax=ax0)
    # cb.set_label('f')

    ax0.set_xlim(sd['xmin'], sd['xmax'])
    ax0.set_ylim(sd['pmin'],sd['pmax'])
    ax0.set_xlabel('x')
    ax0.set_ylabel('p')
    # ax0.set_title(f't={simtime:.03f}')
    plt.tight_layout()
    
    plt.savefig(sim_dir_str + f'logf_{simtime:.2f}.png')
    plt.close()
# end plot_logf


# logf movie
# %%time
def logf_movie(sim_dir, simulation_dictionary, species, simulation_has_run = True, can_do_movie = True, flim=(-12,0)):
    sd = simulation_dictionary
    output_dir = 'simulation_output/' + species + '/'
    sim_dir_str = ''
    if sim_dir is not None:
        sim_dir_str = sim_dir
        if sim_dir[-1] != '/':
            sim_dir_str += '/'
        output_dir = sim_dir_str + output_dir

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

    with writer.saving(fig, sim_dir_str+ 'logf_movie.mp4', dpi=100):



    #     flim = (-8,0)

        panels = np.fromfile(output_dir + 'panels/leaf_point_inds_0',dtype='int32')
        num_panels = int(panels.size/9)
        panels = np.reshape(panels, (num_panels,9))
        panels_fs = np.zeros(4*num_panels)
        xs = np.fromfile(output_dir + 'xs/xs_0')
        ps = np.fromfile(output_dir + 'ps/ps_0')
        fs = np.fromfile(output_dir + 'fs/fs_0')

        patches = []
        for ii, panel in enumerate(panels):
            panel_xs = xs[panel]
            panel_ps = ps[panel]
            panel_fs = fs[panel]

            p0 = [0,1,4,3]
            panels_fs[4*ii] = .25*sum(panel_fs[p0])
            rect_pts = np.vstack([panel_xs[p0],panel_ps[p0]]).T
            patches.append(Polygon(rect_pts))

            p1 = [1,2,5,4]
            panels_fs[4*ii+1] = .25*sum(panel_fs[p1])
            rect_pts = np.vstack([panel_xs[p1],panel_ps[p1]]).T
            patches.append(Polygon(rect_pts))

            p2 = [3,4,7,6]
            panels_fs[4*ii+2] = .25*sum(panel_fs[p2])
            rect_pts = np.vstack([panel_xs[p2],panel_ps[p2]]).T
            patches.append(Polygon(rect_pts))

            p3 = [4,5,8,7]
            panels_fs[4*ii+3] = .25*sum(panel_fs[p3])
            rect_pts = np.vstack([panel_xs[p3],panel_ps[p3]]).T
            patches.append(Polygon(rect_pts))

        p = PatchCollection(patches, cmap=plt.cm.jet)

        p.set_array(np.log10(panels_fs))
        p.set_clim(flim)
        ax.add_collection(p)
        cb = fig.colorbar(p, ax=ax)
        cb.set_label('f')

        ax.set_xlim(sd["xmin"], sd["xmax"])
        ax.set_ylim(sd["pmin"], sd["pmax"])
        ax.set_xlabel('x')
        ax.set_ylabel('p')
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
            ps = np.fromfile(output_dir + f'ps/ps_{iter_num}')
            fs = np.fromfile(output_dir + f'fs/fs_{iter_num}')

            patches = []
            for ii, panel in enumerate(panels):
                panel_xs = xs[panel]
                panel_ps = ps[panel]
                panel_fs = fs[panel]

                p0 = [0,1,4,3]
                panels_fs[4*ii] = .25*sum(panel_fs[p0])
                rect_pts = np.vstack([panel_xs[p0],panel_ps[p0]]).T
                patches.append(Polygon(rect_pts))

                p1 = [1,2,5,4]
                panels_fs[4*ii+1] = .25*sum(panel_fs[p1])
                rect_pts = np.vstack([panel_xs[p1],panel_ps[p1]]).T
                patches.append(Polygon(rect_pts))

                p2 = [3,4,7,6]
                panels_fs[4*ii+2] = .25*sum(panel_fs[p2])
                rect_pts = np.vstack([panel_xs[p2],panel_ps[p2]]).T
                patches.append(Polygon(rect_pts))

                p3 = [4,5,8,7]
                panels_fs[4*ii+3] = .25*sum(panel_fs[p3])
                rect_pts = np.vstack([panel_xs[p3],panel_ps[p3]]).T
                patches.append(Polygon(rect_pts))

    #             panel_ps = ps[panel]
    #             panel_fs = fs[panel]
    #             panels_fs[ii] = np.log10(.25*sum(panel_fs))
    #             rect_pts = np.vstack([panel_xs,panel_ps]).T
    #             patches.append(Polygon(rect_pts))

            p = PatchCollection(patches, cmap=plt.cm.jet)

            p.set_array(np.log10(panels_fs))
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
    #             plt.savefig(mya.sim_dir + f'logf_phase_space_image_t_{iter_num*dt}.svg')
                plt.savefig(sim_dir_str + f'logf_phase_space_image_t_{iter_num*sd["dt"]}.png')
            print_update_counter += 1

        plt.savefig(sim_dir_str + f'logf_phase_space_final_image_t_{iter_num*sd["dt"]}.png')
    print('done with logf movie!')
    t2 = time.time()
    print(f'logf movie took {t2 - t1:.3f}s')
    plt.close()
# end logf movie

def logf_movie_standard_tree(simulation_dictionary,species, root_dir=None,simulation_has_run = True, can_do_movie = True, flim=(-8,0)):

    sim_dir, directories_found = make_dirs.generate_standard_names_dirs(simulation_dictionary,root_dir)
    logf_movie(sim_dir, simulation_dictionary, species)
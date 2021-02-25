
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

plt.rcParams.update({'font.size': 18})

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
    # mymap = plt.cm.get_cmap('gist_rainbow_r',ncolors)
    mymap = plt.cm.get_cmap('jet',ncolors)
    

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
    plt.savefig(sim_dir_str + f'panel_height_image_t_{iter_num*sd["dt"]}.png')
    print('panel height plot done!')
    t2 = time.time()
    print(f'panels plot took {t2-t1:.3f}s')
    plt.close()
# end plot height


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

    fig, ax = plt.subplots(figsize=(12,10))

    with writer.saving(fig, sim_dir_str + 'panel_heights'+ ".mp4", dpi=100):



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
        mymap = plt.cm.get_cmap('jet',ncolors)
        p = PatchCollection(patches, mymap)
        p.set_array(panels_fs)
        p.set_clim(height_range[0] - 0.5, height_range[1] + 0.5)
        ax.add_collection(p)
        # cb = fig.colorbar(p, ax=ax)
        plt_ticks = np.arange(height_range[0],height_range[1]+1)
        cb = fig.colorbar(p, ax=ax, ticks=plt_ticks)
        cb.set_label('panel height')

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
            cb = fig.colorbar(p, ax=ax)
            cb = fig.colorbar(p, ax=ax, ticks=plt_ticks)
            cb.set_label('panel height')
            ax.set_title(f't={iter_num*sd["dt"]:.3f}')


            fig.canvas.draw()
            writer.grab_frame()


            if print_update_counter == print_update_frequency:
                print(f'Movie is about {iter_num/(sd["num_steps"]+1)*100 :0.0f}% complete')
                print_update_counter = 0
#                 plt.savefig(mya.sim_dir + f'phase_space_image_t_{iter_num*dt}.svg')
                plt.savefig(sim_dir_str + f'panel_height_image_t_{iter_num*sd["dt"]}.png')
            print_update_counter += 1

    plt.savefig(sim_dir_str + f'panel_height_image_t_{iter_num*sd["dt"]}.png')
    print('panel height movie done!')
    t2 = time.time()
    print(f'panels movie took {t2-t1:.3f}s')
    plt.close()
# end panels movie

def panel_height_movie_standard_tree(simulation_dictionary,root_dir=None,height_range=(7,11), simulation_has_run=True, can_do_movie=True):

    sim_dir, directories_found = generate_standard_names_dirs(simulation_dictionary,root_dir)
    panel_height_movie(sim_dir, simulation_dictionary, height_range=height_range, simulation_has_run=simulation_has_run, can_do_movie=can_do_movie)


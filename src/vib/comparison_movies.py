"""
Compare four simulations with different time step in the
same movie.
Problem: ffmpeg/avconv will only make a movie of the first
8 periods. scitools movie can make html player that plays
all 30 periods.
"""

from vib_undamped import solver, visualize_front
from math import pi
import os, shutil, glob

def run_simulations(N, dt0, num_periods):
    """
    Run N simulations where the time step is halved in each
    simulation, starting with dt0.
    Make subdirectories tmp_case0, tmp_case1, etc with plot files
    for each simulation (tmp_*.png).
    """
    for i in range(N):
        dt = dt0/2.0**i
        u, t = solver(I=1, w=2*pi, dt=dt, T=num_periods)
        # visualize_front removes all old plot files :)
        visualize_front(u, t, I=1, w=2*pi, savefig=True,
                        skip_frames=2**i)
        # skip_frames is essential: for N=4 we have to store
        # only each 2**4=16-th file to get as many files
        # as for the dt0 simulation!

        # Move all plot files tmp_*.png for movie to a
        # separate directory. Delete that directory if it
        # exists and recreate it.
        dirname = 'tmp_case%d' % i
        if os.path.isdir(dirname):
            shutil.rmtree(dirname)  # remove directory (tree)
        os.mkdir(dirname)           # make new directory
        for filename in glob.glob('tmp_*.png'):
            # Move file to subdirectory dirname
            os.rename(filename, os.path.join(dirname, filename))

def make_movie(N):
    """
    Combine plot files in subdirectories tmp_case0,
    tmp_case1, ..., tmp_caseN, with 2 plots per row, in a movie.
    """
    # With skip_frames set correctly, there should be equally many
    # plot files in each directory.
    plot_files = []
    for i in range(N):
        frames = glob.glob(os.path.join(
            'tmp_case%d' % i, 'tmp_*.png'))
        frames.sort()
        plot_files.append(frames)
    num_frames = len(plot_files[0])
    # Consistency check that all cases have the same number of frames
    for i in range(1, len(plot_files)):
        if len(plot_files[i]) != num_frames:
            raise ValueError(
                'tmp_case%d has %d frames, tmp_case0 has %d'
                % (i, len(plot_files[i]), num_frames))
    combinedir = 'tmp_combined'
    if os.path.isdir(combinedir):
        shutil.rmtree(combinedir)
    os.mkdir(combinedir)
    for i in range(num_frames):
        frame_files = ' '.join([plot_files[j][i] for j in
                                range(len(plot_files))])
        # Output files must be numbered from 1 and upwards
        cmd = 'montage -background white -geometry 100%% '\
              '-tile 2x %s %s' \
              % (frame_files, os.path.join(
                 combinedir, 'tmp_%04d.png' % i))
        print cmd
        os.system(cmd)
    os.chdir(combinedir)
    cmd = 'ffmpeg -r 2 -i tmp_%04d.png -c:v flv movie.flv'
    os.system(cmd)

if __name__ == '__main__':
    N = 4
    run_simulations(N, 0.25, 30)
    make_movie(N)

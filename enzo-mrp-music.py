import os
import sys
import glob
import yt
import ConfigParser as cp
import multiprocessing as mp
from get_halo_initial_extent import *
from particle_only_mask import *
#
# Location of MUSIC executable
#
music_exe_dir = "/Users/jwise/codes/music"
#
# This name will be used to create the MUSIC configuration file and
# simulation run directory name, both of which will have the level.
#
simulation_name = "auto-wrapper"
#
# Template MUSIC configuration file.  These parameters will be
# supplemented with parameters that describe the zoom-in setup.
#
template_config = "template.conf"
#
# Original MUSIC configuration file.  If named (simulation_name)-L0,
# set to None.
#
original_config = None
#
# Base run directory.  Simulation initial conditions will be moved
# into a subdirectory called (simulation_time)-L(level)
#
simulation_run_directory = "/Users/jwise/runs/13Oct15_wrapper/"
#
# Number of cores to use with MUSIC.  If none, then use all cores.
#
num_cores = None
#
# Find the Lagrangian volume of some halo (TODO: extent to any
# selector function besides a sphere) The routine either accepts the
# radius and its units or mass in solar masses.  Two examples are
# below.  A "redshift" keyword can be given to specify the redshift
# of the target halo.  If not given, it is assumed to be in the last
# dataset that was created.
#
#halo_info = dict(center = [0.5, 0.5, 0.5],
#                 rvir = 10.0, r_units = "kpc")
#halo_info = dict(center = [0.957747,    0.13981756,  0.84212896],
#                 mass = 4.13e11)  # most massive
#halo_info = dict(center = [0.47444206,  0.48449989,  0.50042353],
#                 mass=6.05e11)
halo_info = dict(center = [0.47097584,  0.48067752,  0.5347055],
                 redshift=1.0,
                 mass=6.64e11)

#
# Safety factor to increase the radius of the sphere in units of the virial radius.
#
radius_factor = 3.0
#
# Shape type to use when calculating the Lagrangian region.  Can be
# box, ellipsoid, convex_hull, or exact.
#
shape_type = "exact"
#
########################################################################
########################### end of options #############################
########################################################################
# Obtain the next level from the command line
#
if len(sys.argv) != 2:
    raise RuntimeError("usage: %s level\n"
                       "\t level: 0-based level of the next set of ICs" % \
                       (sys.argv[0]))
level = int(sys.argv[-1])

# Error check
if level == 0:
    raise RuntimeError("level must be >0. "
                       "Please run the unigrid simulation first.")
files_to_check = ["%s/MUSIC" % (music_exe_dir),
                  template_config,
                  simulation_run_directory]
if original_config != None: files_to_check += original_config
for f in files_to_check:                  
    if not os.path.exists(f):
        raise RuntimeError("File/directory not found: %s" % (f))

# Set simulation directories
prev_sim_dir = os.path.join(simulation_run_directory, "%s-L%d" %
                            (simulation_name, level-1))
sim_dir = os.path.join(simulation_run_directory,
                       "%s-L%d" % (simulation_name, level))
#
# Obtain the maxlevel of the original run
if original_config == None:
    original_config_file = "%s-L0.conf" % (simulation_name)
else:
    original_config_file = original_config
music_cf0 = cp.ConfigParser()
music_cf0.read(original_config_file)
initial_min_level = music_cf0.getint("setup", "levelmin")
initial_max_level = music_cf0.getint("setup", "levelmax")

# Obtain the shift of the Lagrangian region from the previous zoom-in
# (or unigrid) simulation
region_shift = [0, 0, 0]
prev_config_logfile = "%s-L%d.conf_log.txt" % (simulation_name, level-1)
with open(prev_config_logfile) as fp:
    for l in fp.readlines():
        if l.find("setup/shift_x") >= 0:
            region_shift[0] = int(l.split("=")[1])
        if l.find("setup/shift_y") >= 0:
            region_shift[1] = int(l.split("=")[1])
        if l.find("setup/shift_z") >= 0:
            region_shift[2] = int(l.split("=")[1])
        if l.find("setup/levelmin") >= 0:
            region_point_levelmin = int(l.split("=")[1])

# Rounding factor for the Lagrangian region if using a rectangular
# prism.
round_factor = 2**initial_max_level

#
# Get the inital dataset of the simulation and either
# the final dataset or the dataset at the specified redshift.
#
sim_par_file = os.path.join(prev_sim_dir, "%s-L%d.enzo" %
                            (simulation_name, level-1))
es = yt.simulation(sim_par_file, "Enzo", find_outputs=True)

enzo_initial_fn = es.all_outputs[0]["filename"]
if "redshift" in halo_info:
    es.get_time_series(redshifts=[halo_info["redshift"]])
    ds = es[0]
    enzo_final_fn = os.path.join(ds.fullpath, ds.basename)
else:
    enzo_final_fn = es.all_outputs[-1]["filename"]

particle_output_format = None if shape_type == "box" else "txt"
region_center, region_size, lagr_particle_file = \
               get_center_and_extent(halo_info,
                                     enzo_initial_fn, enzo_final_fn,
                                     round_size = round_factor,
                                     radius_factor = radius_factor,
                                     output_format = particle_output_format)

#
# Read the zoom-in MUSIC file, modify/add zoom-in parameters, and write out.
#
music_cf1 = cp.ConfigParser()
# Turn-on case-sensitive for config files
music_cf1.optionxform = str

music_cf1.read(template_config)
# Delete some options if they exist.  If we need them, we'll create them again.
for option in ["ref_offset", "ref_center", "ref_extent"]:
    if music_cf1.has_option("setup", option):
        music_cf1.remove_option("setup", option)

music_cf1.set("setup", "levelmax", "%d" % (initial_min_level + level))
music_cf1.set("output", "filename", "%s-L%d" % (simulation_name, level))
music_cf1.set("setup", "region",
              "convex_hull" if shape_type == "exact" else shape_type)
if shape_type == "box":
    music_cf1.set("setup", "ref_center", "%f, %f, %f" % \
                  (region_center[0], region_center[1], region_center[2]))
    music_cf1.set("setup", "ref_extent", "%f, %f, %f" % \
                  (region_size[0], region_size[1], region_size[2]))
else:
    music_cf1.set("setup", "region_point_file", lagr_particle_file)
    music_cf1.set("setup", "region_point_shift",
                  "%d, %d, %d" % (region_shift[0], region_shift[1], region_shift[2]))
    music_cf1.set("setup", "region_point_levelmin", "%d" % (initial_min_level))
    
new_config_file = "%s-L%d.conf" % (simulation_name, level)
with open(new_config_file, "wb") as fp:
    music_cf1.write(fp)

# Set the number of OpenMP threads and run MUSIC
if num_cores == None:
    num_cores = mp.cpu_count()
os.environ["OMP_NUM_THREADS"] = "%d" % (num_cores)
os.system("%s/MUSIC %s" % (music_exe_dir, new_config_file))

# If we require the exact Lagrangian region, then we directly modify
# the RefinementMask file that's written by MUSIC.
#
# smooth_edges: further smooth the CIC interpolation of the particles
# in the Lagrangian region with a Gaussian over a 3x3x3 cell volume.
#
# backup: Copy original file with the suffix .bak
if shape_type == "exact":
    particle_only_mask(new_config_file, smooth_edges=True, backup=True)

# Modify the skeleton Enzo parameter file created by MUSIC to include
# the parameters for must-refine particles.
ic_dir = music_cf1.get("output", "filename")
fp = open("%s/parameter_file.txt" % (ic_dir), "a")
fp.write("\n"
         "#\n"
         "# must-refine particle parameters\n"
         "# *** must also include method 8 in CellFlaggingMethod ***\n"
         "# *** do NOT include the RefineRegion parameters above ***\n"
         "#\n"
         "MustRefineParticlesCreateParticles = 3\n"
         "MustRefineParticlesRefineToLevel   = %d\n"
         "CosmologySimulationParticleTypeName          = RefinementMask\n" % (level))
fp.close()

# Copy initial conditions directory to the simulation run directory
print "Moving initial conditions to %s" % (sim_dir)
os.rename(ic_dir, sim_dir)

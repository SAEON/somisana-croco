import crocotools_py.plotting as crocplt
import numpy as np

# input files
croco_blk_file = 'croco_blk_WASA3_Y2000M12.nc'
croco_grd = '/home/gfearon/code/somisana-croco/configs/algoa_01/croco_v1.3.1/GRID/croco_grd.nc'

crocplt.plot_blk(croco_grd,
                croco_blk_file,
                tstep=0, # time index in file (not going to worry about decoding the actual dates here)
                var='wspd',
                figsize=(6,6), # (hz,vt)
                ticks = np.linspace(0,15,num=16), # the ticks to plot
                cmap = 'Spectral_r',
                cbar_loc = [0.9, 0.2, 0.02, 0.6],
                add_vectors = True,
                scale_uv = 150,
                skip_uv = 10,
                jpg_out=None,
                write_jpg=False,
                gif_out=None,
                write_gif=False,
                tstep_end=None
                )

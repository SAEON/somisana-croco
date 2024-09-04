from datetime import datetime
import crocotools_py.preprocess as pre
    
wasa_grid = '/media/external_1/geo_em.d03.nc'
wasa_zarr_dir = '/media/external_1/uv_ds_10m_v2'
croco_grd = '/home/gfearon/code/somisana-croco/configs/sa_west_02/croco_v1.3.1/GRID/croco_grd.nc'
croco_blk_dir = '/media/external_1/somisana-croco/configs/sa_west_02/croco_v1.3.1/ERA5'
out_wasa_dir = '/media/external_1/somisana-croco/configs/sa_west_02/croco_v1.3.1/WASA3'
ref_date = datetime(1993,1,1)

pre.make_WASA3_from_blk(wasa_grid, 
                 wasa_zarr_dir, 
                 croco_grd,
                 croco_blk_dir,
                 out_wasa_dir,
                 ref_date,
                 croco_atmos='ERA5',
                 interp_method='linear'
                 )    

# after running make_WASA3_from_blk we found some missing values in the WASA
# data for the month of Y2010M04, so in this next step we just fill in 
# these missing data so that the model can run
# pre.fill_blk(croco_grd, out_wasa_dir + '/croco_blk_WASA3_Y2010M04_old.nc', out_wasa_dir + '/croco_blk_WASA3_Y2010M04.nc')

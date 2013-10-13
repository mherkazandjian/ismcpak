#merging all Dbs into one (particular application to grids with different amount of CR rates)
#--------------------------------------------------------------------------------------------

import meshUtils
import time

merged_db_dir_path = '/home/mher/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-sweep/'

#the base DB
arxv = meshUtils.meshArxv(dirPath = '/home/mher/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-1_solar/', readDb=True)
arxv.infoAll['parms'][:,3] = arxv.parms_used_in_PDR_models['zeta'] 


db_paths_add = [ '/home/mher/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-10_solar/',
                 '/home/mher/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-100_solar/',
                 '/home/mher/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-1000_solar/',
                 '/home/mher/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-10000_solar/',
                ]

for path in db_paths_add:
    arxv_new = meshUtils.meshArxv(dirPath = path, readDb=True)
    arxv_new.infoAll['parms'][:,3] = arxv_new.parms_used_in_PDR_models['zeta'] 
    
    arxv.mergeDbs(new_arxv=arxv_new, outDirPath=None)

arxv.writeDb(merged_db_dir_path)
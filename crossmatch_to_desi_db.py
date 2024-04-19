import pandas as pd
import psycopg2
import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
import sys

def main():
    if (len(sys.argv) != 5):
        print("crossmatch_to_desi_db.py - Crossmatch a catalog with observations from DESI.")
        print("Usage: python crossmatch_to_desi_db.py [catalag .fits filename with columns 'RA' and 'DEC'] [DESI data reduction] [search radius in arcsec] [output file]]")
        print("Example: crossmatch_to_desi_db.py /my/path/my_transients.fits daily 1.6 /my/path/output.fits")
        sys.exit(0)

    filename = sys.argv[1]
    data_release = sys.argv[2]
    radius = float(sys.argv[3])
    outfile = sys.argv[4]

    # Read in catalog
    catalog = pd.read_csv(filename)

    # connect to database
    db = psycopg2.connect(host='decatdb.lbl.gov', database='desidb', user='desi', password = "5kFibers!", port="5432")
    cursor = db.cursor()

    # Query parameters
    radius = radius / 3600.0 # arcseconds

    obs_status = []  
    z_lst = []
    zerr_lst = []
    spectype_lst = []
        
    for idx, (ra, dec) in enumerate(zip(catalog["RA"], catalog["DEC"])):
        observed = "N"
        if np.isnan(ra) or np.isnan(dec):
            match = np.empty((6,))
            match[:] = np.nan
        else:
            # Query database
            query = 'SELECT f.target_ra, f.target_dec, r.z, r.zerr, r.zwarn, r.spectype\n' \
                        f'FROM {data_release}.tiles_fibermap f\n' \
                        f'INNER JOIN {data_release}.cumulative_tiles c ON f.cumultile_id=c.id\n' \
                        f'INNER JOIN {data_release}.tiles_redshifts r ON r.cumultile_id=c.id AND r.targetid=f.targetid\n' \
                        f'WHERE q3c_radial_query(f.target_ra, f.target_dec, {ra}, {dec}, {radius});'
            cursor.execute(query)
            targets = cursor.fetchall()


            if len(targets) >= 1:
                targets = np.array(targets)
                zwarn = '1'

                while zwarn != '0' and len(targets) > 0:
                    ra_lst, dec_lst = targets[:,0].astype(dtype=np.float64).tolist(), targets[:,1].astype(dtype=np.float64).tolist()
                    match_idx = match_coordinates_sky(SkyCoord(ra << u.deg, dec << u.deg), SkyCoord(ra_lst << u.deg, dec_lst << u.deg), nthneighbor=1)[0]
                    match = targets[match_idx]
                    zwarn = match[4]
                    if zwarn != '0':
                        targets = np.delete(targets, match_idx, 0)

            if len(targets) == 0:
                match = np.empty((6,))
                match[:] = np.nan
            else:
                observed = "Y"
            
        z_lst.append(match[2])
        zerr_lst.append(match[3])
        spectype_lst.append(match[5])   
        obs_status.append(observed)

    catalog["Z"] = z_lst
    catalog["Z_ERR"] = zerr_lst
    catalog["SPECTYPE"] = spectype_lst
    catalog["OBSERVATION_STATUS"] = obs_status

    n_match = len(catalog["Z"]) - catalog["Z"].isna().sum()
    n_tot = len(catalog["Z"])
    print(f"{n_match} out of {n_tot} objects have a DESI redshift.")

    catalog.to_csv(outfile)

if __name__ == "__main__":
    main()





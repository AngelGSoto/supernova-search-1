import splusdata
import getpass
import pandas as pd

# Connecting with SPLUS database

username = str(input("Login: "))
password = getpass.getpass("Password: ")
conn = splusdata.connect(username, password)

# Reading the csv table
df = pd.read_csv('./results/candidates.csv')

# Query with our criteria nad join the tables
Query = f"""SELECT detection.ID, detection.RA, detection.DEC, r.FWHM_r, 
            u.U_auto, J0378.J0378_auto, J0395.J0395_auto, J0410.J0410_auto, J0430.J0430_auto, g.G_auto, 
            J0515.J0515_auto, r.R_auto, J0660.J0660_auto, i.I_auto, J0861.J0861_auto, z.Z_auto, u.e_U_auto, 
            J0378.e_J0378_auto, J0395.e_J0395_auto, J0410.e_J0410_auto, J0430.e_J0430_auto, g.e_G_auto, 
            J0515.e_J0515_auto, r.e_R_auto, J0660.e_J0660_auto, i.e_I_auto, J0861.e_J0861_auto, z.e_Z_auto,
            pz.zml  
            FROM TAP_UPLOAD.upload as tap 
            LEFT OUTER JOIN idr3.detection_image as detection ON tap.ID= detection.ID 
            LEFT OUTER JOIN idr3.u_band as u ON tap.ID=u.ID 
            LEFT OUTER JOIN idr3.J0378_band as J0378 ON tap.ID=J0378.ID 
            LEFT OUTER JOIN idr3.J0395_band as J0395 ON tap.ID=J0395.ID 
            LEFT OUTER JOIN idr3.J0410_band as J0410 ON tap.ID=J0410.ID 
            LEFT OUTER JOIN idr3.J0430_band as J0430 ON tap.ID=J0430.ID 
            LEFT OUTER JOIN idr3.g_band as g ON tap.ID=g.ID 
            LEFT OUTER JOIN idr3.J0515_band as J0515 ON tap.ID=J0515.ID 
            LEFT OUTER JOIN idr3.r_band as r ON tap.ID=r.ID 
            LEFT OUTER JOIN idr3.J0660_band as J0660 ON tap.ID=J0660.ID 
            LEFT OUTER JOIN idr3.i_band as i ON tap.ID=i.ID 
            LEFT OUTER JOIN idr3.J0861_band as J0861 ON tap.ID=J0861.ID 
            LEFT OUTER JOIN idr3.z_band as z ON tap.ID=z.ID
            JOIN "idr3_vacs"."photoz" AS pz ON tap.ID=pz.ID"""

# Applying query
result = conn.query(Query, df)

# Converting the astropy table into pandas and saving
df_result = result.to_pandas()
df_result.to_csv('./results/candidates-vacs.csv', index=False)

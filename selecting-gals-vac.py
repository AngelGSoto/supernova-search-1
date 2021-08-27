import splusdata
import getpass
import pandas as pd

# Connecting with SPLUS database

username = input(prompt="Login: ")
password = getpass.getpass("Password: ")
conn = splusdata.connect(username, password)

# Reading the csv table
df = pd.read_csv('selected-gals.csv')
df.head(3)

# Query with our criteria nad join the tables
Query = f"""SELECT upl.ID, upl.RA, upl.DEC, upl.FWHM_R, sgq.PROB_STAR, sgq.PROB_QSO, sgq.PROB_GAL 
                 FROM TAP_UPLOAD.upload AS upl LEFT OUTER JOIN "idr3_vacs"."star_galaxy_quasar" AS sgq 
                 ON upl.id = sgq.id WHERE "CLASS" = 2"""

# Applien query
result = conn.query(Query, df)

#Converting the astropy table into pandas and saving
df_result = result.to_pandas()
df_result.to_csv("selected-gals-class-liliane.csv")
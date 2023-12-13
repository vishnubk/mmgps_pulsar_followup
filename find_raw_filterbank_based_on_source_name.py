import numpy as np
import pandas as pd 
from sqlalchemy import create_engine
import sys
import argparse
import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Postgres username, password, and database name
POSTGRES_ADDRESS = os.getenv("DB_HOST") ## INSERT YOUR DB ADDRESS IF IT'S NOT ON PANOPLY
POSTGRES_PORT = os.getenv("DB_PORT")
POSTGRES_USERNAME = os.getenv("DB_USERNAME")
POSTGRES_PASSWORD = os.getenv("DB_PASSWORD") ## CHANGE THIS TO YOUR PANOPLY/POSTGRES PASSWORD
POSTGRES_DBNAME = 'trapum_web' ## 
# A long string that contains the necessary Postgres login information
postgres_str = ('mysql+pymysql://{username}:{password}@{ipaddress}:{port}/{dbname}'
.format(username=POSTGRES_USERNAME,
password=POSTGRES_PASSWORD,
ipaddress=POSTGRES_ADDRESS,
port=POSTGRES_PORT,
dbname=POSTGRES_DBNAME))
# Create the connection
cnx = create_engine(postgres_str)

parser = argparse.ArgumentParser(description='Extract Raw fil Files for a given Source name and beam number')
parser.add_argument('-b', '--beam', help='beam name. for eg cfbf00000', type=str)
parser.add_argument('-t', '--source_name', help='MGPS Source Name', default='MSGPS_L_0243', type=str)

args = parser.parse_args()


source_name = args.source_name
file_type_pattern = 'filterbank-raw%'

if args.beam:
    beam_name = args.beam
    query_get_filterbank_files = '''
    SELECT DISTINCT p.id AS pointing_id, t.source_name as target, b.id AS beam_id, b.name AS beam_num, p.utc_start, 
    concat(dp.filepath, "/", dp.filename) AS filterbank_files, ft.name AS file_type
    FROM data_product dp
    JOIN pointing p ON dp.pointing_id = p.id
    JOIN beam b ON dp.beam_id = b.id
    JOIN target t ON p.target_id = t.id
    JOIN file_type ft ON dp.file_type_id = ft.id
    WHERE t.source_name LIKE %s AND b.name LIKE %s AND ft.name LIKE %s order by p.utc_start desc
    '''

    printable_query = query_get_filterbank_files % (repr(source_name), repr(beam_name), repr(file_type_pattern))
   
    print(printable_query)
    df = pd.read_sql_query(query_get_filterbank_files, con=cnx, params=(source_name, beam_name, file_type_pattern))
  
else:
    query_get_filterbank_files = '''
    SELECT DISTINCT p.id AS pointing_id, t.source_name as target, b.id AS beam_id, b.name AS beam_num, p.utc_start, 
    concat(dp.filepath, "/", dp.filename) AS filterbank_files, ft.name AS file_type
    FROM data_product dp
    JOIN pointing p ON dp.pointing_id = p.id
    JOIN beam b ON dp.beam_id = b.id
    JOIN target t ON p.target_id = t.id
    JOIN file_type ft ON dp.file_type_id = ft.id
    WHERE t.source_name LIKE %s AND ft.name LIKE %s order by p.utc_start desc, b.name asc
    '''

    printable_query = query_get_filterbank_files % (repr(source_name), repr(file_type_pattern))
    print(printable_query)
    df = pd.read_sql_query(query_get_filterbank_files, con=cnx, params=(source_name, file_type_pattern))
   

df.to_csv('filterbank_files.csv', index=False)
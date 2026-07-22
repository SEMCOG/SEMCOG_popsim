import os

import pandas as pd
import yaml
from sqlalchemy.orm import close_all_sessions

CREDENTIALS_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "configs",
    "db_connection.yaml",
)


def get_db_connection_str():
    """ Resolve the database connection string

  Checked in order: the POPSIM_DB_URL environment variable, then the gitignored
  configs/db_connection.yaml. The credential is deliberately kept out of
  configs/sql.yaml, which is tracked in a public repository.

  Returns:
    Connection string
  """
    from_env = os.environ.get("POPSIM_DB_URL")
    if from_env:
        return from_env
    if os.path.exists(CREDENTIALS_PATH):
        with open(CREDENTIALS_PATH) as f:
            conn = (yaml.safe_load(f) or {}).get("db_connection_str")
        if conn:
            return conn
    raise RuntimeError(
        "No database connection string found. Set POPSIM_DB_URL, or copy "
        "configs/db_connection.yaml.example to configs/db_connection.yaml and "
        "fill it in. Do not put the credential in configs/sql.yaml - that file "
        "is tracked in a public repository."
    )


def list_tables(sql_config):
    """ Get a list of tables in the database
  Args:
    sql_config:   sql config object
  
  Returns:
    Table names in Pandas Series
  """
    table = pd.read_sql(sql_config["sql_list_tables"], get_db_connection_str())
    return table["tablename"]


def load_from_sql(sql_config, hdf):
    """ Load tables from sql to target hdf
  Args: 
    sql_config:        list of sql config items
    hdf:               target hdf object
  
  Returns: 
    hdf
  """
    sql_tables_to_load = sql_config["sql_tables_to_load"]
    for sql_table_config in sql_tables_to_load:
        table_name = sql_table_config["name"]
        print("sql table:", table_name)
        table_sql = sql_table_config["sql"]
        table_index_col = sql_table_config["index_col"]
        pd.read_sql(
            table_sql, get_db_connection_str(), index_col=table_index_col
        ).to_hdf(hdf, table_name)
    close_all_sessions()
    return hdf


def load_pop_synthetic_csv(df):
    df = df.rename(columns={
        "household_id": "household_id",
        "PUMA": "PUMA",
        "TRACT": "tract",
        "BLKGRP": "block group",
        "hh_id": "hh_id",
        "NP": "persons",
        "VEH": "cars",
        "HINCP": "HINCP",
        "R18": "children",
        "AGEHOH": "age_of_head",
        "HRACE": "race_id",
        "HHISP": "hhisp",
        "HWORKERS": "workers",
        "ADJINC": "adjinc",
        "income": "income",
        "NOC": "NOC",
        "TYPE": "type",
        "YBL": "ybl",
        "BLD": "bld",
        "VALP": "valp",
        "GRNTP": "rent",
        "ADJHSG": "ADJHSG",
        "TEN": "tenure",
        "HHT": "hht",
    })
    return df.set_index("household_id")


def load_from_files(files_config, hdf):
    """ Load tables from files to target hdf
  Args: 
    files_config:      list of sql config items
    hdf:               target hdf object
  
  Returns: 
    hdf
  """
    files_to_load = files_config["csv_files_to_load"]
    for file_table_config in files_to_load:
        table_name = file_table_config["name"]
        table_path = file_table_config["path"]
        # optional
        table_kwargs = file_table_config.get("kwargs", {})
        # optional
        table_col_rename_mapping = file_table_config.get("rename", {})
        try:
            excel = file_table_config["excel"]
        except KeyError:
            excel = False

        if excel == True:
            df = pd.read_excel(table_path, **table_kwargs)
        else:
            df = pd.read_csv(table_path, **table_kwargs)

        df.rename(columns=table_col_rename_mapping)
        post_processing_steps = file_table_config.get("post_processing_steps", [])
        for post_processing_step in post_processing_steps:
            eval(post_processing_step)
        df.to_hdf(hdf, table_name)
        print("Finishing loading %s to hdf" % (table_name))
    return hdf


import pandas as pd
from sqlalchemy.orm import close_all_sessions


def list_tables(sql_config):
    """ Get a list of tables in the database
  Args:
    sql_config:   sql config object
  
  Returns:
    Table names in Pandas Series
  """
    table = pd.read_sql(sql_config["sql_list_tables"], sql_config["db_connection_str"])
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
            table_sql, sql_config["db_connection_str"], index_col=table_index_col
        ).to_hdf(hdf, table_name)
    close_all_sessions()
    return hdf


def load_pop_synthetic_csv(df):
    df.columns = [
        "household_id",
        "PUMA",
        "tract",
        "block group",
        "hh_id",
        "persons",
        "cars",
        "HINCP",
        "children",
        "age_of_head",
        "race_id",
        "hhisp",
        "workers",
        "adjinc",
        "income",
        "NOC",
        "type",
        "ybl",
        "bld",
        "valp",
        "rent",
        "ADJHSG",
        "tenure",
        "hht",
    ]
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


# hdf related functions/decorators 
  
def add_table(hdf, name):
  def inner_decorator(f):
    def wrapped(*args, **kwargs):
      table = f(*args, **kwargs) 
      table.to_hdf(hdf, name) 
      return table 
    return wrapped
  return inner_decorator

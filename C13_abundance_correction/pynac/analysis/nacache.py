# Custom Packages
from pynac.core.nalookup import *
from pynac.core          import isotopelabels as IsotopeLabels

class NACache(object):
  def __init__(self, isotopes):
    self.Isotopes   = isotopes
    self.Tables     = {}
  
  def _check2(self, isotope, max):
    if(not isinstance(isotope, str)):
      return IndexError("first index must be a string specifying an isotope symbol")
    elif(not IsotopeLabels.IsSupported(isotope)):
      return IndexError("cannot use unsupported isotope symbol '"+isotope+"' as index")
    elif(not self.Isotopes.count(isotope)):
      return IndexError("table not indexed for isotope symbol '"+isotope+"'")
    elif(not isinstance(max, int) and not isinstance(max, float)):
      return IndexError("isotope maximum must be an integer")
    elif(max < 1):
      return IndexError("isotope maximum must be greater than or equal to 1")
    return None 
    
  def __getitem__(self, index):
    error = self._check2(index[0], index[1])
    if(error): raise error
    key = index[0] + str(index[1])
    try: return self.Tables[key]
    except Exception as exception: raise exception
  
  def GetNALookup(self, maximums):
    tables = {}
    for isotope,max in maximums.iteritems():
      tables[isotope] = self[isotope, max]
    return NALookup(tables)

  def BuildTable(self, isotope, max): 
    # TODO: throw custom exception if status != initialized
    # TODO: throw custom exception if status == finalized
    error = self._check2(isotope, max)
    if(error): raise error
    key = isotope + str(max)
    if(not self.Tables.has_key(key)):
      na = IsotopeLabels.NA[isotope]
      self.Tables[key] = NABCTables(max, na)

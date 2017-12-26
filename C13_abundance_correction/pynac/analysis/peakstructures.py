#!/usr/bin/python
# Standard Packages
import re
import math
# Custom Packages
from pynac.core import isotopelabels as IsotopeLabels

class PeakSet(object):
  PTRegEx = re.compile("(?P<percent>0*.\d+)(?P<statistic>MAX|MIN|AVG)")
  def __init__(self, maximums, formula):
    self.__maximums  = maximums
    self.__order = []
    for m in maximums: 
      if(maximums[m] > 0): 
        self.__order.append(m)
    self.__order     = tuple(self.__order)
    self.__shape     = tuple(map(lambda i: maximums[i]+1, self.__order))
    self.__peaks     = {}
    self.__formula   = formula
    self.__rawdata   = {}
    self.__threshold = float('inf')
    self.__minpeak   = None
    self.__maxpeak   = None
    self.__sum       = 0
  
  def __len__(self):
    return len(self.__peaks)
  
  @property
  def Formula(self):
    return self.__formula
  
  @property
  def Shape(self):
    return self.__shape
    
  @property
  def RawSet(self):
    return self.__peaks
  
  @property
  def TotalPeaks(self):
    return len(self.__peaks)
  
  @property
  def LabelOrder(self):
    return self.__order
  
  @property
  def RawData(self):
    return self.__rawdata
    
  @property
  def MaximumPeak(self):
    return self.__minpeak
  
  @property
  def MinimumPeak(self):
    return self.__maxpeak

  @property
  def IsotopeMaximums(self):
    return self.__maximums.copy()
    
  @property
  def Threshold(self):
    if(isinstance(self.__threshold, float) or
       isinstance(self.__threshold, int)):
      return self.__threshold
    match = PeakSet.PTRegEx.match(self.__threshold)  
    if(match):
      percent   = float(match.group("percent"))
      statistic = match.group("statistic")
      if(statistic == "MAX"):
        return percent*self.MaximumPeak.Intensity
      elif(statistic == "MIN"):
        return percent*self.MinimumPeak.Intensity
      elif(statistic == "AVG"):
        return percent*self.TotalIntensity/len(self.Peaks)
    else:
      return float('inf')

  @Threshold.setter
  def Threshold(self, value):
    # TODO: Input hardening (string validation, type checking)
    self.__threshold = value
  
  def __checkpeak(self, peak, index):
    #Peak already in set, must be a duplicate
    if(self.__peaks.has_key(index)):
      e1msg = "duplicated: peak set already contains peak for isotope counts "
      for isotope, count in peak.IsotopeCounts.iteritems():
        e1msg += "{0} {1}, ".format(count, IsotopeLabels.Name[isotope])
      return PeakError(e1msg.strip(','))
    #Isotope count errors
    for isotope in self.__order:
      name    = IsotopeLabels.Name[isotope]
      element = IsotopeLabels.Element[isotope]
      #No count for expected isotope label 
      if(isotope not in peak.IsotopeCounts):
        e2msg = "peak contains no isotope count for {0}"
        return PeakError(e2msg.format(name))
      #Count for a label greater than molecule's max for label's element
      if(peak.IsotopeCounts[isotope] > self.__maximums[isotope]):
        e3msg = "isotope count for {0} above molecule's count for {1}" 
        return PeakError(e3msg.format(name,element))
    return False
  
  def __buildindex(self, peak):
    index = []
    for isotope in self.__order:
      index.append(peak.IsotopeCounts[isotope])
    return tuple(index)
  
  def IsotopeMaximum(self, isotope):
    return self.__maximums[isotope]
  
  def AddPeak(self, peak):
    index = self.__buildindex(peak)
    error = self.__checkpeak(peak,index)
    if(not error):
      self.__peaks[index]   = peak
      self.__rawdata[index] = peak.Intensity
      if(peak.Intensity != 0):
        self.__sum += peak.Intensity
        if(self.MinimumPeak == None or peak < self.MinimumPeak):
          self.__minpeak = peak
        if(self.MaximumPeak == None or peak > self.MaximumPeak):
          self.__maxpeak = peak   
    else: raise error
    
class Peak(object):
  def __init__(self, intensity, isocounts):
    self.__intensity = intensity
    self.__isotopecounts = isocounts  
    self.Comment = ""
    
  @property
  def Intensity(self):
    return self.__intensity
  
  @property
  def IsotopeCounts(self):
    return self.__isotopecounts
 
  def GetIsotopeCount(self,isotope):
    try            : return self.__isotopecounts[isotope]
    except KeyError: return None

  def copy(self):
    copy = Peak(self.Intensity, self.IsotopeCounts)
    copy.Comment = self.Comment
    return copy

  def __lt__(self, other):
    if(isinstance(other, Peak)):
      return self.Intensity < other.Intensity
    else:
      raise TypeError("unsupported operand type(s) for <=: 'Peak' and '" + type(other) + "'")
  
  def __le__(self, other):
    if(isinstance(other, Peak)):
      return self.Intensity <= other.Intensity
    else:
      raise TypeError("unsupported operand type(s) for <=: 'Peak' and '" + type(other) + "'")

  def __ge__(self, other):
    if(isinstance(other, Peak)):
      return self.Intensity > other.Intensity
    else:
      raise TypeError("unsupported operand type(s) for >: 'Peak' and '" + type(other) + "'")
  
  def __gt__(self, other):
    if(isinstance(other, Peak)):
      return self.Intensity >= other.Intensity
    else:
      raise TypeError("unsupported operand type(s) for >=: 'Peak' and '" + type(other) + "'")

  def __str__(self):
    return "{0} : {1}".format(str(self.IsotopeCounts),str(self.Intensity))
 
class PeakError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

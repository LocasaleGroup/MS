# Standard Packages
import ConfigParser as cp
import re, os, sys, numpy
# Custom Packages 
from pynac.core           import isotopelabels as IsotopeLabels  
from pynac.core.nalookup  import *
from pynac.core.nacmath   import *
from peakstructures       import *
from nacache              import * 
from datasource           import *
############################################################
# Default output structures
ALL_OUTPUT = {
    # Analysis Results
    "TotalNARemoved"        : "Total NA Removed",
    "Predicted"             : "Predicted",
    "Corrected"             : "Corrected",
    "Renormalized"          : "Renormalized",
    # Error and warning information
    "Comments"              : "Comments",
    # Algorithm Performance Information
    "Iterations"            : "Iterations",
    "SampleGroupKey"        : "Sample Group Key"}
PERFORMANCE_OUTPUT = {
    # Algorithm Performance Information
    "Iterations"            : "Iterations",
    "SampleGroupKey"        : "Sample Group Key"}
RESULTS_OUTPUT = {
    # Analysis Results
    "Predicted"             : "Predicted",
    "Corrected"             : "Corrected",
    "Renormalized"          : "Renormalized",
    "TotalNARemoved"        : "Total NA Removed",
    # Error and warning information
    "Comments"              : "Comments"}
DEFAULT_OUTPUT = RESULTS_OUTPUT
############################################################

class DataExtractor(object):
  def ConstructGroupKey(self, data, grouping):
    key = []
    grouping = tuple(grouping) #self.Configuration["GroupColumns"]
    for column in grouping:
      key.append(str(data[column]))
    return tuple(key)
    
  def ConstructPeak(self, data, intcol, isocols):
    counts = {}
    #isocols = self.Configuration["IsotopeColumns"]
    #intcol  = self.Configuration["IntensityColumn"]
    for isotope,column in isocols.iteritems():
      counts[isotope] = int(data[column])
    intensity = float(data[intcol])
    return Peak(intensity, counts)

  def ConstructPeakSet(self, data, isotopes, fcolumn):
    maximums = {}
    #isotopes = self.Configuration["IsotopeColumns"].keys()
    #fcolumn  = self.Configuration["FormulaColumn"]
    for isotope in isotopes:
      formula = data[fcolumn]
      maximum = IsotopeLabels.Maximum(isotope, formula)
      if(maximum != 0):
        maximums[isotope] = IsotopeLabels.Maximum(isotope, formula)
    return PeakSet(maximums,formula)
  
  def DataSourceMagic(self, filepath):
    basename,extension = os.path.splitext(filepath)
    if(not os.path.exists(filepath)):
      raise RuntimeError("Data file '{0}' does not exist".format(filepath))
    if(extension.lower() == ".csv"):
      return CSVDataSource(filepath)
    else:
      raise RuntimeError("Unsupported file type '{0}'".format(extension))

class Analysis(DataExtractor):
  PTRegEx = re.compile("(?P<percent>0*.\d+)(?P<statistic>GMAX|GMIN|GAVG|GAMIN|GAMAX)")
  def __init__(self):
    # Global Statistics
    self.GlobalMinimum = None
    self.GlobalMaximum = None
    self.GlobalSum     = 0
    self.GoodPeakCount = 0
    # Related files paths    
    self.DataSource    = None  
    # Analysis configuration 
    self.Configuration = None
    self.Initialized   = False
    # Analysis data structures
    self.BadPeaks      = []      
    self.PeakSets      = {}      
    self.PeakRowMap    = {}      
    self.Results       = None      
    self.NACache       = None
    self.SetErrors     = []

  def ConfigureFromFile(self, file):
    s = "PyNACOptions"
    c = {}
    p = cp.ConfigParser()
    p.read(file)
    ##########################################################
    # Required options
    try:
      c["FormulaColumn"  ] = int(p.get(s,"FormulaColumn"))
      c["GroupColumns"   ] = eval(p.get(s,"GroupColumns"))
      c["IsotopeColumns" ] = eval(p.get(s,"IsotopeColumns"))
      c["IntensityColumn"] = int(p.get(s,"IntensityColumn"))
    except ValueError:
      raise RuntimeError("Required option is invalid")
    except cp.NoOptionError, cp.NoSectionError:
      raise RuntimeError("Missing required option in configuration file")
    ###########################################################
    # Options that are not required
    try    : c["SkipRows"     ] = eval(p.get(s,"SkipRows"))
    except : pass
    try    : c["Threshold"    ] = str(p.get(s,"Threshold"))
    except : pass
    try    : c["DataFile"     ] = str(p.get(s,"DataFile"))
    except : pass
    try    : c["OutputFile"   ] = str(p.get(s,"OutputFile"))
    except : pass
    try    : c["Output"       ] = eval(p.get(s,"Output"))
    except : pass
    ###########################################################
    self.Configuration = self.__fillconfdefaults(c)
    self.Configuration["Filled"] = True

  def Initialize(self):
    if(not self.Configuration):
      raise RuntimeError("Analysis not configured")
    else:
      if(not self.Configuration.has_key("Filled")):
        self.Configuration = self.__fillconfdefaults(self.Configuration)
        
    if(not self.DataSource):
      try:
        self.DataSource = self.DataSourceMagic(self.Configuration["DataFile"])
      except KeyError: 
        raise RuntimeError("No data source has been specified")
        
    isotopes  = self.Configuration["IsotopeColumns"].keys()
    threshold = self.Configuration["Threshold"]
    self.Results = {}

    # Read in all the data, build all the NA lookups we will need
    self.NACache = NACache(isotopes)
    self.DataSource.ReadAll(subscriber=self)
    # Initialization done, now apply threshold to peak sets
    self.Initialized = True
    threshold = self.__resolvethreshold(threshold)
    for set in self.PeakSets.itervalues():
      set.Threshold = threshold

  def CorrectAll(self):
    if(not self.Initialized):
      raise RuntimeError("Analysis has not been initialized")
    # Enqueue all sets to run the correction analysis
    for key,set in self.PeakSets.iteritems():
      if(len(set.Shape) == 0):
        emsg = "Bad set with key {0}".format(str(key))
        sys.stderr.write(emsg+"\n")
        self.SetErrors.append(emsg)
        continue
      if(len(set.RawData) == 0):
        emsg = "Zero length set with key {0}".format(str(key))
        sys.stderr.write(emsg+"\n")
        self.SetErrors.append(emsg)
        continue
      self.Results[key] = self.CorrectNA(set)

  def CorrectNA(self, set):
    if(not self.Initialized):
      raise RuntimeError("Analysis has not been initialized")
    # Build the appropriate NATables from the cache
    lookup = {}
    ptables = []
    stables = []
    for axis in xrange(len(set.LabelOrder)):
      isotope = set.LabelOrder[axis]
      # Don't add tables if this set does not incorporate the isotope
      if set.IsotopeMaximum(isotope) > 0:
        tables = self.NACache[isotope, set.IsotopeMaximum(isotope)]
        ptables.append(tables[0]) # Product Lookups
        stables.append(tables[1]) # Sum Product Lookups
    dimensions = len(set.Shape)
    if(dimensions == 1): 
      data = [0]*(set.Shape[0])
      slookup = stables[0]
      plookup = ptables[0]
    else:
      data = numpy.zeros(shape=set.Shape, dtype=float)
      slookup = NASumProduct(stables)
      plookup = NAProduct(ptables)
    for coordinates in set.RawData:
      if(dimensions == 1): 
        data[coordinates[0]] = set.RawData[coordinates]
      else: 
        data[coordinates] = set.RawData[coordinates]
    analysis = NACorrector(set.LabelOrder, set.Shape, slookup, plookup)
    analysis.Run(data)
    return analysis.FilterFormatResults(set.Threshold)

  def WriteOutput(self, filepath=None):
    if(not self.Configuration):
      raise RuntimeError("Analysis not configured")
    if(not self.Results):
      raise RuntimeError("No results to output")
    if(filepath == None and not self.Configuration.has_key("OutputFile")):
      raise RuntimeError("No output file specified")
    elif(filepath == None):
      filepath = self.Configuration["OutputFile"]
    # Add headers and get their column numbers
    columnnums = self.__fillheaders()
    #########################################################################
    # Iterate through set results and write
    for setkey in self.Results:
      setresults  = self.Results[setkey]
      peakresults = self.Results[setkey]["#PeakResults"]
      set         = self.PeakSets[setkey]
      #########################################################################
      # Iterate through peak results and write
      for coordinates in peakresults:
        observed = set.RawSet.has_key(coordinates)
        if(not observed):
          comment = "Note: A value for this peak was predicted but not found in the original data"
          rownum  = self.DataSource.MaxRow+1
          self.__fillunobserved(rownum, coordinates, setkey)
        else:
          peak    = set.RawSet[coordinates]
          comment = peak.Comment
          rownum  = self.PeakRowMap[peak]
        #########################################################################
        # Special output options
        if(columnnums.has_key("SampleGroupKey")):
          self.DataSource.AddCell(rownum, columnnums["SampleGroupKey"], setkey)
        if(columnnums.has_key("Comments")):
          self.DataSource.AddCell(rownum, columnnums["Comments"], comment)
        #########################################################################
        # Regular peak level output options
        for option in peakresults[coordinates]:
          if(columnnums.has_key(option)):
            self.DataSource.AddCell(rownum, columnnums[option], peakresults[coordinates][option])
        #########################################################################
        # Regular set level output options
        for option in setresults:
          # All keys begining with '#' in front are skipped
          if(option[0] == "#"): continue
          if(columnnums.has_key(option)):
            self.DataSource.AddCell(rownum, columnnums[option], setresults[option])
      #########################################################################
      # Iterate through bad peaks and write their information
      for peak in self.BadPeaks:
        rownum = self.PeakRowMap[peak]
        comment = peak.Comment
        if(columnnums.has_key("SampleGroupKey")):
          self.DataSource.AddCell(rownum, columnnums["SampleGroupKey"], setkey)
        if(columnnums.has_key("Comments")):
          self.DataSource.AddCell(rownum, columnnums["Comments"], comment)
    #########################################################################
    # Finish and write everything to the file specified by filepath
    self.DataSource.Write(filepath=filepath)

  def __fillconfdefaults(self,conf): 
    if(not conf.has_key("SkipRows")):
      conf["SkipRows"     ] = [0]
    if(not conf.has_key("Threshold")):
      conf["Threshold"    ] = "0.99MIN"
    if(not conf.has_key("Output")):
      conf["Output"       ] = DEFAULT_OUTPUT 
    return conf 

  def __fillheaders(self):
    options = self.Configuration["Output"]
    columnnums = {}
    for option in options:
      columnname = options[option]
      columnnums[option] = self.DataSource.MaxColumn+1
      self.DataSource.AddCell(0, columnnums[option], columnname)
    return columnnums  

  def __fillunobserved(self, rownum, coordinates, setkey):
    set = self.PeakSets[setkey]
    # Write grouping information
    grouping = tuple(self.Configuration["GroupColumns"])
    for i in xrange(len(grouping)):
      colnum = grouping[i]
      colval = setkey[i]
      self.DataSource.AddCell(rownum, colnum, colval)
    # Fill in intensity column
    colnum = self.Configuration["IntensityColumn"]
    colval = "not observed"
    self.DataSource.AddCell(rownum, colnum, colval)
    # Fill in formula
    colnum = self.Configuration["FormulaColumn"]
    colval = set.Formula
    self.DataSource.AddCell(rownum, colnum, colval)
    # Fill in isotope information
    for i in range(len(set.LabelOrder)):
      isotope = set.LabelOrder[i]
      colnum  = self.Configuration["IsotopeColumns"][isotope]
      colval  = coordinates[i]
      self.DataSource.AddCell(rownum, colnum, colval)

  def __resolvethreshold(self, threshold):
    if(not self.Initialized):
      raise RuntimeError("Analysis has not been initialized")
    match = Analysis.PTRegEx.match(str(threshold))
    # TODO: Do better threshold checking using PeakSet.PTRegEx
    # Global Statistic Threshold
    if(match):
      percent   = float(match.group("percent"))
      statistic = match.group("statistic")
      sets      = 0
      minsum    = 0
      maxsum    = 0
      if(statistic == "GAMAX" or statistic == "GAMIN"):  
        for set in self.PeakSets.itervalues():
          sets   += 1
          minsum += set.MinimumPeak.Intensity
          maxsum += set.MaximumPeak.Intensity
      if(statistic == "GMAX"):
        rv = percent*self.GlobalMaximum.Intensity
      elif(statistic == "GMIN"):
        rv = percent*self.GlobalMinimum.Intensity    
      elif(statistic == "GAVG"):
        rv = percent*self.GlobalSum/self.GoodPeakCount
      elif(statistic == "GAMAX"):
        rv = percent*maxsum/sets
      elif(statistic == "GAMIN"):
        rv = percent*minsum/sets
    # Flat Value or Local Statistic Thresholds
    else: 
      rv = threshold
    return rv

  def OnRowRead(self, rown, data):
    skip = self.Configuration["SkipRows"]
    # Don't bother doing anything if we are skipping this row
    if(rown in skip):
      return
    ##########################################################
    # Configuration data
    fcol   = self.Configuration["FormulaColumn"  ]
    isos   = self.Configuration["IsotopeColumns" ].keys()
    gcols  = self.Configuration["GroupColumns"   ]
    iscols = self.Configuration["IsotopeColumns" ]
    intcol = self.Configuration["IntensityColumn"]
    ##########################################################
    key   = self.ConstructGroupKey(data,gcols)
    peak  = self.ConstructPeak(data,intcol,iscols)
    # This peak set is new, add it
    if(not self.PeakSets.has_key(key)):
      set = self.ConstructPeakSet(data,isos, fcol)
      self.PeakSets[key] = set
      set.Threshold = self.Configuration["Threshold"]
      # Build NA Look up tables in parallel 
      for isotope, max in set.IsotopeMaximums.iteritems():
        if(max > 0): 
          self.NACache.BuildTable(isotope, max)
    # Try to add peak, record peak as bad if there is a peak error
    try:
      self.PeakSets[key].AddPeak(peak)
      # Collect global statistics if this peak added with no errors
      self.GlobalSum += peak.Intensity
      self.GoodPeakCount += 1
      if(self.GlobalMaximum == None):
        self.GlobalMaximum = peak
      elif(self.GlobalMaximum.Intensity < peak.Intensity):
        self.GlobalMaximum = peak
      if(self.GlobalMinimum == None):
        self.GlobalMinimum = peak
      elif(self.GlobalMinimum.Intensity > peak.Intensity):
        self.GlobalMinimum = peak
    except PeakError as p:
      peak.Comment = "Error: " + p.value
      self.BadPeaks.append(peak)
    finally:
      self.PeakRowMap[peak] = rown

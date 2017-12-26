#!/usr/bin/python
# Standard Packages
from   time  import clock
import numpy
# Custom Packages
import isotopelabels as IsotopeLabels
import naciter       as NI 
from   nalookup      import NASumProduct,NAProduct

class NACorrector(object):
  """
  NACorrector Object
  ==================
    Description
    -----------
    Corrects for the natural abudance of isotopes used as labels in stable 
    isotope resolved mass spectrometry experiments
    
    Example
    -------
    # Example of setup to correct the natural abudance in samples
    # of C17H27N3O17P2 where carbon-13, nitrogen-15 and deuterium were used as
    # labeling isotopes. It is assumed that "my_data" has been intialized. If
    # the correction algorithm is configured for a single-label only then
    # "my_data" can be a standard python list rather than a numpy ndarray.
    #
    # To see what isotope labels are supported please look at the
    # pynac.core.istopelabels module which contains various functions and
    # constants used for generating the NASumProduct and NAProduct tables
    # required by this algorithm.
    
    >>> my_shape = (  18,   4,   28 ) # For just carbon-13 use  (18,)
    >>> my_order = ("13C","15N","2H") # For just carbon-13 use ("13C",)
    >>> my_analysis = NACorrector(my_order, my_shape)
    >>> my_analaysis.Run(my_data) # my_data must have the same shape as my_shape
    >>> my_results_dictionary = my_analysis.FilterFormatResults()    
  """
  def __init__(self, order, shape, slookup=None, plookup=None):
    """
    NACorrector Constructor
    =======================
      Description
      -----------
      Initializes a Multi-Labeled NA Correction algorithm.
      Note that if length of the parameter "order" does not equal the length of
      the parameter "shape" this function will fail.      
            
      Arguments
      ---------
        - order : tuple of strings specifying the isotopes used to label the
                  data this algorithm will operate on
        - shape : tuple of integers specifying the shape of the data this 
                  algorithm will operate on, this shape must be the maximum
                  number of atoms of a labeling isotope that can appear in the
                  molecule intensities being corrected by this algorithm.  For
                  example, if carbon-13 is the labeling isotope, and the molecule
                  being examined has 30 carbons then the proper order and shape
                  would be ("13C",) and (30+1,)                  
      
      Keyword Arguments
      -----------------
        - slookup : NASumProduct object to be used by this correction algorithm.
                    If it is not supplied it will be generated based on the
                    arguments "order" and "shape".
        - plookup : NAProduct object to be used by this correction algorithm.
                    If it is not supplied it will be generated based on the
                    arguments "order" and "shape".                    
    """
    assert(len(order) == len(shape))
    self.__dims  = len(shape)
    self.__order = order
    if(slookup is None):
      maximums = {}
      for i in xrange(len(order)): maximums[order[i]] = shape[i] - 1
      slookup  = NASumProduct()
      slookup.BuildTables(order,maximums)
      if(self.__dims == 1): slookup = slookup._tables[0]
    if(plookup is None):
      maximums = {}
      for i in xrange(len(order)): maximums[order[i]] = shape[i] - 1
      plookup  = NAProduct()
      plookup.BuildTables(order,maximums)
      if(self.__dims == 1): plookup = plookup._tables[0]
      
    self.__slookup   = slookup
    self.__plookup   = plookup
    self.__results   = None
    # Initialize critical methods and members based on dimensionality
    if(self.__dims == 1): 
      self.__shape     = shape[0]
      self.__maximums  = shape[0] - 1
      self.__minimums  = -1 
      self.__ynacrange = xrange
      self.__xnacrange = xrange
      self.__sum       = sum
      self.__sumdiff   = lambda a,b:sum(map(lambda i: abs(a[i]-b[i]),xrange(len(a)))) 
      self.__copy      = lambda arg:arg[:]
    else:
      self.__shape     = shape
      self.__maximums  = tuple(map(lambda c: c - 1, shape))
      self.__minimums  = self.__dims*(-1,)
      self.__ynacrange = NI.UltimateNACIter
      self.__xnacrange = NI.PenultimateNACIter
      self.__sum       = lambda arg:arg.sum()
      self.__sumdiff   = lambda a,b:(abs(a - b)).sum()
      self.__copy      = lambda arg:arg.copy()

  @property
  def RawResults(self):
    """
    RawResults Property
    ===================
      Description
      -----------
      Reference to the last set of results obtained by calling the Run function
      of this object.  It is recommened that you use FilterFormatResults to
      obtain a more useful output structure      
    """   
    return self.__results
  
  @property
  def Shape(self):
    """
    Shape Property :
    ==============
      Description
      -----------
      This property is a tuple of integers specifying the shape of the data
      that will be corrected with this object. Set by constructor.      
    """
    return self.__shape
    
  @property
  def Order(self):
    """
    Order Property
    ==============
      Description
      -----------
      This property is a tuple of strings specifying which isotopes
      correspond to which axes in this object's Shape property. 
    """
    return self.__order
    
  def Run(self, data):
    """
    Run Function
    ============
      Description
      -----------
      Runs this natural abundance correction algorithm on the parameter "data"
      If "data" is not the same shape (or the same length) as 
      
      Arguments
      ---------
        - data : numpy ndarray
          A numpy n-dimensional array with the same shape as this object
          containg intensity values for isotopologues specified in this
          object's Order property.
          
      Returns
      -------
        A dictionary containg results of the analysis with the following
        structure:
        {
          "OriginalData"    : The original input data (ndarray or list)
          "Iterations"      : The number of iterations the correction
                              algorithm ran (integer)
          "CPUTime"         : The CPU Time consumed by this analysis (float)
          "Corrected"       : The corrected values (ndarray or list)
          "TotalNARemoved"  : The total of the data after it was corrected for
                              natural abundance, ie the sum of "Corrected" 
                              (float)
          "Predicted"       : The predicted peak intensities(ndarray or list)
          "Renormalized"    : "Corrected" normalized with the ratio obtained
                              by dividing the sum of the original data by the
                              sum of the data after its zero-valued intensites
                              are supplimented with values in "Predicted".
                              (ndarray or list)                              
        }
    """
    if(self.__dims == 1) : assert(len(data)  == self.__shape)
    else                 : assert(data.shape == self.__shape)
    started = clock()
    corrected, predicted, iterations = self.RemoveNA(data)
    stopped = clock()
    rfactor = self.__sum(data)/self.__sum(self.ReplaceNegatives(data, predicted))
    renormalized = numpy.multiply(corrected,rfactor)
    self.__results = {}
    # Book keeping
    self.__results["OriginalData"  ] = data
    self.__results["Iterations"    ] = iterations
    self.__results["CPUTime"       ] = stopped - started
    # The good stuff
    self.__results["TotalNARemoved"] = self.__sum(corrected)
    self.__results["Predicted"     ] = predicted
    self.__results["Corrected"     ] = corrected 
    self.__results["Renormalized"  ] = renormalized
    return self.__results

  def FilterFormatResults(self, threshold=float("inf")):
    """
    FilterFormatResults Function
    ============================
      Description
      -----------
      This function is used to organize the results data in a more output
      oriented format.  It can also be used to find high predicted
      intensities for peaks that were not observed in the original data
      set.
      
      Arguments
      ---------
        - threshold : float
          A value that the predicted intensity of an unobserved peak must be
          greater than to be included in the results. Default is inf
          (no predictions are included in the results)
    
      Returns
      -------
      A multi-layer dictionary with the following structure:
      { 
        "TotalNARemoved"      : The sum of the corrected intensities
        "CPUTime"             : The total CPU time used to remove natural
                                abundance from the data passed to Run
        "Iterations"          : Iterations required to remove natural 
                                abundance from the data passed to Run
        "AccpetedPredictions" : The number of unobserved peaks that whose
                                predicted values were above threshold which
                                were included in the results                                
        
        The "#PeakResults" element is special, its value is a dictionary of
        of dictionaries where the first level of keys are tuples specifying
        the peak's coordinate in the data set and the second level are strings
        defining the result.  It has the following structure:
        "#PeakResults"        :                             
        {
          (Peak 0's Coordinates):  
          {
            "Predicted"    : The predicted intensity for this peak.
            "Corrected"    : The corrected intensity for this peak.
            "Renormalized" : The renormalized intensity for this peak.
            "Unobserved"   : True if this peak's predicted value was above
                             threshold but was not included in the original 
                             data. False otherwise.
          }
          .
          .
          .
          (Peak n's Cordinates): { ... }
        }        
      }
    """
    assert(self.__results != None)
    accepted = 0
    data     = self.__results["OriginalData"]
    filtered = {"TotalNARemoved"      : self.__results["TotalNARemoved"],
                "CPUTime"             : self.__results["CPUTime"       ],
                "Iterations"          : self.__results["Iterations"    ],
                "AcceptedPredictions" : 0,
                "#PeakResults"        : {}}
    for i in self.__ynacrange(self.__shape):
      observed = (data[i] != 0)
      if(observed or self.__results["Predicted"][i] > threshold):
        if(not observed): accepted += 1
        peakresult = {"Predicted"    : self.__results["Predicted"   ][i],
                      "Corrected"    : self.__results["Corrected"   ][i],
                      "Renormalized" : self.__results["Renormalized"][i],
                      "Unobserved"   : not observed}
        # Tuplize the index for 1D data to keep interface uniform
        index = (i,) if(self.__dims == 1) else i 
        filtered["#PeakResults"][index] = peakresult
    filtered["AcceptedPredictions"] = accepted
    return filtered
    
  def SubtractNA(self, data):
    """
    SubtractNA Function
    ===================
      Description
      -----------
      Performs one round of Natural Abudance subtraction on the data passed in
      
      Arguments
      ---------
        - data : numpy.ndarray or list of floats
          This argument defines the data set to which predicted values of
          Natural Abudance will be subtracted from.  It must have the same 
          shape as this object
    """
    calc      = self.__copy(data) # Copy Data
    xnacrange = self.__xnacrange  # Iteration function 
    plookup   = self.__plookup    # Binomial Coefficient Product Table
    slookup   = self.__slookup    # Binomial Coefficient Sum Product Table
    for y in self.__ynacrange(self.__shape):
      yvalue = calc[y]
      for x in xnacrange(y):
        yvalue -= calc[x] * plookup[x, y]
      yvalue /= slookup[y]
      calc[y] = yvalue if yvalue > 0 else 0 # Improves accuracy
    return calc

  def AddNA(self, data):
    """
    AddNA Function
    ==============
      Description
      -----------
      Performs one round Natural Abudance addition on the data passed in
      
      Arguments
      ---------
        - data : numpy.ndarray or list of floats
          This argument defines the data set to which predicted values of
          Natural Abudance will be added to.  It must have the same shape
          as this object
    """
    calc      = self.__copy(data) # Copy data
    xnacrange = self.__xnacrange  # Iteration function
    plookup   = self.__plookup    # Binomial Coefficient Product Table
    slookup   = self.__slookup    # Binomial Coefficient Sum Product Table
    for y in self.__ynacrange(self.__maximums, self.__minimums, -1):
      yvalue = calc[y] * slookup[y]
      for x in xnacrange(y):
        yvalue += calc[x] * plookup[x, y]
      calc[y] = yvalue
    return calc
    
  def ReplaceNegatives(self, a, b):
    """
    ReplaceNegatives Function
    =========================
      Arguments
      ---------
        - a : numpy.ndarray or list of floats
          Base matrix to be used in replacement
        - b : numpy.ndarray or list of floats
          Elements that will replace elements in a
    
      Returns
      -------
      A numpy.ndarray or list copy of 'a' where elements in 'a' that are 
      smaller than or equal to 0 will be replaced by the corresponding element
      in 'b'
    """ 
    r = self.__copy(a)
    for coordinates in self.__ynacrange(self.__shape):
      if(a[coordinates] <= 0):
        r[coordinates] = b[coordinates]
    return r
  
  def RemoveNA(self, data):
    AddNA            = self.AddNA
    SubtractNA       = self.SubtractNA
    ReplaceNegatives = self.ReplaceNegatives
    Sum              = self.__sum
    Scale            = numpy.multiply
    SumDifference    = self.__sumdiff

    targetsum  = Sum(data)
    subtracted = SubtractNA(data)
    subsum     = Sum(subtracted)
    added      = AddNA(Scale(subtracted,targetsum/subsum))
    targetsum  = Sum(ReplaceNegatives(data, added))
    addedsum   = Sum(added)
    diff       = SumDifference(data, added)
    lastdiff   = diff * 2  
    iterations = 0

    while(diff < lastdiff):
      lastdiff   = diff
      nextdata   = ReplaceNegatives(data, added)
      nextsum    = Sum(nextdata)
      subtracted = SubtractNA(Scale(nextdata,targetsum/nextsum))
      subsum     = Sum(subtracted)
      added      = AddNA(Scale(subtracted,targetsum/subsum))
      addedsum   = Sum(added)
      targetsum  = Sum(ReplaceNegatives(data, added))
      diff       = SumDifference(data, added)
      iterations += 1
    return [subtracted, added, iterations] 

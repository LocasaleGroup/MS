class MDIter(object):
  def __init__(self, dimensions, start, stop, sign):
    self._dims    = dimensions
    self._start   = start if stop else [0]*dimensions
    self._current = [s for s in start if True] if stop else [0]*dimensions
    self._stop    = stop if stop else start
    self._sign    = sign

class NACIter(MDIter):
  def __init__(self, start, stop=None, step=1):
    if(step == 0): raise ValueError("NACIter() arg 3 must not be zero")
    dimensions = len(start)
    sign = cmp(step,0)
    self._step = step
    MDIter.__init__(self, dimensions, start, stop, sign)
    self._current[0] -= step

class UltimateNACIter(NACIter):
  def __iter__(self):
    return self

  def next(self):
    # Resolving to an object takes longer than using a local reference, so
    # we do this to speed up when self.__dims starts to get large 
    last    = self._dims - 1
    axis    = 0
    step    = self._step
    current = self._current
    stop    = self._stop
    start   = self._start
    sign    = self._sign
  
    current[axis] += step
    while(current[axis]*sign >= stop[axis]*sign and axis != last):
      current[axis] = start[axis]
      axis += 1
      current[axis] += step
    
    if(current[last]*sign >= stop[axis]*sign):
      raise StopIteration

    return tuple(current)


class PenultimateNACIter(NACIter): 
  def __iter__(self):
    return self

  def next(self):
    # Resolving to an object takes longer than using a local reference, so
    # we do this to speed up when self.__dims starts to get large 
    last    = self._dims - 1
    axis    = 0
    step    = self._step
    current = self._current
    stop    = self._stop
    start   = self._start
    sign    = self._sign
  
    current[axis] += step
    while(current[axis]*sign > stop[axis]*sign and axis != last):
      current[axis] = start[axis]
      axis += 1
      current[axis] += step
    
    if(current[last] == stop[last]):
      stpcrit = True
      for axis in xrange(self._dims):
        stpcrit = stpcrit and (current[axis] == stop[axis])
      if(stpcrit):
        raise StopIteration

    return tuple(current)

class ExtendedNACIter(MDIter):
  def __init__(self, start, stop=None, step=None):
    dimensions = len(start)
    sign = cmp(step[0],0) if step else 1
    MDIter.__init__(self, dimensions, start, stop, sign)
    nstep = False
    pstep = False
    # Check step for nonsense and apply
    if(step):
      for s in step:
        if(s == 0): raise ValueError("a step value cannot equal 0")
        if(s <  0): nstep = True
        else      : pstep = True
        if(nstep and pstep): raise ValueError("all step values must have the same sign")
      self.__step = step
    # Or initialize default step
    else: 
      self.__step = tuple([1]*dimensions)
    self._current[0] -= self.__step[0]

  def __iter__(self): 
    return self
  
  def next(self):
    # Resolving to an object takes longer than using a local reference, so
    # we do this to speed up when self.__dims starts to get large 
    last   = self._dims - 1
    axis    = 0
    step    = self.__step
    current = self._current
    stop    = self._stop
    start   = self._start
    sign    = self._sign

    current[axis] += step[axis]
    while((current[axis] - stop[axis])*sign > 0 and axis != last):
      current[axis] = start[axis]
      axis += 1
      current[axis] += step[axis]
    if((current[last] - stop[last])*sign > 0):
      raise StopIteration
    return tuple(current)

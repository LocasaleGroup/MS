#!/usr/bin/python
# Standard Packages
import re
# The supported isotope labels
Supported = [
'13C' ,
'2H'  ,
'15N'
]
# Element names of the supported isotope labels  
Element = {
'13C' : "Carbon",
'2H'  : "Hydrogen",
'15N' : "Nitrogen"
}
# Element symbols
ElementSymbol = {
'13C' : "C",
'2H'  : "H",
'15N' : "N"
}
# Full name of supported isotope labels
Name = {
'13C' : "Carbon-13",
'2H'  : "Deuterium",
'15N' : "Nitrogen-15"
}
# Natural abundance of the supported isotope labels
NA = {
'13C' : 0.01109,
'2H'  : 0.00015,
'15N' : 0.00370
}
# Regular expressions for formula extraction
__fregexes = {
'13C' : re.compile(".*(C|c)(?P<max>\d+).*"),
'2H'  : re.compile(".*(H|h)(?P<max>\d+).*"),
'15N' : re.compile(".*(N|n)(?P<max>\d+).*")
}

def Maximum(isotope, formula):
  if(IsSupported(isotope)):
    match   = __fregexes[isotope].match(formula)
    maximum = 0
    if(match): maximum = int(match.group("max"))
    return maximum
  else: raise Exception("Invalid isotope label '"+str(isotope)+"'")

def IsSupported(label):
  if(Supported.count(label) == 0):
    return False
  else:
    return True
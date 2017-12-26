#!/usr/bin/python
# Standard Packages
import csv

class DataSource(object):
  def __init__(self):
    self.Rows      = []
    self.MaxRow    = -1
    self.MaxColumn = -1
  
  def AddCell(self, rown, coln, value):
    # Pad to sepcified row if required
    if(rown > self.MaxRow):
      d = rown - self.MaxRow
      for i in range(d):
        self.Rows.append([])
        self.MaxRow += 1
    # Get the row specified by rown    
    row = self.Rows[rown]
    # Pad to specified column if required
    if(coln >= len(row)):
      d = coln - len(row) + 1
      for i in range(d):
        row.append(" ")
      if(len(row) - 1 > self.MaxColumn):
        self.MaxColumn = len(row) - 1
    # Set the value at rown and coln
    # print "DEBUG : {0}:{1}".format(rown,coln)
    # print "DEBUG : {0}:{1}".format(self.MaxRow,len(row))
    self.Rows[rown][coln] = value

  def DumpSheet(self):
    for row in self.Rows:
      print " ".join(row) + "\n"

class CSVDataSource(DataSource):
  def __init__(self, path):
    DataSource.__init__(self)
    self.FilePath  = path
  
  def Write(self, filepath=None):
    if(filepath == None):
      filepath = self.FilePath
    ofile = open(filepath, "wb")
    writer = csv.writer(ofile)
    for row in self.Rows:
      writer.writerow(row)
    ofile.close()
  
  def ReadAll(self, limit=0, subscriber=None):
    ifile  = open(self.FilePath, "rb")
    reader = csv.reader(ifile)
    for row in reader:
      self.Rows.append(row)
      self.MaxRow += 1
      if(limit != 0 and self.MaxRow == limit):
        break;
      if(subscriber != None):
        subscriber.OnRowRead(self.MaxRow, row)
      if(len(row) - 1 > self.MaxColumn):
        self.MaxColumn = len(row) - 1
    ifile.close() 

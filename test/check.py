from sys import *
from math import *

class Curve:
	def read(self,file,x,y):
		line = file.readline()
		while line:
			record = line.split()
			if record[0] != "#":
				self.l.append((float(record[x-1]),float(record[y-1])))
			line = file.readline()
	def __init__(self,f=None,x=1,y=2):
		self.l = []
		if f:
			self.read(open(f,'r'),x,y)
	def value(self,x):
		l,u = 0,len(self.l)-1
		if u == l:
			return self.l[u][1]
		while u - l > 1:
			i = (u + l)/2
			if self.l[i][0] > x:
				u = i
			else:
				l = i
		x0,x1,y0,y1 = self.l[l][0],self.l[u][0],self.l[l][1],self.l[u][1]
		return y0 + (x - x0)/(x1 - x0)*(y1 - y0)
	def __sub__(self,other):
		c = Curve()
		if len(self.l) < len(other.l):
			for p in self.l:
				c.l.append((p[0], p[1] - other.value(p[0])))
		else:
			for p in other.l:
				c.l.append((p[0], self.value(p[0]) - p[1]))
		return c
	def sum(self,f=lambda x: x):
		s = 0.
		p1 = None
		for p in self.l:
			if p1:
				s += (p[0] - p1[0])*f((p[1] + p1[1])/2.)
			p1 = p
		return s
	def max(self,f=lambda x: x):
		m = None
		for p in self.l:
			v = f(p[1])
			if not m or v > m:
				m = v
		return m
	def min(self,f=lambda x: x):
		m = None
		for p in self.l:
			v = f(p[1])
			if not m or v < m:
				m = v
		return m
	def mean(self):
		return self.sum()/self.sum(lambda x: 1.)
	def norm1(self):
		return self.sum(lambda x: abs(x))/self.sum(lambda x: 1.)
	def norm2(self):
		return sqrt(self.sum(lambda x: x*x)/self.sum(lambda x: 1.))
	def normi(self):
		return self.max(lambda x: abs(x))
	def write(self):
		for p in self.l:
			print p[0],p[1]

#!/usr/bin/env python -i

from __future__ import print_function

infile = 'in.simple'
me = 0

from lammps import lammps
lmp = lammps()

# run infile one line at a time

lines = open(infile,'r').readlines()
for line in lines: lmp.command(line)

lmp.command("run 10")
x = lmp.gather_atoms("x",1,3)
epsilon = 0.1
x[0] += epsilon
lmp.scatter_atoms("x",1,3,x)
lmp.command("run 1");

f = lmp.extract_atom("f",3)
print ("Force on 1 atom via extract_atom: ",f[0][0])

fx = lmp.extract_variable("fx","all",1)
print ("Force on 1 atom via extract_variable:",fx[0])


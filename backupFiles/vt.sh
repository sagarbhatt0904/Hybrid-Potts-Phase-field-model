#!/bin/sh

for (( i = 4600; i <= 5000; i=i+100 )); do
	./mmsp2vtk Ti64_LabTry.$i.dat Ti64_LabTry.$i.vtk
done

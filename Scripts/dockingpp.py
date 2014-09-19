#!/bin/python

"""
Functions used for post-docking analysis
"""

from optparse import OptionParser
import prody
from sys import exit

import numpy as np
from math import sin, cos
from subprocess import check_call
import csv

def normalize_list( l ):
	m = min( l )
	r = max( l ) - m
	l = [ ( x - m )/ r for x in l ]
	return l

def initialize_res_score ( model ):
	d = {}
	for k in model.iterResidues():
		#d[ k ] = [ k.getResnum(), getCent( k ), 0 ]
		d[ k.getResnum() ] = { 'coord':getCent( k ), 'score':0 } 
	return d
		

def getInitLig(lig, rec_init_cent, lig_init_cent):
	""" rec_init_cent, lig_init_cent should be numpy arraay 
	lig: ligand structure data"""
	init_translation = rec_init_cent - lig_init_cent
	t = prody.Transformation( np.identity( 3 ), init_translation )
	initlig = lig.copy()
	t.apply( initlig )
	return initlig

def getMeanCoord(structure):
	"""Returns a mean of *structure* coordinate in 3x1 NumPy arrays.
	*structure* canbe Prody instances of Molecule, AtomGroup, Chain, or Residue."""
	return np.mean( structure.getCoords(), axis=0 )

def getCent(structure):
	"""Returns a center coordinates of *structure* coordinate in 3x1 NumPy arrays.
	*structure* canbe Prody instances of Molecule, AtomGroup, Chain, or Residue."""
	max = np.max( structure.getCoords(), axis=0 )
	min = np.min( structure.getCoords(), axis=0 )
	cent = (max - min) / 2.0 + min
	return cent 

def getTranslation( trans ):
	i=0
	f = lambda x: x >= grid_num/2.0 and x - grid_num or x 
	trans = [ f(x) for x in trans ]
     
	x = -trans[0] * grid_size
	y = -trans[1] * grid_size
	z = -trans[2] * grid_size
	return np.array( [x, y, z] )

def preprocess_single( model, chain_name ):

	hv = model.getHierView()
	for chain in hv.iterChains():
		chain.setChids( chain_name )
		
	prody.writePDB( 'tmp.pdb', model )
	model = prody.parsePDB( 'tmp.pdb' )
	check_call( ["rm",  "tmp.pdb"] )

	pos=1
	apos=1

	hv = model.getHierView()
	for chain in hv.iterChains():
		#print "chain : " , chain
		for res in chain.iterResidues():
			res.setResnums( pos )
			for atom in res:
				atom.setSerial( apos )
				apos += 1
			if len( res.getIcode() ) == 0:
				pos += 1
	hv.update()
	return model

def getRotationMatrix( theta ):
	R = np.identity(3)

        cx = cos(theta[0])
	cy = cos(theta[1])
	cz = cos(theta[2])
	sx = sin(theta[0])
	sy = sin(theta[1])
	sz = sin(theta[2])

	R[0,0] = cx*cz - sx*cy*sz
	R[0,1] = -cx*sz - sx*cy*cz
	R[0,2] =  sx*sy

	R[1,0] = sx*cz + cx*cy*sz
	R[1,1] = -sx*sz + cx*cy*cz
	R[1,2] = -cx*sy

	R[2,0] = sy*sz
	R[2,1] = sy*cz
	R[2,2] = cy
	
	return R


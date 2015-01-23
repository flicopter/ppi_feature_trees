#!/bin/python

"""
fill
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

if __name__ == '__main__':

	usage = "usage: %prog [options] docking_out_file" 
	p = OptionParser(usage=usage)
	(options, args) = p.parse_args()

	if len(args) != 2:
		p.error( "incorrect number of arguments" )

	prody.confProDy( auto_show=False )
	
	outFileLines = open( args[0], 'r' ).read().splitlines()
	grid_num = int(outFileLines[0].split('\t')[0])
	grid_size = float(outFileLines[0].split('\t')[1])
	recFile = outFileLines[2].split('\t')[0]
	rec_pdbid = recFile[0:4]
	ligFile = outFileLines[3].split('\t')[0]
	lig_pdbid = ligFile[0:4]
	rec_init_cent = np.array([ float(x) for x in outFileLines[2].strip().split('\t')[1:] if not len(x) == 0])
	lig_init_cent = np.array([ float(x) for x in outFileLines[3].strip().split('\t')[1:] if not len(x) == 0])
	
	partnerPdbFile = args[1]
	partner = prody.parsePDB( partnerPdbFile )
	partner = partner.select( 'protein' )
	partner = preprocess_single( partner, 'C' )
	rec = prody.parsePDB( recFile )
	rec = rec.select('protein')
	rec = preprocess_single( rec, 'A' )
	lig = prody.parsePDB( ligFile )
	lig = lig.select('protein')
	lig = preprocess_single( lig, 'B' )

	interface_cent = []
	l_interface_cent = []
	interface_resnums = []
	l_interface_resnums = []
	interface_residue_pairs = []
	r_int_score = {}
	l_int_score = {}

	rec_dict = initialize_res_score( rec )
	lig_dict = initialize_res_score( lig )
	for line in outFileLines[4:254]:
		decoylig = lig.copy()

		rotation = [ float(x) for x in line.split('\t')[0:3] ] 
		rotation = getRotationMatrix( rotation ) 

		translation = [ float(x) for x in line.split('\t')[3:6] ]
		translation = getTranslation( translation )
		score = float( line.split('\t')[-1] )
	
		decoylig.setCoords( decoylig.getCoords() - lig_init_cent )

		# apply rotation
		t = prody.Transformation(rotation, np.zeros(3) ) 
		t.apply(decoylig) 
		
		decoylig.setCoords( decoylig.getCoords() + translation + rec_init_cent )
        
		mobile = rec + decoylig

		decoy_contacts_rec = prody.Contacts( rec )
		decoy_contacts_rec = decoy_contacts_rec.select( 4, decoylig )
		decoy_contacts_lig = prody.Contacts( decoylig )
		decoy_contacts_lig = decoy_contacts_lig.select( 4, rec )
	
		#lig_contact_atoms = prody.matchChains( lig, decoy_contacts_lig, overlap=0.01, subset='all' )[0][0]
		interface_resnums += ( list(set(decoy_contacts_rec.getResnums())) )
		l_interface_resnums += ( list(set(decoy_contacts_lig.getResnums()) ))
		interface_residue_pairs.append ( [set(decoy_contacts_rec.getResnums()) , set(decoy_contacts_lig.getResnums()), score] )
		crec = set( decoy_contacts_rec.getResnums() )
		clig = set( decoy_contacts_lig.getResnums() )
		for r in crec:
			rec_dict[ r ][ 'score' ] += score
		for r in clig:
			lig_dict[ r ][ 'score' ] += score
		interface_cent.append( getCent( decoy_contacts_rec ) )

	import itertools

	residue_pairs = itertools.product( range( 1, rec.numResidues() + 1), range( 1, lig.numResidues() + 1 ) )
	rlifp = dict( [ [ i, [] ] for i in residue_pairs ] )
	for recr, ligr, score in interface_residue_pairs:
		pairs = itertools.product( recr, ligr )
		for i in pairs:
			rlifp[i].append( score )

	residue_pairs = itertools.product( range( 1, rec.numResidues() + 1), range( 1, lig.numResidues() + 1 ) )
	rlifp_file = open( args[0][0:17] + ".rlifp.csv", 'w' )
	for i in residue_pairs:
		rlifp_file.write( str(i[0]) + ", " +  str(i[1]) + ", " + str( sum( rlifp[i] ) ) + "\n" )
	rlifp_file.close()

	import matplotlib
	matplotlib.use( 'Agg' )
	import matplotlib.pyplot as plt
	
	from mpl_toolkits.mplot3d import Axes3D

	fig=plt.figure()
	ax = fig.gca( projection='3d')

	prody.showProtein(rec, lw=1, A='blue')
	rec_res = sorted( rec_dict.keys() )
	rec_ires = [ x for x in rec_res if rec_dict[ x ][ 'score' ] > 0 ]
	rec_ires_zeroscore = [ x for x in rec_res if rec_dict[ x ][ 'score' ] == 0 ]
	lig_res = sorted( rec_dict.keys() )

	xs=[ rec_dict[i]['coord'][0] for i in rec_ires ]
	ys=[ rec_dict[i]['coord'][1] for i in rec_ires ]
	zs=[ rec_dict[i]['coord'][2] for i in rec_ires ]
	scores = [ rec_dict[i]['score'] for i in rec_ires]
	scores = normalize_list( scores ) 

	scores_all = [rec_dict[i]['score'] for i in rec_res ]
	scores_all_n = normalize_list( scores_all ) 
		
	ax.scatter( xs, ys, zs, s=50, c=scores, cmap=plt.cm.OrRd )
	for n, txt in enumerate(rec_ires):
		ax.annotate( rec_dict[txt]['score'], (xs[n],ys[n]) )

	xs=[ rec_dict[i]['coord'][0] for i in rec_ires_zeroscore ]
	ys=[ rec_dict[i]['coord'][1] for i in rec_ires_zeroscore ]
	zs=[ rec_dict[i]['coord'][2] for i in rec_ires_zeroscore ]
	ax.scatter( xs, ys, zs, s=50, c='black', marker="x" )

	# only applicable for non-randomized benchmark data docking 
	correct_contacts_rec = prody.Contacts( rec )
    	correct_contacts_lig = prody.Contacts( partner )
	craa = list(set( correct_contacts_rec.select( 4, partner ).getResnums() ))
	craaplus = []	
	for i in craa:
		oneaa = rec.select( 'resnum ' + str(i) )
		neighbors = list( set( prody.Contacts( rec ).select( 4, oneaa ).getResnums() ) )
		
		craaplus.append( i )
		craaplus += neighbors
	craaplus = list( set( craaplus ) )
	craaplus.sort()
	craa.sort()

	ifmarker = [ 1 if x in craaplus else 0 for x in rec_res ]
	recif_filename = args[0].replace(".out", ".recint.score.csv" )
	writer = csv.writer( file( recif_filename, 'w' ) )
	for i,j,k,l in zip( rec_res, scores_all_n, scores_all, ifmarker):
		writer.writerow( [i,j,k,l] )
		
	xs=[ rec_dict[i]['coord'][0] for i in craa ]
	ys=[ rec_dict[i]['coord'][1] for i in craa ]
	zs=[ rec_dict[i]['coord'][2] for i in craa ]

    	ax.scatter( xs, ys, zs, s=50, color='black', marker='*' )
	xs=[ rec_dict[i]['coord'][0] for i in craaplus if not i in craa]
	ys=[ rec_dict[i]['coord'][1] for i in craaplus if not i in craa]
	zs=[ rec_dict[i]['coord'][2] for i in craaplus if not i in craa]
    	ax.scatter( xs, ys, zs, s=50, color='green', marker='*' )
    	cl = getCent( correct_contacts_lig.select( 4, rec ) )
	claa = list(set( correct_contacts_lig.select( 4, rec).getResnums() ))

	plt.savefig( args[0].replace(".out",".recint.score.png") , format='png' ) 


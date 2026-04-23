#!/usr/bin/python
# description: a program for computing iterative tutte layouts of planar graphs
#              see also http://arxiv.org/abs/1708.06449
#              (c) Manfred Scheucher 2017-2023


from scipy.spatial import ConvexHull
from sys import argv
from ast import literal_eval
import xml.etree.ElementTree as ET

from scipy import optimize

def graph_2_ipe(G2,filepath,colormap):
	points = G2.get_pos()

	ipestyle = 'ipestyle.txt'
	g = open(filepath,'w')
	g.write("""<?xml version="1.0"?>
		<!DOCTYPE ipe SYSTEM "ipe.dtd">
		<ipe version="70005" creator="Ipe 7.1.4">
		<info created="D:20150825115823" modified="D:20150825115852"/>
		""")
	with open(ipestyle) as f:
		for l in f.readlines():
			g.write("\t\t"+l)
	g.write("""<page>
		<layer name="alpha"/>
		<layer name="beta"/>
		<view layers="alpha beta" active="alpha"/>\n""")
	
	# normalize
	x0 = min(x for (x,y) in points.values())
	y0 = min(y for (x,y) in points.values())
	x1 = max(1,max(x for (x,y) in points.values())-x0,1)
	y1 = max(1,max(y for (x,y) in points.values())-y0,1)
	maxval = max(x1,y1)
	
	#scale 
	M = 392
	points = {i:(100+float((points[i][0]-x0)*M)/maxval,100+float((points[i][1]-y0)*M)/maxval) for i in points}

	# write edges	
	edges = G2.edges()
	drawn = set()
	for i,j,color in edges:
		pen_text    = ' pen="heavier"'
		dashed_text = '' # ' dash="dashed"'
		layer_text  = ' layer="alpha"' if color != None else ' layer="beta"'

		x,y = points[i]
		u,v = points[j]

		#color = "black"
		col = "black" if color in colormap else colormap[color]
		g.write('<path stroke="'+col+'"'+pen_text+dashed_text+layer_text+'>\n')
		g.write(str(x)+' '+str(y)+' m\n')
		g.write(str(u)+' '+str(v)+' l\n')
		g.write('</path>\n')
		
	# write points
	for (x,y) in points.values():
		g.write('<use layer="alpha" name="mark/disk(sx)" pos="'+str(x)+' '+str(y)+'" size="normal" stroke="black"/>\n')
	
	g.write("""</page>\n</ipe>""")
	g.close()
	print ("finished ",filepath)


def tutte_layout(G,outer_face,weights):
	V = G.vertices()
	pos = dict()
	l = len(outer_face)

	a0 = pi/l+pi/2
	for i in range(l):
		ai = a0+pi*2*i/l
		pos[outer_face[i]] = (cos(ai),sin(ai))
	
	n = len(V)
	M = zero_matrix(RR,n,n)
	b = zero_matrix(RR,n,2)

	for i in range(n):
		v = V[i]
		if v in pos:
			M[i,i] = 1
			b[i,0] = pos[v][0]
			b[i,1] = pos[v][1]
		else:
			nv = G.neighbors(v)
			s = 0
			for u in nv:
				j = V.index(u)
				wu = weights[u,v]
				s += wu
				M[i,j] = -wu
			M[i,i] = s

	sol = M.pseudoinverse()*b
	return {V[i]:sol[i] for i in range(n)}


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("filename",type=str,help="input file")
parser.add_argument("informat",choices=['list','s6','g6'],help="input file is plain text. each line encodes a graph. lines either encode edge list (python list), sparse6 or graph6 format")
parser.add_argument("outformat",choices=['png','pdf','ipe'],help="output format is either png, pdf or ipe/xml")

args = parser.parse_args()
print("args",args)


filename = args.filename
informat = args.informat
outformat = args.outformat


#G2,outer_vertices = read_graph_from_ipe(filename)
ct = 0
drawall = 0

for l in open(filename):
	ct += 1
	print("graph#",ct,":",l.replace("\n",""))
	if args.informat == 'list':
		G2 = Graph(literal_eval(l))
	elif args.informat in ['s6','g6']:
		G2 = Graph(l)
	else:
		exit("format not implemented")

	G2.set_pos(G2.layout_planar())


	drawn_graphs = set()

	candidates = []
	maxsym = 0

	for outer_face in G2.faces():
		outer_face = [e[0] for e in outer_face]

		this_graph = Graph(G2)
		for v in outer_face: 
			this_graph.add_edge((-1-v,v))

		sym = this_graph.automorphism_group().order()
		if not drawall:
			if sym < maxsym: continue
			if sym > maxsym:
				maxsym = sym
				candidates = []

		gstr = this_graph.canonical_label().sparse6_string()
		if gstr in drawn_graphs: continue
		drawn_graphs.add(gstr)

		candidates.append(outer_face)

	print ("maxsym:",maxsym)

	outer_face = candidates[0]
	for C in candidates:
		if len(C) >= 20: outer_face=C


	if outer_face == None:
		print ("no outer face given! please mark outer vertices by squares!")
		exit()

	weights = dict()
	for u,v in G2.edges(labels=None):
		weights[u,v] = weights[v,u] = 0.000001

	F = G2.faces() # 2
	F = [tuple(e[0] for e in f) for f in F]

	
	weights = dict()
	for (u,v) in G2.edges(labels=None):
		weights[u,v] = weights[v,u] = 1
		
	G2.set_pos(tutte_layout(G2,outer_face,weights))

	n = len(G2)
	pos0 = G2.get_pos()


	colormap = {c:c for c in G2.edge_labels()}

	plotfile = argv[1]+"_"+str(ct)+".tutte."+outformat

	if outformat in ['pdf','png']:
		G2.plot(
			#edge_colors=G2._color_by_label(colormap),
			#edge_thickness=2,
			#vertex_size=0,
			#vertex_labels=None
			).save(plotfile)

	if outformat == 'ipe':
		graph_2_ipe(G2,plotfile,colormap)	

	print("wrote visualization to",plotfile)
	print()

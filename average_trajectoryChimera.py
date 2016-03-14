# find the trajectory
from Movie.gui import MovieDialog
from chimera.extension import manager
movies = [inst for inst in manager.instances if isinstance(inst, MovieDialog)]
if len(movies) != 1:
	raise AssertionError("not exactly one MD Movie")
movie = movies[0]

if hasattr(movie, 'findCoordSet'):
	findCoordSet = movie.findCoordSet
	LoadFrame = movie._LoadFrame
else:
	findCoordSet = movie.model.findCoordSet
	LoadFrame = movie.model.LoadFrame

# load all frames
for frameNum in range(movie.startFrame, movie.endFrame+1, 1):
	if not findCoordSet(frameNum):
		movie.status("loading frame %d" % frameNum)
		LoadFrame(frameNum, makeCurrent=False)

# fetch coords into arrays
movie.status("fetching coords")
atoms = movie.model._mol.atoms
coords = []
coordSets = []
from chimera.match import _coordArray
for fn in range(movie.startFrame, movie.endFrame+1, 1):
	cs = findCoordSet(fn)
	coordSets.append(cs)
	coords.append(_coordArray(atoms, coordSet=cs))

# create an averaged equivalent
movie.status("computing averages")
window = 2 # number of adjacent frames to average
averaged = []
import numpy
for i in range(len(coords)):
	weightTot = 0
	avg = numpy.zeros(coords[i].shape, type(coords[i][0]))
	for j in range(i-window, i+window+1):
		if j < 0 or j >= len(coords):
			continue
		weight = window + 1 - abs(i-j)
		weightTot += weight
		avg += weight * coords[j]
	avg /= weightTot
	averaged.append(avg)

# set coords to averaged
movie.status("setting coords")
from chimera import Point
for cs, avg in zip(coordSets, averaged):
	for a, data in zip(atoms, avg):
		a.setCoord(Point(*data), cs)
	
movie.status("trajectory averaged")

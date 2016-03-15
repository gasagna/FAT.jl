# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
import Base: call
export PolyHedralCell, volume, centre, facesIDs, HexaCell

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Polyhedral cell definition ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
immutable PolyHedralCell{T<:Real, M} 
	facesIDs::NTuple{M, UInt32}  # Tuple of face IDs. The parameter M defines the number
								 # of faces of the cell, i.e. 4 for tetrahedra, 6 for hexahedra, ...
	centre::Point3D{T}           # The centre of mass of the cell. The parameter T defines
								 # the precision used for all calculations
	volume::T                    # Volume of the cell.
	# Note that we need to provide the centre and volume ourselves, since we only provide the 
	# point IDs and not the points themselves. These calculations are done by appropriate methods
end

# outer constructor
PolyHedralCell{T, M}(faces::NTuple{M, UInt32}, centre::Point3D{T}, volume::T) = 
	PolyHedralCell{T, M}(faces, centre, volume)

@inline volume(c::PolyHedralCell) = c.volume
@inline centre(c::PolyHedralCell) = c.centre
@inline facesIDs(c::PolyHedralCell) = c.facesIDs

""" Compute centroid and volume of any polyhedral cell from its faces.

	Notes
	-----
	No verification is made that faces are congruent, i.e. they
	are defined by the same set of points. We trust OpenFOAM here.

	References
	----------
	src/OpenFOAM/meshes/primitiveMesh/primitiveMeshCellCentresAndVols.C
	see page 159 in: The Finite Volume Method in Computational 
					 Fluid Dynamics: An Advanced Introduction with 
					 OpenFOAM® and Matlab
"""
function _centreAndVolume{T, M}(faces::NTuple{M, PolygonalFace{T}}) 
	# estimated cell centre, based on average face centre
	cEst = Point3D(zero(T), zero(T), zero(T))
	for face in faces
		cEst += centre(face)
	end
	cEst /= M

	# now compute volume weighted average of pyramid centroids
	vol = zero(T)
	ctr = Point3D(zero(T), zero(T), zero(T))
	for face in faces
		pyrvol = area(face)*distance(centre(face), cEst)
		pyrctr = centre(face)*0.75 + cEst*0.25
		vol += pyrvol
		ctr += pyrctr*pyrvol
	end
	ctr /= vol
	ctr, vol/3.0
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
for (name, M) in zip([:TetraCell, :PentaCell, :HexaCell], [4, 5, 6])
	# create type 
	@eval typealias $name{T} PolyHedralCell{T, M}

	# overload call function to allow using a constructor without type parameter
	@eval call{T}(q::Type{$name}, 
				  faces::NTuple{M, UInt32}, 
				  centre::Point3D{T}, 
				  volume::T) = PolyHedralCell{T, M}(faces, centre, volume)
end
# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
import Base: call

export PolyHedralCell,
	   volume,
	   centre,
	   facesIDs,
	   HexaCell

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Polyhedral cell definition ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
	A polyhedral cell

    Parameters
    ----------
    T : datatype at which computations will be done
    M : number of faces defining the cell

    Fields
    ------
	facesIDs : tuple of face IDs
	centre   : centre of mass of the cell
	volume   : volume of the cell.
"""
immutable PolyHedralCell{T<:Real, M} 
	facesIDs::NTuple{M, UInt32} 
	centre::Point{T}            
	volume::T                   
end

#= 
	Outer constructor. Note that we need to provide the centre 
	and volume ourselves, since we only provide the face IDs and 
	not the faces themselves.
=#
PolyHedralCell{T, M}(faces::NTuple{M, UInt32}, 
					 centre::Point{T}, volume::T) = 
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
	cEst = Point(zero(T), zero(T), zero(T))
	for face in faces
		cEst += centre(face)
	end
	cEst /= M

	# now compute volume weighted average of pyramid centroids
	vol = zero(T)
	ctr = Point(zero(T), zero(T), zero(T))
	for face in faces
		pyrvol = area(face)*distance(centre(face), cEst)
		pyrctr = centre(face)*3/4 + cEst*1/4
		vol += pyrvol
		ctr += pyrctr*pyrvol
	end
	ctr /= vol
	ctr, vol/3.0
end
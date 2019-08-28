# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
import Statistics

""" Compute centroid and volume of any polyhedral cell based on the 
    areas and centres of its faces.

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
function _centreAndVolume(areas::AbstractVector{T},
                        centres::AbstractVector{<:Point},
                              N::Integer) where {T<:Real}
    # estimated cell centre, based on average face centre
    cEst = Point(zero(T), zero(T), zero(T))
    for i = 1:N
        cEst += centres[i]
    end
    cEst /= N

    # now compute volume weighted average of pyramid centroids
    vol = zero(T)
    ctr = Point(zero(T), zero(T), zero(T))
    for i = 1:N
        pyrvol = areas[i]*distance(centres[i], cEst)
        pyrctr = centres[i]*3/4 + cEst*1/4
        vol += pyrvol
        ctr += pyrctr*pyrvol
    end

    return ctr/vol, vol/3.0
end
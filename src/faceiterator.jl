# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

""" 
        faceiterator(m, :region)

    Return appropriate ranges to iterate over a given mesh region.

    This function returns two range objects that are used to perform
    operations over a specific region of the mesh, e.g. a boundary face,
    or the internal mesh faces. 

        - The first range is used to point to vectors such as fowners, 
          fsvecs, ..., i.e. it corresponds to the face ID. Example:
          the first element of the first range of `faceiterator(m, :inlet)` 
          is the ID of the first face that stay on the `inlet` patch.
        - The second range is used to point to the boundaryField vector, 
          used to contain the solution on the domain boundary in the
          ScalarField object in field.jl. 

    An example of usage of this is in the function `der!` in src/field.jl.
""" 
function faceiterator(m::Mesh, region::Symbol)
    # Define the faceIDs range, which points to vectors such as 
    # fowners, fsvecs, ...
    if haskey(m.patches, region)
        patch = m.patches[region]
        faceIDs = firstfaceID(patch):lastfaceID(patch)
    elseif region == :internal 
        faceIDs = 1:ninternalfaces(m) 
    elseif region == :boundary 
        faceIDs = (ninternalfaces(m)+1):nfaces(m)
    else
        error("`region` not understood")
    end
    # Define the fieldIDs, which is used to point to the boundaryField
    # vector, used to contain the solution on the domain boundary in the
    # ScalarField object in field.jl. 
    # If we are not on the :internal region, we need to subtract
    # the number of internal faces, as the boundary faces come after 
    # the internal faces in vectors such as fsvecs. 
    fieldIDs = region == :internal ? facesIDs : facesIDs - ninternalfaces(m)
    return zip(faceIDs, fieldIDs)
end
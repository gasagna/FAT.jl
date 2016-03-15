# Generic iterator over faces
immutable FaceIterator{I}
	facesVec::I
	first::Int # index of the first
	last::Int  # index of the last
	function FaceIterator(facesVec::I, first::Int, last::Int)
		last >= first && first >=1 && last <= length(facesVec) || error("wrong face iterator specification")
		new(facesVec, first, last)
	end
end
# efficient iterators over the internal and boundary faces
internalfaces(m::HexaMesh) = FaceIterator{typeof(m.faces)}(m.faces, 1, m.ninternalfaces)
boundaryfaces(m::HexaMesh) = FaceIterator{typeof(m.faces)}(m.faces, m.ninternalfaces+1, nfaces(m))


# iterate over the faces of a particular patch
faces(m::HexaMesh, p::Symbol) = 
	FaceIterator{typeof(m.faces)}(m.faces, 
								  firstfaceID(m.patches[p]), 
								  lastfaceID(m.patches[p]))

start(fi::FaceIterator) = fi.first
next(fi::FaceIterator, i::Int) = fi.facesVec[i], i+1
done(fi::FaceIterator, i::Int) = i == fi.last + 1
length(fi::FaceIterator) = fi.last - fi.first + 1
eltype{I}(::Type{FaceIterator{I}}) = eltype(I)

# Enumerate over the faces of the FaceIterator.
# This will return the appropriate face ID and the face, as using the 
# standard enumerate would start from 1. Do not get confused.
# the following lines are copied from iterator.jl in Julia Base
immutable EnumerateFaces{I}
    faceitr::I
end
enumerateFaces(faceitr::FaceIterator) = EnumerateFaces(faceitr)

length(e::EnumerateFaces) = length(e.faceitr)
start(e::EnumerateFaces) = (e.faceitr.first, start(e.faceitr))
function next(e::EnumerateFaces, state)
    n = next(e.faceitr, state[2])
    (state[1], n[1]), (state[1]+1, n[2])
end
done(e::EnumerateFaces, state) = done(e.faceitr, state[2])
eltype{I}(::Type{EnumerateFaces{I}}) = Tuple{Int, eltype(I)}
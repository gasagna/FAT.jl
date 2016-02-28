import Base: isempty

# A patch is just an index of the faces that form that patch.
# The OpenFoam notation of the `boundary` files is used here.
immutable Patch
	name::Symbol
	isempty::Bool
	nfaces::UInt32
	startface::UInt32
end

show(io::IO, p::Patch) = 
	print(io, "Patch '$(patchname(p))': first face ID $(firstfaceID(p)), last face ID $(lastfaceID(p)), #faces $(nfaces(p))\n")

" Is is an empty patch "
isempty(p::Patch) = p.isempty

" Patch name "
patchname(p::Patch) = p.name

" Number of faces in a Patch "
@inline nfaces(p::Patch) = Int(p.nfaces)

" ID of first face - This is 1-based "
@inline firstfaceID(p::Patch) = Int(p.startface)

" ID of last face - This is 1-based "
@inline lastfaceID(p::Patch) = Int(p.startface + p.nfaces - 1) 


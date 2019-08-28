# ------------------------------------------------------------------- #
# Copyright 2015-2019, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module Constants

# Openfoam types for boundary conditions and patches
const BC_FIXEDVALUE = 0x00000001   # The value at the face is given and can be loaded from file
const BC_EMPTY = 0x00000002        # The face value is never used

end
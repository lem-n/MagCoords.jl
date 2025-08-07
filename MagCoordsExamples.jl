
# If needed...
# import Pkg; Pkg.add("Interpolations")

using MagCoords
using Dates

# The really quick, cut-and-paste examples:
#-----------------------------------------

# calculate Poker Flat ISR MLAT/MLON (doesn't change very much with time, only secular variation)
# functions cartesian and spherical convert *to* those coordinates. All rotation matrices work on cartesian vectors.
# Lat/Lon numbers for PFISR are from Google Earth
mag_coords = spherical(MAG_from_GEO(DateTime(2013, 3, 17, 0, 0, 0, 0))*cartesian([1.0, 65.13, -147.471]))
dipole_mapped_mlat(mag_coords[2], 0.0, 100)   # dipole-mapped latitude of Poker Flat at at RCM grid altitude
dipole_mapped_mlat(mag_coords[2], 0.0, 850)   # dipole-mapped latitude of Poker Flat at DMSP altitude

# magnetic coordinates of a GEO spacecraft at 0 deg longitude
mag_coords = spherical(MAG_from_GEO(DateTime(2018, 1, 1, 0, 0, 0, 0))*cartesian([6.6, 0.0, 0.0]))
L = L_shell(cartesian(mag_coords))

# Somewhat more depth/explanation:
#---------------------------------

# if we want to transform coordinates at noon UT during the March 17 2013 storm ...
at_time = DateTime(2013, 3, 17, 12, 0, 0, 0)

# 6 coordinate systems implemented: GEO, MAG, GSM, GSE, SM, GEI
# with function names like GEO_from_GEI, GSM_from_GSE, etc.
#
# all matrix operations work in cartesian coordinates (x,y,z) (arbitrary units)
# function "spherical" converts from (x,y,z) to (r, θ, ϕ)
#    r is radius, θ is latitude (-90 to 90), ϕ is azimuth (0-360)
#
# function "cartesian" converts from (r, θ, ϕ) to (x, y, z)
#    all cartesian units are arbitrary (meters, RE, lightyears, etc)
#        conversion functions "degrees" and "radians" are available

# rotation matrix to go from SM to GEO coordinates
rotation_matrix = GEO_from_SM(at_time)

# for the same time, calculate geographic coordinates of north magnetic pole x,y,z = (0, 0, 1)
geo_coords_cartesian = GEO_from_MAG(at_time) * [0, 0, 1]

# if MAG location is given in spherical (r, θ, ϕ) ...
geo_coords_cartesian = GEO_from_MAG(at_time) * cartesian([1, 90, 0])
geo_coords_cartesian = GEO_from_MAG(at_time) * degrees(cartesian([1, π/2, 0]))

# convert result to spherical coordinates (r, θ, ϕ)
geo_coords_spherical =         spherical(geo_coords_cartesian)
geo_coords_spherical = radians(spherical(geo_coords_cartesian))

#---------------
# Array of times (take advantage of the "dot" operator for functions)
# --------------

# set up times array, every 10 minutes on 2015-03-17
times = DateTime(2015, 3, 17, 0, 0, 0, 0):Dates.Minute(10):DateTime(2015, 3, 20, 0, 0, 0, 0)

# the compact (but not so easy to read) way to calculate SM coordinates (r, θ, MLT) of Poker Flat at each time in times
PokerSM = [MLTspherical(mat*cartesian([1, 65.1367, -147.4472])) for mat in SM_from_GEO.(times)]

times = DateTime(2019, 4, 5, 19, 0, 0, 0):Dates.Minute(10):DateTime(2019, 4, 6, 0, 0, 0, 0)
AndøyaSM = [MLTspherical(mat*cartesian([1, 69.2944, 16.0198])) for mat in SM_from_GEO.(times)]

# The easier-to-read way to calculate SM coordinates of Poker Flat at each time in times
rotation_matrices = SM_from_GEO.(times)  # calculate the rotation matrix (SM to GEO) for each time in times
PokerGEOcart = cartesian([1, 65.1367, -147.4472])                           # returns [x, y, z]_GEO (arbitrary units)
PokerSMcart = [rot_mat*PokerGEOcart for rot_mat in rotation_matrices]      # returns [x, y, z]_SM  (arbitrary units)
PokerSM = MLTspherical.(PokerSMcart)                                       # returns [r, θ, MLT]_SM (units [arbitrary, degrees, hours])


times = DateTime(2013, 3, 17, 0, 0, 0, 0):Dates.Minute(10):DateTime(2013, 3, 20, 0, 0, 0, 0)
PokerSM = [MLTspherical(mat*cartesian([1, 65.1367, -147.4472])) for mat in SM_from_GEO.(times)]

for t in eachindex(times)
    (r, θ, MLT) = PokerSM[t]
    print(times[t], "  ", θ, ' ', MLT, "\n")
end

times = DateTime(2015, 3, 17, 0, 0, 0, 0):Dates.Minute(10):DateTime(2015, 3, 20, 0, 0, 0, 0)
PokerSM = [MLTspherical(mat*cartesian([1, 65.1367, -147.4472])) for mat in SM_from_GEO.(times)]

mlt2aloct(mlt)   = ((mlt + 12) % 24)*π/12  # RCM-E local time coordinate
for t in eachindex(times)
    (r, θ, MLT) = PokerSM[t]
    print(times[t], "  ", θ, ' ', MLT, (90 - θ)*π/180, 3 + 150*mlt2aloct(MLT)/2π, "\n")

end


times = DateTime(2013, 3, 17, 0, 0, 0, 0):Dates.Minute(10):DateTime(2013, 3, 20, 0, 0, 0, 0)
MillstoneSM = [MLTspherical(mat*cartesian([1, 42.82, -71.5])) for mat in SM_from_GEO.(times)]

for t in eachindex(times)
    (r, θ, MLT) = PokerSM[t]
    print(times[t], "  ", θ, ' ', MLT, "\n")
end

times = DateTime(2015, 3, 17, 0, 0, 0, 0):Dates.Minute(10):DateTime(2015, 3, 20, 0, 0, 0, 0)
MillstoneSM = [MLTspherical(mat*cartesian([1, 42.82, -71.5])) for mat in SM_from_GEO.(times)]
MillstoneMAG = [MLTspherical(mat*cartesian([1, 42.82, -71.5])) for mat in MAG_from_GEO.(times)]

mlt2aloct(mlt)   = ((mlt + 12) % 24)*π/12  # RCM-E local time coordinate
for t in eachindex(times)
    (r, θ, MLT) = MillstoneSM[t]
    print(times[t], "  ", θ, ' ', MLT, (90 - θ)*π/180, 3 + 150*mlt2aloct(MLT)/2π, "\n")

end





# You might need to...
# import Pkg; Pkg.add("Interpolations")


module MagCoords
using Dates
using Interpolations
using LinearAlgebra

export GEO_from_GEI, GSE_from_GEI, GSE_from_GEO, GSM_from_GSE, GSM_from_GEI, GSM_from_GEO,
        SM_from_GSM,  SM_from_GSE,  SM_from_GEI,  SM_from_GEO,  SM_from_MAG, MAG_from_GEO,
       MAG_from_GEI, MAG_from_GSE, MAG_from_GSM, GEI_from_GEO, GEI_from_GSE, GEI_from_GSM,
       GEI_from_SM,  GEI_from_MAG, GEO_from_GSE, GEO_from_GSM, GEO_from_SM,  GEO_from_MAG,
       GSE_from_GSM, GSE_from_SM,  GSE_from_MAG, GSM_from_SM,  GSM_from_MAG, MAG_from_SM,
       cartesian, spherical, degrees, radians, LTspherical, MLTspherical,
       cartesian_from_LTspherical, cartesian_from_MLTspherical, L_shell, dipole_mapped_mlat, RE_km

# IGRF coefficients; earlier years are available
# 2030 is made up for extrapolation purposes, using the "SV" (secular variation) numbers from the latest IGRF coefficient.
# If you're here because you're trying to run a date in 2030 or beyond, and it didn't work, it's easy to update the coefficients below
# https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt
# https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field
years      =  [  1970.0,   1975.0,   1980.0,   1985.0,   1990.0,   1995.0,   2000.0,    2005.0,    2010.0,   2015.0,      2020.0,   2025.0,             2030.0,]
g1_0_array =  [-30220.0, -30100.0, -29992.0, -29873.0, -29775.0, -29692.0, -29619.4, -29554.63, -29496.57, -29441.46,  -29403.41, -29350.0,  -29350.0 + 5*12.6,]
g1_1_array =  [ -2068.0,  -2013.0,  -1956.0,  -1905.0,  -1848.0,  -1784.0,  -1728.2,  -1669.05,  -1586.42,  -1501.77,   -1451.37,  -1410.3,   -1410.3 + 5*10.0,]
h1_1_array =  [  5737.0,   5675.0,   5604.0,   5500.0,   5406.0,   5306.0,   5186.1,   5077.99,   4944.26,   4795.99,    4653.35,   4545.5,    4545.5 + 5*(-21.5),]

# interpolated IGRF parameters
g1_0 = scale(interpolate(g1_0_array, BSpline(Quadratic(Natural(OnGrid())))), years[1]:5:years[end])
g1_1 = scale(interpolate(g1_1_array, BSpline(Quadratic(Natural(OnGrid())))), years[1]:5:years[end])
h1_1 = scale(interpolate(h1_1_array, BSpline(Quadratic(Natural(OnGrid())))), years[1]:5:years[end])

# magnetic dipole properties as calculated from interpolated IGRF coefficients
y(T₀) = 2000 + 0.5/365.25 + 100*T₀             # decimal year to be used for IGRF coefficient interpolation
earth_magnetic_moment(T₀) = norm([ h1_1(y(T₀)), g1_1(y(T₀)), g1_0(y(T₀)) ])  # Fränz and Harper Equation 21
magnetic_Npole_glon(T₀) = (180/π)*atan(-h1_1(y(T₀)), -g1_1(y(T₀)))
magnetic_Npole_glat(T₀) = 90 - acosd(-g1_0(y(T₀))/earth_magnetic_moment(T₀))      # also sometimes calculated from arctan


# Earth Magnetosphere coordinate transformation matrices
# most everything here is done in degrees because that's what astronomical measurements are generally made in (deg/min/arcsec)
ΔT = 67.6439/3600/24/36525   # this is actually the value for 2015, but it's a miniscule correction anyway

# Time (in Julian centuries) since the J2000 epoch (2012-1-1, 12:00:00) is used as the time variable for most functions
J2000century(datetime) = (Dates.datetime2julian(datetime) - 2451545.0)/36525
d₀(T₀) = T₀*36525.0                             # days since J2000.0
r360(angle) = (angle % 360 + 360) % 360         # return a value equivalent to the argument, but between 0 and 360 degrees
fday(T₀) = d₀(T₀) - floor(d₀(T₀))               # fraction of a day (between 0.0 and 1.0)

# Greenwich Mean Siderial Time; Capitaine, 2003, equation B.2 (assuming the UT1 and TT offset is zero)
# See also IERS Conventions, Gerard Petit and Brian Luzum, 2010; Section 5.7
GMST(T₀) = r360(180 + 360*(fday(T₀) + (307.47716*ΔT + @evalpoly T₀ 24110.5493771 8640184.7945360 0.0931118 -0.0000062 0.0000013)/86400))
Λ(T₀) = r360(GMST(T₀) - fday(T₀)*360) # Mean Longitude of the Sun

# Mean Anomaly of the Sun; IERS Conventions, Gerard Petit and Brian Luzum, 2010
M(T₀) = r360(357.52910918 + (@evalpoly T₀ 0 129596581.0481 -0.5532 0.000136 -0.00001149)/3600)

# Mean Longitude of the Sun; Simon et al., 1994, Numerical Expressions for precession formulae, 5.9.3
#Λ(T₀) = r360(280.46645683 + (@evalpoly T₀ 0 12960277.1103429 -109.15809 0.07207 0.23530 0.00180 0.00020)/3600)

# Hilton et al., Celestial Mechanics and Dynamical Astronomy (2006) 94:351–367  (Table I, ϵA)
sun_mean_obliquity(T₀) = (@evalpoly T₀ 84381.406000 -46.836769 -0.0001831 0.00200340 -5.76e-7 -4.34e-8)/3600

# From Meeus: http://www.geoastro.de/elevaz/basics/meeus.htm
# Meeus calls the sun ecliptic longitude the true longitude
sun_ecliptic_longitude(T₀) = r360(Λ(T₀) + (1.914600 - 0.004817*T₀ - 0.000014*T₀*T₀)*sind(M(T₀)) +
                                       (0.019993 - 0.000101*T₀)*sind(2*M(T₀)) + 0.000290*sind(3*M(T₀)))

# single-dimension rotation matrices needed to compute coordinate transformation matrices
Rx(angle) = [1       0            0     ;
             0  cosd(angle)  sind(angle);
             0 -sind(angle)  cosd(angle)]

Ry(angle) = [cosd(angle)  0  sind(angle);
                  0       1       0     ;
            -sind(angle)  0  cosd(angle)]

Rz(angle) = [cosd(angle)  sind(angle)  0;
            -sind(angle)  cosd(angle)  0;
                  0            0       1]

# single-dimension transformation matrices from Hapgood, 1992 (Equations 2-11); Also see Fränz and Harper 2002, section 3.3.3
T1(T₀) = Rz(GMST(T₀))
T2(T₀) = Rz(sun_ecliptic_longitude(T₀))*Rx(sun_mean_obliquity(T₀))

Qg(T₀) = [cosd(magnetic_Npole_glat(T₀))*cosd(magnetic_Npole_glon(T₀)),
       cosd(magnetic_Npole_glat(T₀))*sind(magnetic_Npole_glon(T₀)),
       sind(magnetic_Npole_glat(T₀))]
Qe(T₀) = T2(T₀)*inv(T1(T₀))*Qg(T₀)
xₑ(T₀) = Qe(T₀)[1];    yₑ(T₀) = Qe(T₀)[2];    zₑ(T₀) = Qe(T₀)[3]

ψ(T₀) = atand(yₑ(T₀) / zₑ(T₀))
μ(T₀) = atand(xₑ(T₀)) / √(yₑ(T₀)^2 + zₑ(T₀)^2)

T3(T₀) = Rx(-ψ(T₀))
T4(T₀) = Ry(-μ(T₀))
T5(T₀) = Ry(magnetic_Npole_glat(T₀) - 90.0)*Rz(magnetic_Npole_glon(T₀))


# coordinate transformation matrices from Hapgood, 1992, Table 3
GEO_from_GEI(T₀) = T1(T₀)                 # relies only on GMST
GSE_from_GEI(T₀) = T2(T₀)                 # This is "true" GEI, not "mean"
GSE_from_GEO(T₀) = T2(T₀)*inv(T1(T₀))     # relies on all sun location parameters
GSM_from_GSE(T₀) = T3(T₀)                 # relies on all sun parameters + IGRF coefficients
GSM_from_GEI(T₀) = T3(T₀)*T2(T₀)
GSM_from_GEO(T₀) = T3(T₀)*T2(T₀)*inv(T1(T₀))
SM_from_GSM(T₀)  = T4(T₀)
SM_from_GSE(T₀)  = T4(T₀)*T3(T₀)
SM_from_GEI(T₀)  = T4(T₀)*T3(T₀)*T2(T₀)
SM_from_GEO(T₀)  = T4(T₀)*T3(T₀)*T2(T₀)*inv(T1(T₀))
SM_from_MAG(T₀)  = T4(T₀)*T3(T₀)*T2(T₀)*inv(T1(T₀))*inv(T5(T₀))
MAG_from_GEO(T₀) = T5(T₀)
MAG_from_GEI(T₀) = T5(T₀)*T1(T₀)
MAG_from_GSE(T₀) = T5(T₀)*T1(T₀)*inv(T2(T₀))
MAG_from_GSM(T₀) = T5(T₀)*T1(T₀)*inv(T2(T₀))*inv(T3(T₀))
GEI_from_GEO(T₀) = inv(GEO_from_GEI(T₀))
GEI_from_GSE(T₀) = inv(GSE_from_GEI(T₀))
GEI_from_GSM(T₀) = inv(GSM_from_GEI(T₀))
GEI_from_SM(T₀)  = inv( SM_from_GEI(T₀))
GEI_from_MAG(T₀) = inv(MAG_from_GEI(T₀))
GEO_from_GSE(T₀) = inv(GSE_from_GEO(T₀))
GEO_from_GSM(T₀) = inv(GSM_from_GEO(T₀))
GEO_from_SM(T₀)  = inv( SM_from_GEO(T₀))
GEO_from_MAG(T₀) = inv(MAG_from_GEO(T₀))
GSE_from_GSM(T₀) = inv(GSM_from_GSE(T₀))
GSE_from_SM(T₀)  = inv( SM_from_GSE(T₀))
GSE_from_MAG(T₀) = inv(MAG_from_GSE(T₀))
GSM_from_SM(T₀)  = inv(SM_from_GSM(T₀))
GSM_from_MAG(T₀) = inv(MAG_from_GSM(T₀))
MAG_from_SM(T₀)  = inv( SM_from_MAG(T₀))

# wrapper for calling the transformation matrix functions with a DateTime object
GEO_from_GEI(t::DateTime) = GEO_from_GEI(J2000century(t))
GSE_from_GEI(t::DateTime) = GSE_from_GEI(J2000century(t))
GSE_from_GEO(t::DateTime) = GSE_from_GEO(J2000century(t))
GSM_from_GSE(t::DateTime) = GSM_from_GSE(J2000century(t))
GSM_from_GEI(t::DateTime) = GSM_from_GEI(J2000century(t))
GSM_from_GEO(t::DateTime) = GSM_from_GEO(J2000century(t))
 SM_from_GSM(t::DateTime) =  SM_from_GSM(J2000century(t))
 SM_from_GSE(t::DateTime) =  SM_from_GSE(J2000century(t))
 SM_from_GEI(t::DateTime) =  SM_from_GEI(J2000century(t))
 SM_from_GEO(t::DateTime) =  SM_from_GEO(J2000century(t))
 SM_from_MAG(t::DateTime) =  SM_from_MAG(J2000century(t))
MAG_from_GEO(t::DateTime) = MAG_from_GEO(J2000century(t))
MAG_from_GEI(t::DateTime) = MAG_from_GEI(J2000century(t))
MAG_from_GSE(t::DateTime) = MAG_from_GSE(J2000century(t))
MAG_from_GSM(t::DateTime) = MAG_from_GSM(J2000century(t))
GEI_from_GEO(t::DateTime) = GEI_from_GEO(J2000century(t))
GEI_from_GSE(t::DateTime) = GEI_from_GSE(J2000century(t))
GEI_from_GSM(t::DateTime) = GEI_from_GSM(J2000century(t))
GEI_from_SM( t::DateTime) = GEI_from_SM( J2000century(t))
GEI_from_MAG(t::DateTime) = GEI_from_MAG(J2000century(t))
GEO_from_GSE(t::DateTime) = GEO_from_GSE(J2000century(t))
GEO_from_GSM(t::DateTime) = GEO_from_GSM(J2000century(t))
GEO_from_SM( t::DateTime) = GEO_from_SM( J2000century(t))
GEO_from_MAG(t::DateTime) = GEO_from_MAG(J2000century(t))
GSE_from_GSM(t::DateTime) = GSE_from_GSM(J2000century(t))
GSE_from_SM( t::DateTime) = GSE_from_SM( J2000century(t))
GSE_from_MAG(t::DateTime) = GSE_from_MAG(J2000century(t))
GSM_from_SM( t::DateTime) = GSM_from_SM( J2000century(t))
GSM_from_MAG(t::DateTime) = GSM_from_MAG(J2000century(t))
MAG_from_SM( t::DateTime) = MAG_from_SM( J2000century(t))

# Spherical<-->Cartesian transformation functions
function cartesian(vect)  # returns [x, y, z]
    r, θ, ϕ = vect   # theta is latitude, not colatitude
    [ r*sind(90-θ)*cosd(ϕ),
      r*sind(90-θ)*sind(ϕ),
      r*cosd(90-θ)]
end

function spherical(vect)  # returns r, θ, ϕ, where θ is latitude
    x, y, z = vect
    [ norm(vect),
      90 - (180.0/π)*atan(norm([x, y]), z),
           (180.0/π)*atan(y, x)]
end

# convert cartesian to spherical, but with LT/MLT instead of azimuth (r, θ, MLT)
# Need to make it do 0 to 24, maybe? (currently -12 to 12)
LTspherical(vec) = (spherical(vec) .* [1, 1, 12/180] + [0, 0, 12]) .% [1e50, 1000, 24]  # assumes input is cartesian coordinates
MLTspherical = LTspherical

# conversion from LT/MLT spherical (r, θ, MLT) to cartesian (x,y,z)
cartesian_from_LTspherical(vec) = cartesian(vec .* [1, 1, 180/12] - [0, 0, 180])  # assumes input is spherical with LT as azimuth coordinate
cartesian_from_MLTspherical = cartesian_from_LTspherical

radians(angle) = angle*π/180
degrees(angle) = angle*180/π
radians(vect::AbstractArray) = vect .* [1, π/180, π/180]  # assumes (r, θ, ϕ) vector
degrees(vect::AbstractArray) = vect .* [1, 180/π, 180/π]  # assumes (x, y, z) vector

RE_km = 6371.2  # Earth radius in km

# These need tests
dipole_mapped_mlat(mlat, from_alt_km, to_alt_km) = sign(mlat+1.0e-15)*acosd(cosd(mlat)*√((RE_km + to_alt_km)/(RE_km + from_alt_km)))

L_shell(vect) = sum(vect.^2).^1.5 / sum(vect[1:2].^2)

# Maybe we need an invarint latitude
# Conversion to/from (r,L) or (lambda, L); these can use calls to the L_shell function, for example? Or how best to compose them?
# Maybe the dipole_mapped_mlat should convert to L-shell, then change latitude or radial distance, then convert back? (is that useful? Just want to keep the number of functions to a minimum somehow)
# In other words, maybe there is a function to convert to/from SM coordinates to (r, L, MLT) or (r, L, phi)

end

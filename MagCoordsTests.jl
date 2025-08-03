using Base.Test.@test
using MagCoords
import MagCoords.fday, MagCoords.r360, MagCoords.d₀, MagCoords.GMST, MagCoords.sun_mean_obliquity, MagCoords.M
import MagCoords.Λ, MagCoords.sun_ecliptic_longitude, MagCoords.magnetic_Npole_glat, MagCoords.magnetic_Npole_glon
import MagCoords.J2000century, MagCoords.Rx, MagCoords.Ry, MagCoords.Rz


# Formulas from GEOPACK08 and Hapgood 1992 are used to test against those used in MagCoords.jl
# The expectation is that the formulas in MagCoords.jl are more accurate, so we just test to ensure they are close
# The purpose of such a test is therefore primarily as a "sanity check"
sun_mean_anomaly_H1992(T₀)     = r360(357.528    + 35999.050*floor(T₀*36525)/36525 + 0.04107*24*fday(T₀))
sun_mean_anomaly_geopack08(T₀) = r360(358.475845 + 0.9856002670*d₀(1+T₀))
sun_mean_anomaly_AA12(T₀)      = r360(357.528    + 0.98560028*d₀(T₀))    # The Astronomical Almanac, 2012, page C5

sun_mean_longitude_H1992(T₀)     = r360(280.460    + 36000.772*floor(T₀*36525)/36525 + 0.04107*24*fday(T₀))
sun_mean_longitude_geopack08(T₀) = r360(279.696678 + 0.9856473354*d₀(1+T₀))   # mean longitude of the sun (~annual cycle)
sun_mean_longitude_AA12(T₀)      = r360(280.460    + 0.9856474*d₀(T₀))     # The Astronomical Almanac, 2012, page C5
sun_mean_long_meeus98(T₀)        = r360(280.46646  + (36000.76983 + 0.0003032*T₀)*T₀)

# Hapgood1992, Equation 5
sun_ecliptic_longitude_H1992(T₀) = r360(sun_mean_longitude_H1992(T₀) + (1.915 - 0.0048*T₀)*sind(sun_mean_anomaly_H1992(T₀)) +
                                        0.020*sind(2*sun_mean_anomaly_H1992(T₀)))

sun_ecliptic_longitude_geopack08(T₀) = r360(sun_mean_longitude_geopack08(T₀)
                                       + (1.91946 - 0.004789*(1+T₀))*sind(sun_mean_anomaly_geopack08(T₀))
                                       + 0.020094*sind(2*sun_mean_anomaly_geopack08(T₀)))

sun_longitude_functions = [sun_mean_longitude_geopack08, sun_mean_longitude_H1992, sun_mean_long_meeus98, sun_mean_longitude_AA12]

# evidently this should actually be UT1 J2000 epoch-based, but very small difference (see Fränz and Harper, 2002)
GMST82(T₀)     = r360(180 + 360*(fday(T₀) + (24110.54841 + (8640184.812866 + (0.093104 - 0.0000062*T₀)*T₀)*T₀)/86400))
GMST_FH02(T₀)  = r360(280.46061837 + 360.98564736629*d₀(T₀) + 0.0003875*T₀^2 - 2.6E-8*T₀^3)
GMST_H1992(T₀) = r360(280.461 + 36000.770*floor(T₀*36525)/36525 + 15.04107*(24*fday(T₀))) # Hapgood 1992 equation 3.
GMST_geopack08(T₀) = r360(279.690983 + 0.9856473354*d₀(1+T₀) + 360*fday(1+T₀))   # Greenwich Mean Siderial Time

sun_mean_obliquity_FH2002(T₀) = (23 + 26/60 + 21.448/3600) - 0.013004167*T₀ - 0.000000164*T₀^2 + 0.000000504*T₀^3
sun_mean_obliquity_geopack08(T₀) = 23.45229 - 0.0130125*(1 + T₀)  # mean obliquity of the ecliptic (changes very slowly, with precession)

# Nutation correction to obliquity, from Fränz and Harper, 2002; here for reference as estimate of error on sun-Earth line
# Typically, magnetospheric coordinate systems are defined with the x-axis along the *mean* (not nutation corrected) Sun-Earth line
# However, it's good to confirm that the correction is tiny, in case you think a wobbly Earth might affect your result
obliquity_nutation_correction(T₀) = 0.0026*cos(125.0-0.05295*d₀(T₀)) + 0.0002*cos(200.9 + 1.97129*d₀(T₀))

# make sure rotations rotate in the correct direction, maintain magnitude, that 90 degrees is a quarter turn (units correct), etc
@test Rx(-90)*[1.0, 0.0, 0.0] == [1.0,  0.0,  0.0]
@test Ry(-90)*[1.0, 0.0, 0.0] == [0.0,  0.0,  1.0]
@test Rz(-90)*[1.0, 0.0, 0.0] == [0.0,  1.0,  0.0]
@test Rx( 90)*[1.0, 0.0, 0.0] == [1.0,  0.0,  0.0]
@test Ry( 90)*[1.0, 0.0, 0.0] == [0.0,  0.0, -1.0]
@test Rz( 90)*[1.0, 0.0, 0.0] == [0.0, -1.0,  0.0]

@test Rx(-90)*[0.0, 1.0, 0.0] == [0.0,  0.0,  1.0]
@test Ry(-90)*[0.0, 1.0, 0.0] == [0.0,  1.0,  0.0]
@test Rz(-90)*[0.0, 1.0, 0.0] == [-1.0, 0.0,  0.0]
@test Rx( 90)*[0.0, 1.0, 0.0] == [0.0,  0.0, -1.0]
@test Ry( 90)*[0.0, 1.0, 0.0] == [0.0,  1.0,  0.0]
@test Rz( 90)*[0.0, 1.0, 0.0] == [1.0,  0.0,  0.0]

@test Rx(-90)*[0.0, 0.0, 1.0] == [0.0, -1.0,  0.0]
@test Ry(-90)*[0.0, 0.0, 1.0] == [-1.0, 0.0,  0.0]
@test Rz(-90)*[0.0, 0.0, 1.0] == [0.0,  0.0,  1.0]
@test Rx( 90)*[0.0, 0.0, 1.0] == [0.0,  1.0,  0.0]
@test Ry( 90)*[0.0, 0.0, 1.0] == [1.0,  0.0,  0.0]
@test Rz( 90)*[0.0, 0.0, 1.0] == [0.0,  0.0,  1.0]

# This is just a sanity check that I've defined the inverses properly
@test GEO_from_GEI(0.0) * GEI_from_GEO(0.0) ≈ eye(3)
@test GSE_from_GEI(0.0) * GEI_from_GSE(0.0) ≈ eye(3)
@test GSE_from_GEO(0.0) * GEO_from_GSE(0.0) ≈ eye(3)
@test GSM_from_GSE(0.0) * GSE_from_GSM(0.0) ≈ eye(3)
@test GSM_from_GEI(0.0) * GEI_from_GSM(0.0) ≈ eye(3)
@test GSM_from_GEO(0.0) * GEO_from_GSM(0.0) ≈ eye(3)
@test MAG_from_GEO(0.0) * GEO_from_MAG(0.0) ≈ eye(3)
@test MAG_from_GEI(0.0) * GEI_from_MAG(0.0) ≈ eye(3)
@test MAG_from_GSE(0.0) * GSE_from_MAG(0.0) ≈ eye(3)
@test MAG_from_GSM(0.0) * GSM_from_MAG(0.0) ≈ eye(3)
@test  SM_from_GSM(0.0) * GSM_from_SM(0.0)  ≈ eye(3)
@test  SM_from_GSE(0.0) * GSE_from_SM(0.0)  ≈ eye(3)
@test  SM_from_GEI(0.0) * GEI_from_SM(0.0)  ≈ eye(3)
@test  SM_from_GEO(0.0) * GEO_from_SM(0.0)  ≈ eye(3)
@test  SM_from_MAG(0.0) * MAG_from_SM(0.0)  ≈ eye(3)
@test  SM_from_GEO(0.0) * GEO_from_GSM(0.0) * GSM_from_GEI(0.0) * GEI_from_MAG(0.0) * MAG_from_GSE(0.0) * GSE_from_SM(0.0) ≈ eye(3)
@test GEO_from_GEI(0.0) * GEI_from_GSE(0.0) * GSE_from_GSM(0.0) * GSM_from_MAG(0.0) * MAG_from_SM(0.0)  * SM_from_GEO(0.0) ≈ eye(3)

# Compare formulas used to other similar formulas to see that they are somewhat close (not expecting them to be in precise agreement)
# they are being compared over a range of years from 1900 to 2100
@test maximum(abs.(GMST.(-1.0:1/25:1.0) - 280.46)) < 1  # GMST should be almost the same every 4 years
@test maximum(abs.(GMST82.(-1.0:1/1000:1.0)                           - GMST.(-1.0:1/1000:1.0))) < 0.001
@test maximum(abs.(GMST_FH02.(-1.0:1/1000:1.0)                        - GMST.(-1.0:1/1000:1.0))) < 0.001
@test maximum(abs.(GMST_H1992.(-1.0:1/1000:1.0)                       - GMST.(-1.0:1/1000:1.0))) < 0.001
@test maximum(abs.(GMST_geopack08.(-1.0:1/1000:1.0)                   - GMST.(-1.0:1/1000:1.0))) < 0.01
@test maximum(abs.(sun_mean_obliquity_FH2002.(-1.0:0.001:1.0)         - sun_mean_obliquity.(-1.0:0.001:1.0)))      < 1.0e-4
@test maximum(abs.(sun_mean_obliquity_geopack08.(-1.0:0.001:1.0)      - sun_mean_obliquity.(-1.0:0.001:1.0)))      < 1.0e-4
@test maximum(abs.(sun_mean_anomaly_H1992.(-1.0:1/1000:1.0)           - M.(-1.0:1/1000:1.0)))       < 0.0015
@test maximum(abs.(sun_mean_anomaly_AA12.(-1.0:1/1000:1.0)            - M.(-1.0:1/1000:1.0)))       < 0.0015
@test maximum(abs.(sun_mean_anomaly_geopack08.(-1.0:1/1000:1.0)       - M.(-1.0:1/1000:1.0)))       < 0.01

for long in sun_longitude_functions
    println(maximum(abs.(long.(-1.0:1/1000:1.0)     - Λ.(-1.0:1/1000:1.0))))
end

using PyPlot
times = collect(-1.0:1/1000:1.0)
for long in sun_longitude_functions
    plot(times, (long.(times)     - Λ.(times)))
end
legend(('1', '2', '3', '4'))

@test maximum(abs.(sun_mean_longitude_geopack08.(-1.0:1/1000:1.0)     - Λ.(-1.0:1/1000:1.0)))     < 0.01
@test maximum(abs.(sun_mean_longitude_H1992.(-1.0:1/1000:1.0)         - Λ.(-1.0:1/1000:1.0)))     < 0.01
@test maximum(abs.(sun_mean_long_meeus98.(-1.0:1/1000:1.0)            - Λ.(-1.0:1/1000:1.0)))     < 0.01
@test maximum(abs.(sun_mean_longitude_AA12.(-1.0:1/1000:1.0)          - Λ.(-1.0:1/1000:1.0)))     < 0.01
@test maximum(abs.(sun_ecliptic_longitude_H1992.(-1.0:1/1000:1.0)     - sun_ecliptic_longitude.(-1.0:1/1000:1.0))) < 0.01
@test maximum(abs.(sun_ecliptic_longitude_geopack08.(-1.0:1/1000:1.0) - sun_ecliptic_longitude.(-1.0:1/1000:1.0))) < 0.01

# http://wdc.kugi.kyoto-u.ac.jp/igrf/gggm/ (online MAG-to-GEO tranformation)
# Uses geodetic coordinates: http://www.geomag.bgs.ac.uk/data_service/models_compass/coord_calc.html
# https://www.ngdc.noaa.gov/geomag-web/?model=igrf
@test abs(magnetic_Npole_glat(J2000century(DateTime(2010, 1, 1, 0, 0, 0, 0))) - 80.0) < 0.1
@test abs(magnetic_Npole_glon(J2000century(DateTime(2010, 1, 1, 0, 0, 0, 0))) + 72.2) < 0.1

test_vectors_cartesian = ([ 1.0,  0.0,  0.0],
                          [-1.0,  0.0,  0.0],
                          [ 0.0,  1.0,  0.0],
                          [ 0.0, -1.0,  0.0],
                          [ 0.0,  0.0,  1.0],
                          [ 0.0,  0.0, -1.0])

# These should be equivalent to the above, for testing purposes
test_vectors_spherical = ([1.0,   0.0,   0.0],
                          [1.0,   0.0, 180.0],
                          [1.0,   0.0,  90.0],
                          [1.0,   0.0, -90.0],
                          [1.0,  90.0,   0.0],
                          [1.0, -90.0,   0.0])

test_vectors_LTspherical = ([1.0,  0.0, 12.0],
                            [1.0,  0.0,  0.0],
                            [1.0,  0.0, 18.0],
                            [1.0,  0.0,  6.0],
                            [1.0, 90.0, 12.0],
                            [1.0,-90.0, 12.0])

[@test spherical(cart)   == sph  for (cart, sph) in zip(test_vectors_cartesian, test_vectors_spherical)]
[@test cartesian(sph)    == cart for (cart, sph) in zip(test_vectors_cartesian, test_vectors_spherical)]
[@test LTspherical(cart) == sph  for (cart, sph) in zip(test_vectors_cartesian, test_vectors_LTspherical)]
[@test cartesian_from_LTspherical(sph) == cart for (cart, sph) in zip(test_vectors_cartesian, test_vectors_LTspherical)]

@test LTspherical === MLTspherical
@test cartesian_from_LTspherical == cartesian_from_MLTspherical

@test radians.(0:1:360) ≈ 0:2π/360:2π
@test degrees.(0:2π/360:2π) ≈ 0:1:360
@test radians([1, 1, 1]) == [1, π/180, π/180]  # assumes (r, θ, ϕ) vector
@test degrees([1, 1, 1]) == [1, 180/π, 180/π]


# conversion from LT/MLT spherical (r, θ, MLT) to cartesian (x,y,z)
cartesian_from_MLTspherical = cartesian_from_LTspherical


# Table 8 from Fränz and Harper, 2002; Aug 28, 1996, 16:46:00
# Here are some transformed vectors that we will use to test our coordinate transformations
Ttime = J2000century(DateTime(1996, 8, 28, 16, 46, 0, 0))

GEOT     = [ 6.9027400, -1.6362400,  1.9166900] # technically GEOD is what we use, not GEOT, but they should be very close
GEIT     = [-5.7864335, -4.1039357,  1.9166900] # Compare GEIT and GEID, which should differ by a similar amount to GEOD and GEOT
GEID     = [-5.7864918, -4.1039136,  1.9165612]
GEIJ2000 = [-5.7840451, -4.1082375,  1.9146822]
GSED     = [ 4.0378470,  5.1182566,  3.3908764]
GSMD     = [ 4.0378470,  6.0071917,  1.2681654]
SMD      = [ 3.3601371,  6.0071917,  2.5733108]
MAGD     = [ 3.3344557,  6.0215108,  2.5732497]
HAED     = [-5.7864918, -3.0028771,  3.3908764]
HAEJ2000 = [-5.7840451, -3.0076174,  3.3908496]
HGCJ2000 = [-5.4328785,  4.1138243,  2.7493786]
HEED     = [-4.0378470, -5.1182566,  3.3908764]
HEEQD    = [-4.4132668, -5.1924440,  2.7496187]
HCD      = [-4.3379628,  5.2555187,  2.7496187]
HGRTNE   = [ 4.0360303,  5.1931904, -3.2771992]

# This is where the rubber meets the road
@test norm(GEI_from_MAG(Ttime)*MAGD - GEID) < 0.003*norm(GEID)
@test norm(GEI_from_GSM(Ttime)*GSMD - GEID) < 0.003*norm(GEID)
@test norm(GEI_from_GSE(Ttime)*GSED - GEID) < 0.003*norm(GEID)
@test norm(GEI_from_GEO(Ttime)*GEOT - GEID) < 0.003*norm(GEID)
@test norm(GEI_from_SM( Ttime)*SMD  - GEID) < 0.003*norm(GEID)

@test norm(GEO_from_MAG(Ttime)*MAGD - GEOT) < 0.003*norm(GEOT)
@test norm(GEO_from_GEI(Ttime)*GEID - GEOT) < 0.003*norm(GEOT)
@test norm(GEO_from_GSM(Ttime)*GSMD - GEOT) < 0.003*norm(GEOT)
@test norm(GEO_from_GSE(Ttime)*GSED - GEOT) < 0.003*norm(GEOT)
@test norm(GEO_from_SM( Ttime)*SMD  - GEOT) < 0.003*norm(GEOT)

@test norm(GSM_from_GSE(Ttime)*GSED - GSMD) < 0.003*norm(GSMD)
@test norm(GSM_from_MAG(Ttime)*MAGD - GSMD) < 0.003*norm(GSMD)
@test norm(GSM_from_GEI(Ttime)*GEID - GSMD) < 0.003*norm(GSMD)
@test norm(GSM_from_GEO(Ttime)*GEOT - GSMD) < 0.003*norm(GSMD)
@test norm(GSM_from_SM( Ttime)*SMD  - GSMD) < 0.003*norm(GSMD)

@test norm(GSE_from_MAG(Ttime)*MAGD - GSED) < 0.003*norm(GSED)
@test norm(GSE_from_GEI(Ttime)*GEID - GSED) < 0.003*norm(GSED)
@test norm(GSE_from_GSM(Ttime)*GSMD - GSED) < 0.003*norm(GSED)
@test norm(GSE_from_GEO(Ttime)*GEOT - GSED) < 0.003*norm(GSED)
@test norm(GSE_from_SM( Ttime)*SMD  - GSED) < 0.003*norm(GSED)

@test norm( SM_from_MAG(Ttime)*MAGD - SMD ) < 0.003*norm(SMD)
@test norm( SM_from_GEI(Ttime)*GEID - SMD ) < 0.003*norm(SMD)
@test norm( SM_from_GSM(Ttime)*GSMD - SMD ) < 0.003*norm(SMD)
@test norm( SM_from_GSE(Ttime)*GSED - SMD ) < 0.003*norm(SMD)
@test norm( SM_from_GEO(Ttime)*GEOT - SMD ) < 0.003*norm(SMD)

@test L_shell([0, 0, 5]) == Inf
@test L_shell([5, 0, 0]) == 5.0
# If we made it this far...
println("All tests passed!")
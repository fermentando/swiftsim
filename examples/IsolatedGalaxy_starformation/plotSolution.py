import matplotlib

matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py as h5
from sphviewer.tools import QuickView

# Plot parameters
params = {
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "font.size": 9,
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": True,
    "figure.figsize": (3.15, 3.15),
    "figure.subplot.left": 0.15,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.13,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.12,
    "lines.markersize": 6,
    "lines.linewidth": 2.0,
    "text.latex.unicode": True,
}
rcParams.update(params)
rc("font", **{"family": "sans-serif", "sans-serif": ["Times"]})


filename = "output_0033.hdf5"

f = h5.File(filename, "r")

# Physical constants
k_in_cgs = 1.38064852e-16
mH_in_cgs = 1.6737236e-24
year_in_cgs = 3600.0 * 24 * 365.0
Msun_in_cgs = 1.98848e33
G_in_cgs = 6.67259e-8
pc_in_cgs = 3.08567758e18
Msun_p_pc2 = Msun_in_cgs / pc_in_cgs**2

# Gemoetry info
boxsize = f["/Header"].attrs["BoxSize"]
centre = boxsize / 2.0

# Read units
unit_length_in_cgs = f["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = f["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = f["/Units"].attrs["Unit time in cgs (U_t)"]

# Calculate Gravitational constant in internal units
G = G_in_cgs * ( unit_length_in_cgs**3 / unit_mass_in_cgs / unit_time_in_cgs**2)**(-1)

# Read parameters of the model
KS_law_slope = float(f["/Parameters"].attrs["EAGLEStarFormation:SchmidtLawExponent"])
KS_law_norm = float(f["/Parameters"].attrs["EAGLEStarFormation:SchmidtLawCoeff_MSUNpYRpKPC2"])
KS_thresh_Z0 = float(f["/Parameters"].attrs["EAGLEStarFormation:MetDep_Z0"])
KS_thresh_slope = float(f["/Parameters"].attrs["EAGLEStarFormation:MetDep_SFthresh_Slope"])
KS_thresh_norm = float(f["/Parameters"].attrs["EAGLEStarFormation:thresh_norm_HpCM3"])
KS_gas_fraction = float(f["/Parameters"].attrs["EAGLEStarFormation:fg"])
KS_thresh_max_norm = float(f["/Parameters"].attrs["EAGLEStarFormation:thresh_max_norm_HpCM3"])
KS_gamma_effective = float(
    f["/Parameters"].attrs["EAGLEEntropyFloor:Jeans_gamma_effective"]
)
KS_high_den_thresh = float(f["/Parameters"].attrs["EAGLEStarFormation:SchmidtLawHighDens_thresh_HpCM3"])
KS_law_slope_high_den = float(f["/Parameters"].attrs["EAGLEStarFormation:SchmidtLawHighDensExponent"])
EAGLE_Z = float(f["/Parameters"].attrs["EAGLEChemistry:init_abundance_metal"])
EAGLEfloor_rho_norm = float(f["/Parameters"].attrs["EAGLEEntropyFloor:Jeans_density_threshold_H_p_cm3"])
EAGLEfloor_T = float(f["/Parameters"].attrs["EAGLEEntropyFloor:Jeans_temperature_norm_K"])
EAGLEfloor_cool_rho = float(f["/Parameters"].attrs["EAGLEEntropyFloor:Cool_density_threshold_H_p_cm3"])
EAGLEfloor_cool_temperature_norm_K = float(f["/Parameters"].attrs["EAGLEEntropyFloor:Cool_temperature_norm_K"])
EAGLEfloor_cool_gamma_effective = float(f["/Parameters"].attrs["EAGLEEntropyFloor:Cool_gamma_effective"])

KS_law_norm_cgs = KS_law_norm * Msun_in_cgs / ( 1e6 * pc_in_cgs**2 * year_in_cgs )
gamma = 5./3.
KS_press_norm = k_in_cgs * EAGLEfloor_T * EAGLEfloor_rho_norm

# Read gas properties
gas_pos = f["/PartType0/Coordinates"][:, :]
gas_mass = f["/PartType0/Masses"][:]
gas_rho = f["/PartType0/Density"][:]
gas_T = f["/PartType0/Temperature"][:]
gas_SFR = f["/PartType0/SFR"][:]
gas_XH = f["/PartType0/ElementAbundance"][:, 0]
gas_Z = f["/PartType0/Metallicity"][:]
gas_hsml = f["/PartType0/SmoothingLength"][:]
gas_sSFR = f["/PartType0/sSFR"][:]

# Read the Star properties
stars_pos = f["/PartType4/Coordinates"][:, :]
stars_BirthDensity = f["/PartType4/BirthDensity"][:]
stars_BirthFlag = f["/PartType4/NewStarFlag"][:]
stars_BirthTime = f["/PartType4/Birth_time"][:]
stars_XH = f["/PartType4/ElementAbundance"][:,0]

# Centre the box
gas_pos[:, 0] -= centre[0]
gas_pos[:, 1] -= centre[1]
gas_pos[:, 2] -= centre[2]

stars_pos[:,0] -= centre[0]
stars_pos[:,1] -= centre[1]
stars_pos[:,2] -= centre[2]

# Turn the mass into better units
gas_mass *= unit_mass_in_cgs / Msun_in_cgs

# Turn the SFR into better units
gas_SFR = np.maximum(gas_SFR, np.zeros(np.size(gas_SFR)))
gas_SFR /= unit_time_in_cgs / year_in_cgs
gas_SFR *= unit_mass_in_cgs / Msun_in_cgs

# Make it a Hydrogen number density
gas_nH = gas_rho * unit_mass_in_cgs / unit_length_in_cgs ** 3
gas_nH /= mH_in_cgs
gas_nH *= gas_XH

stars_BirthDensity *= unit_mass_in_cgs / unit_length_in_cgs ** 3
stars_BirthDensity /= mH_in_cgs
stars_BirthDensity *= stars_XH

# Equations of state
eos_cool_rho = np.logspace(-5, 5, 1000)
eos_cool_T = EAGLEfloor_cool_temperature_norm_K * eos_cool_rho **( EAGLEfloor_cool_gamma_effective - 1.0 )
eos_Jeans_rho = np.logspace(-1, 5, 1000)
eos_Jeans_T = EAGLEfloor_T * (eos_Jeans_rho / EAGLEfloor_rho_norm) ** (KS_gamma_effective - 1.0 ) 
# Plot the phase space diagram
figure()
subplot(111, xscale="log", yscale="log")
plot(eos_cool_rho, eos_cool_T, "k--", lw=0.6)
plot(eos_Jeans_rho, eos_Jeans_T, "k--", lw=0.6)
scatter(gas_nH, gas_T, s=0.2)
xlabel("${\\rm Density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
ylabel("${\\rm Temperature}~T~[{\\rm K}]$", labelpad=2)
xlim(1e-4, 3e3)
ylim(500.0, 2e5)
savefig("rhoT.png", dpi=200)

# Plot the phase space diagram for SF gas
figure()
subplot(111, xscale="log", yscale="log")
plot(eos_cool_rho, eos_cool_T, "k--", lw=1)
plot(eos_Jeans_rho, eos_Jeans_T, "k--", lw=1)
scatter(gas_nH[gas_SFR > 0.0], gas_T[gas_SFR > 0.0], s=0.2)
xlabel("${\\rm Density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
ylabel("${\\rm Temperature}~T~[{\\rm K}]$", labelpad=2)
xlim(1e-4, 3e3)
ylim(500.0, 2e5)
savefig("rhoT_SF.png", dpi=200)

########################################################################3

# 3D Density vs SFR
figure()
subplot(111, xscale="log", yscale="log")
scatter(gas_nH, gas_SFR, s=0.2)
plot([1, 100], 2e-5 * np.array([1, 100]) ** 0.266667, "k--", lw=1)
xlabel("${\\rm Density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
ylabel("${\\rm SFR}~[{\\rm M_\\odot~\\cdot~yr^{-1}}]$", labelpad=-7)
xlim(1e-4, 3e3)
ylim(8e-6, 2.5e-4)
savefig("rho_SFR.png", dpi=200)

########################################################################3

star_mask = (
    (stars_pos[:, 0] > -15)
    & (stars_pos[:, 0] < 15)
    & (stars_pos[:, 1] > -15)
    & (stars_pos[:, 1] < 15)
    & (stars_pos[:, 2] < 1.0)
    & (stars_pos[:, 2] > -1.0)
)

stars_BirthDensity = stars_BirthDensity[star_mask] 
stars_BirthFlag = stars_BirthFlag[star_mask]
stars_BirthTime = stars_BirthTime[star_mask]
# Histogram of the birth density
figure()
subplot(111, xscale="linear", yscale="linear")
hist(np.log10(stars_BirthDensity),density=True,bins=20,range=[-2,5])
xlabel("${\\rm Birth Density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
ylabel("Probability", labelpad=-7)
savefig("BirthDensity.png", dpi=200)

# Plot of the specific star formation rate in the galaxy
rhos = 10**np.linspace(-1,np.log10(KS_high_den_thresh),100)
rhoshigh = 10**np.linspace(np.log10(KS_high_den_thresh),5,100)

P_effective = KS_press_norm * ( rhos / EAGLEfloor_rho_norm)**(KS_gamma_effective)
P_norm_high = KS_press_norm * (KS_high_den_thresh  / EAGLEfloor_rho_norm)**(KS_gamma_effective)
sSFR = KS_law_norm_cgs * (Msun_p_pc2)**(-KS_law_slope) * (gamma/G_in_cgs * KS_gas_fraction *P_effective)**((KS_law_slope-1.)/2.)
KS_law_norm_high_den_cgs = KS_law_norm_cgs * (Msun_p_pc2)**(-KS_law_slope) * (gamma/G_in_cgs * KS_gas_fraction * P_norm_high)**((KS_law_slope-1.)/2.)
sSFR_high_den = KS_law_norm_high_den_cgs * ((rhoshigh/KS_high_den_thresh)**KS_gamma_effective)**((KS_law_slope_high_den-1)/2.)


figure()
subplot(111)
hist2d(np.log10(gas_nH), np.log10(gas_sSFR), bins=50,range=[[-1.5,5],[-.5,2.5]])
plot(np.log10(rhos),np.log10(sSFR)+np.log10(year_in_cgs)+9.,'k--',label='sSFR low density EAGLE')
plot(np.log10(rhoshigh),np.log10(sSFR_high_den)+np.log10(year_in_cgs)+9.,'k--',label='sSFR high density EAGLE')
xlabel("${\\rm Density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
ylabel("${\\rm sSFR}~[{\\rm Gyr^{-1}}]$", labelpad=-7)
savefig("density-sSFR.png", dpi=200)

########################################################################3

# Select gas in a pillow box around the galaxy
mask = (
    (gas_pos[:, 0] > -15)
    & (gas_pos[:, 0] < 15)
    & (gas_pos[:, 1] > -15)
    & (gas_pos[:, 1] < 15)
    & (gas_pos[:, 2] < 1.0)
    & (gas_pos[:, 2] > -1.0)
)
gas_pos = gas_pos[mask, :]
gas_SFR = gas_SFR[mask]
gas_nH = gas_nH[mask]
gas_rho = gas_rho[mask]
gas_T = gas_T[mask]
gas_mass = gas_mass[mask]
gas_Z = gas_Z[mask]
gas_hsml = gas_hsml[mask]


# Make a crude map of the gas
figure()
subplot(111)
scatter(gas_pos[:, 0], gas_pos[:, 1], s=0.1)
xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
xlim(-12, 12)
ylim(-12, 12)
savefig("face_on.png", dpi=200)

figure()
subplot(111)
scatter(gas_pos[:, 0], gas_pos[:, 2], s=0.1)
xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
ylabel("${\\rm Pos}~z~[{\\rm kpc}]$", labelpad=-3)
xlim(-12, 12)
ylim(-12, 12)
savefig("edge_on.png", dpi=200)

# Now a SF map
rcParams.update({"figure.figsize": (4.15, 3.15)})
figure()
subplot(111)
scatter(gas_pos[:, 0], gas_pos[:, 1], s=0.1, c=gas_SFR)
xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
colorbar()
xlim(-12, 12)
ylim(-12, 12)
savefig("SF_face_on.png", dpi=200)


########################################################################3

# Bin the data in kpc-size patches

x_edges = np.linspace(-15, 15, 31)
y_edges = np.linspace(-15, 15, 31)

map_mass, _, _, _ = stats.binned_statistic_2d(
    gas_pos[:, 0], gas_pos[:, 1], gas_mass, statistic="sum", bins=(x_edges, y_edges)
)
map_SFR, _, _, _ = stats.binned_statistic_2d(
    gas_pos[:, 0], gas_pos[:, 1], gas_SFR, statistic="sum", bins=(x_edges, y_edges)
)

qv = QuickView(
    gas_pos,
    mass=gas_mass,
    r="infinity",
    xsize=len(x_edges) - 1,
    ysize=len(y_edges) - 1,
    p=0,  # Viewing angle theta
    roll=0,  # Viewing angle phi
    plot=False,
    logscale=False,
    hsml=gas_hsml,
)

map_mass2 = qv.get_image()
extent_mass = qv.get_extent()

gas_SFR[gas_SFR <= 0] = 1e-10

qv = QuickView(
    gas_pos,
    mass=gas_SFR,
    r="infinity",
    xsize=len(x_edges) - 1,
    ysize=len(y_edges) - 1,
    p=0,  # Viewing angle theta
    roll=0,  # Viewing angle phi
    plot=False,
    logscale=False,
    hsml=gas_hsml,
)

map_SFR2 = qv.get_image()
extent_SFR = qv.get_extent()


# Mass map
figure()
subplot(111)
pcolormesh(x_edges, y_edges, np.log10(map_mass))
colorbar()
xlim(-12, 12)
ylim(-12, 12)
xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
savefig("Map_mass.png", dpi=200)

# Mass map 2
figure()
subplot(111)
imshow(np.log10(map_mass2), extent=extent_mass)
colorbar()
xlim(-12, 12)
ylim(-12, 12)
xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
savefig("Map_mass_SPHVIEWER.png", dpi=200)

# SF map
figure()
subplot(111)
pcolormesh(x_edges, y_edges, np.log10(map_SFR), vmax=-0.5, vmin=-4.5)
colorbar()
xlim(-12, 12)
ylim(-12, 12)
xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
savefig("Map_SFR.png", dpi=200)

plot_map_SFR2 = np.zeros(np.shape(map_SFR2))
plot_map_SFR2[map_SFR2 > 1e-6] = map_SFR2[map_SFR2 > 1e-6]
# SF map 2
figure()
subplot(111)
imshow(np.log10(plot_map_SFR2), extent=extent_SFR, vmax=-0.5, vmin=-4.5)
colorbar()
xlim(-12, 12)
ylim(-12, 12)
xlabel("${\\rm Pos}~x~[{\\rm kpc}]$", labelpad=0)
ylabel("${\\rm Pos}~y~[{\\rm kpc}]$", labelpad=-3)
savefig("Map_SFR_SPHVIEWER.png", dpi=200)

#########################################################################

# Give a minimum SF surface density for the plots
map_SFR[map_SFR < 1e-6] = 1e-6

# Theoretical threshold (assumes all gas has the same Z)
KS_n_thresh = KS_thresh_norm * (gas_Z[0] / KS_thresh_Z0) ** KS_thresh_slope
if np.isfinite(KS_n_thresh) == False:
    KS_n_thresh = KS_thresh_max_norm
KS_sigma_thresh = 29.0 * np.sqrt(KS_gas_fraction) * np.sqrt(KS_n_thresh)

# Theoretical KS law
KS_sigma_mass = np.logspace(-1, 3, 100)
KS_sigma_SFR = KS_law_norm * KS_sigma_mass ** KS_law_slope

# KS relation
rcParams.update({"figure.figsize": (3.15, 3.15), "figure.subplot.left": 0.18})
figure()
subplot(111, xscale="log", yscale="log")
plot(KS_sigma_mass, KS_sigma_SFR, "k--", lw=0.6)
plot([KS_sigma_thresh, KS_sigma_thresh], [1e-8, 1e8], "k--", lw=0.6)
text(
    KS_sigma_thresh * 0.95,
    2.2,
    "$\\Sigma_{\\rm c} = %.2f~{\\rm M_\\odot\\cdot pc^{-2}}$" % KS_sigma_thresh,
    va="top",
    ha="right",
    rotation=90,
    fontsize=7,
)
text(16, 10 ** (-3.5), "$n_{\\rm H,c} = %.3f~{\\rm cm^{-3}}$" % KS_n_thresh, fontsize=7)
text(
    16,
    2e-6,
    "${\\rm K\\textendash S~law}$:\n$\\Sigma_{\\rm SFR} = A \\times \\Sigma_g^n$\n$n=%.1f$\n$A=%.3f\\times10^{-4}~{\\rm M_\\odot / yr^{1} / kpc^{2}}$\n$f_{\\rm g} = %.1f$\n$\gamma_{\\rm eos} = %.3f$\n$Z=%1.4f$"
    % (
        KS_law_slope,
        KS_law_norm * 10 ** 4,
        KS_gas_fraction,
        KS_gamma_effective,
        EAGLE_Z,
    ),
    fontsize=7,
)
scatter(map_mass.flatten() / 1e6, map_SFR.flatten(), s=0.4)
xlim(0.3, 900)
ylim(3e-7, 3)
xlabel("$\\Sigma_g~[{\\rm M_\\odot\\cdot pc^{-2}}]$", labelpad=0)
ylabel(
    "$\\Sigma_{\\rm SFR}~[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$", labelpad=0
)
savefig("KS_law.png", dpi=200)
close()

plot_map_SFR2[plot_map_SFR2 <= 0] = 1e-6

rcParams.update({"figure.figsize": (3.15, 3.15), "figure.subplot.left": 0.18})
figure()
subplot(111, xscale="log", yscale="log")
plot(KS_sigma_mass, KS_sigma_SFR, "k--", lw=0.6)
plot([KS_sigma_thresh, KS_sigma_thresh], [1e-8, 1e8], "k--", lw=0.6)
text(
    KS_sigma_thresh * 0.95,
    2.2,
    "$\\Sigma_{\\rm c} = %.2f~{\\rm M_\\odot\\cdot pc^{-2}}$" % KS_sigma_thresh,
    va="top",
    ha="right",
    rotation=90,
    fontsize=7,
)
text(16, 10 ** (-3.5), "$n_{\\rm H,c} = %.3f~{\\rm cm^{-3}}$" % KS_n_thresh, fontsize=7)
text(
    16,
    2e-6,
    "${\\rm K\\textendash S~law}$:\n$\\Sigma_{\\rm SFR} = A \\times \\Sigma_g^n$\n$n=%.1f$\n$A=%.3f\\times10^{-4}~{\\rm M_\\odot / yr^{1} / kpc^{2}}$\n$f_{\\rm g} = %.1f$\n$\gamma_{\\rm eos} = %.3f$\n$Z=%1.4f$"
    % (
        KS_law_slope,
        KS_law_norm * 10 ** 4,
        KS_gas_fraction,
        KS_gamma_effective,
        EAGLE_Z,
    ),
    fontsize=7,
)
scatter(map_mass2.flatten() / 1e6, plot_map_SFR2.flatten(), s=0.4)
xlim(0.3, 900)
ylim(3e-7, 3)
xlabel("$\\Sigma_g~[{\\rm M_\\odot\\cdot pc^{-2}}]$", labelpad=0)
ylabel(
    "$\\Sigma_{\\rm SFR}~[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$", labelpad=0
)
savefig("KS_law_SPHVIEWER.png", dpi=200)

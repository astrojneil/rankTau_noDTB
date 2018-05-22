import yt
import numpy as np
from yt.visualization.base_plot_types import get_multi_plot
import matplotlib.colorbar as cb
from matplotlib.colors import LogNorm
import trident as tri

#add metallicity to dataset, constant Z = 1 Zsun
def _metallicity(field, data):
    v = data['ones']  #sets metallicity to 1 Zsun
    return data.apply_units(v, "Zsun")


fn = '../../Blob_paper2/Files/T0.3_v1000_chi300_cond/KH_hdf5_chk_0080'
ds = yt.load(fn) # load data
ds.add_field(('gas', 'metallicity'), function=_metallicity, display_name="Metallicity", units='Zsun')
tri.add_ion_fields(ds, ions=['O VI', 'Mg II', 'C IV'])

fn2 = '../../Blob_paper1/Files/T0.3_v1000_chi300/KH_hdf5_chk_0042'
ds2 = yt.load(fn2) # load data
ds2.add_field(('gas', 'metallicity'), function=_metallicity, display_name="Metallicity", units='Zsun')
tri.add_ion_fields(ds2, ions=['O VI', 'Mg II', 'C IV'])


orient = 'horizontal'

#fields = ['O_p5_number_density', 'Mg_p1_number_density', 'C_p3_number_density']
fields = ['O_p5_number_density', 'velocity_x', 'velocity_y', 'temperature']
#

# There's a lot in here:
#   From this we get a containing figure, a list-of-lists of axes into which we
#   can place plots, and some axes that we'll put colorbars.
# We feed it:
#   Number of plots on the x-axis, number of plots on the y-axis, and how we
#   want our colorbars oriented.  (This governs where they will go, too.
#   bw is the base-width in inches, but 4 is about right for most cases.
fig, axes, colorbars = get_multi_plot(3, 2, colorbar=orient, bw = 4)

slc = yt.SlicePlot(ds, 'y', fields=["Mg_p1_number_density","C_p3_number_density", "O_p5_number_density"])
proj = yt.SlicePlot(ds2, 'y', fields=["Mg_p1_number_density","C_p3_number_density", "O_p5_number_density"])

slc_frb = slc.data_source.to_frb((.3, "kpc"), 512)
proj_frb = proj.data_source.to_frb((.3, "kpc"), 512)

dens_axes = [axes[0][0], axes[1][0]]
temp_axes = [axes[0][1], axes[1][1]]
vels_axes = [axes[0][2], axes[1][2]]

for dax, tax, vax in zip(dens_axes, temp_axes, vels_axes) :

    dax.xaxis.set_visible(False)
    dax.yaxis.set_visible(True)
    dax.set_ylabel('x (kpc)')
    dax.set_yticks([15, 85.3, 170.6, 256, 341.3, 426.6, 500])
    dax.set_yticklabels(['-0.15', '-0.1', '-0.05', '0.0', '0.05', '0.1', '0.15'])
    tax.xaxis.set_visible(True)
    tax.set_xlabel('z (kpc)')
    tax.set_xticks([0, 85.3, 170.6, 256, 341.3, 426.6, 512])
    tax.set_xticklabels(['-0.15', '-0.1', '-0.05', '0.0', '0.05', '0.1', '0.15'])
    tax.yaxis.set_visible(False)
    vax.xaxis.set_visible(False)
    vax.yaxis.set_visible(False)

# Converting our Fixed Resolution Buffers to numpy arrays so that matplotlib
# can render them

slc_vel = np.array(slc_frb["O_p5_number_density"])
proj_vel = np.array(proj_frb["O_p5_number_density"])
slc_dens = np.array(slc_frb["Mg_p1_number_density"])
proj_dens = np.array(proj_frb["Mg_p1_number_density"])
slc_temp = np.array(slc_frb["C_p3_number_density"])
proj_temp = np.array(proj_frb["C_p3_number_density"])


plotmin = 1.e-17
plotmax = 1.e-5
ticks = [1.e-17, 1.e-14, 1.e-11, 1.e-8, 1.e-5]

#find the areas that are below 1e12 and set to 0.0 or Nan so they don't plot
slc_vel_2 = np.zeros((len(slc_vel), len(slc_vel[0])))
for i in range(len(slc_vel)):
    for j in range(len(slc_vel[i])):
        if slc_vel[i][j] >= plotmin:
            slc_vel_2[i][j] = slc_vel[i][j]

proj_vel_2 = np.zeros((len(proj_vel), len(proj_vel[0])))
for i in range(len(proj_vel)):
    for j in range(len(proj_vel[i])):
        if proj_vel[i][j] >= plotmin:
            proj_vel_2[i][j] = proj_vel[i][j]

slc_dens_2 = np.zeros((len(slc_dens), len(slc_dens[0])))
for i in range(len(slc_dens)):
    for j in range(len(slc_dens[i])):
        if slc_dens[i][j] >= plotmin:
            slc_dens_2[i][j] = slc_dens[i][j]

proj_dens_2 = np.zeros((len(proj_dens), len(proj_dens[0])))
for i in range(len(proj_dens)):
    for j in range(len(proj_dens[i])):
        if proj_dens[i][j] >= plotmin:
            proj_dens_2[i][j] = proj_dens[i][j]

slc_temp_2 = np.zeros((len(slc_temp), len(slc_temp[0])))
for i in range(len(slc_temp)):
    for j in range(len(slc_temp[i])):
        if slc_temp[i][j] >= plotmin:
            slc_temp_2[i][j] = slc_temp[i][j]

proj_temp_2 = np.zeros((len(proj_temp), len(proj_temp[0])))
for i in range(len(proj_temp)):
    for j in range(len(proj_temp[i])):
        if proj_temp[i][j] >= plotmin:
            proj_temp_2[i][j] = proj_temp[i][j]



plots = [dens_axes[0].imshow(slc_dens_2, origin='native', norm=LogNorm()),
         dens_axes[1].imshow(proj_dens_2, origin='native', norm=LogNorm()),
         temp_axes[0].imshow(slc_temp_2, origin='native', norm=LogNorm()),
         temp_axes[1].imshow(proj_temp_2, origin='native', norm=LogNorm()),
         vels_axes[0].imshow(slc_vel_2, origin='native', norm=LogNorm()),
         vels_axes[1].imshow(proj_vel_2, origin='native', norm=LogNorm())]

plots[0].set_clim((plotmin,plotmax))
plots[0].set_cmap("jet")
plots[1].set_clim((plotmin,plotmax))
plots[1].set_cmap("jet")
plots[2].set_clim((plotmin,plotmax))
plots[2].set_cmap("jet")
plots[3].set_clim((plotmin,plotmax))
plots[3].set_cmap("jet")
plots[4].set_clim((plotmin,plotmax))
plots[4].set_cmap("jet")
plots[5].set_clim((plotmin,plotmax))
plots[5].set_cmap("jet")

titles=[r'$\mathrm{Mg\ II\ Number\ Density}\ (\mathrm{cm^{-3}})$',
        r'$\mathrm{C\ IV\ Number\ Density}\ (\mathrm{cm^{-3}})$',
        r'$\mathrm{O\ VI\ Number\ Density}\ (\mathrm{cm^{-3}})$']

for p, cax, t in zip(plots[0:6:2], colorbars, titles):
    cbar = fig.colorbar(p, cax=cax, orientation=orient)
    cbar.set_ticks(ticks, update_ticks=True)
    cbar.set_label(t)

# And now we're done!
fig.savefig("slices_y_2.pdf", bbox_inches='tight')

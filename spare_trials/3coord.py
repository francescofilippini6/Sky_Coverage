from mw_plot import MWPlot

from galpy.potential import MWPotential2014
from galpy.orbit import Orbit
import numpy as np
from astropy import units as u

# Orbit Integration using galpy for the Sun
op = Orbit([0., 0., 0., 0., 0., 0.], radec=True, ro=8., vo=220.)
ts = np.linspace(0, 0.5, 10000) * u.Gyr
op.integrate(ts, MWPotential2014)
x = op.x(ts) * u.kpc
y = op.y(ts) * u.kpc
z = op.z(ts)

# setup a MWPlot instance with a certain center and radius
plot_instance = MWPlot(center=(-16, -2.5) * u.kpc, radius=5 * u.kpc)

# set up plot title
plot_instance.title = 'Orbit of Sun in 0.5 Gyr using galpy'

# plot, need to subtract 8kpc to shift to galactic coordinates in right hands frame
plot_instance.plot(x - 8. * u.kpc, y, c='r', linewidth=8.0)

# Save the figure
plot_instance.savefig(file='mw_plot_zoomed.png')

# Show the figure
plot_instance.show()

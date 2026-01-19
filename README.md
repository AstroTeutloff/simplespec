# SIMPLESPEC
A tiny (and growing, depending on my needs) python module for analyses of
Spectra in python.

## INSTALLATION
You can install the package via `pip`:
```
pip install git+https://github.com/AstroTeutloff/simplespec.git
```
or by cloning the repository:
```bash
git clone https://github.com/AstroTeutloff/simplespec.git
cd simplespec
pip install .
```

## USE
### PLOTTING
```python
import matplotlib.pyplot as plt
import simplespec as ssp
sdss_obj = ssp.SDSSSpectrum("YOUR_SDSS_SPEC.fits")
_ = sdss_obj.plot_spectrum(
  show_coadd = True,
  show_individual = True,
  show_uncertainty = True
)
plt.show()
```

### FITTING
The method `fit_continuum` fits a weighted spline to the spectral data. As of
the current time, this method is also still fitting the transition lines. In
future versions of this code (TBD) a linelist can be provided to the method to
exclude specific transition lines.

```python
# Continued from above
spline_continuum = sdss_obj.fit_continuum(sdss_obj.coadd, 1e5)
```

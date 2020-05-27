# pyme-fsc
A Fourier shell (ring) correlation recipe module written for [PYME](https://python-microscopy.org/).

Supports calculation from localization or image data.


## System requirements
* Windows, Linux or OS X
* Python 2.7
* PYME (>18.7.18) and dependencies

- Tested on Windows 10 with PYME (18.7.18)


## Installation

1. Clone repository.
2. Run the following in the project folder. 
	```
		python setup.py develop
	```
3. Start `dh5view` or `VisGUI` (PYME).

(Runtime < 5 min)


## Demo

### For localization data:
1. Open this simulated dataset ([wormlike_simulated.h5r](/fsc/example/wormlike_simulated.h5r)) with `VisGUI` (PYME).
2. Load and run this demo recipe ([frc_from_locs.yaml](/fsc/example/frc_from_locs.yaml)).
3. The results will display in a new window.

(Runtime < 5 min)

### For image data:
1. Open this simulated image ([wormlike_simulated_img_pair.h5](/fsc/example/wormlike_simulated_img_pair.h5)) with `dh5view` (PYME).
2. Load and run this demo recipe ([frc_from_images.yaml](/fsc/example/frc_from_images.yaml)).
3. The results will display in a new window.

(Runtime < 5 min)


## Instructions

### For localization data:
1. Refer to [PYME documention](https://python-microscopy.org/doc/index.html) for general use of PYME.
2. Open localization file with `VisGUI` (PYME).
3. Add the module **FSCFromLocs**.
	Detailed description the module and its inputs, outputs and parameters are accessible in PYME.
4. The recipe should autorun.
5. The results will display in a new window.

### For image data:
1. Refer to [PYME documention](https://python-microscopy.org/doc/index.html) for general use of PYME.
2. Open image file with `dh5view` (PYME).
3. Add the module **FSCFromImages**.
	Detailed description the module and its inputs, outputs and parameters are accessible in PYME.
4. Need to at least provide the file path of the 2nd image (`image_b_path`).
5. Run recipe.
6. The results will display in a new window.


## To do's
* Find method to calculate resolution per (pair of) dimension.
* Implement support for running multiply times to estimate errors.

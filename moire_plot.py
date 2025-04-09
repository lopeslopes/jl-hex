import numpy as np
from pathlib import Path
import holoviews as hv
from holoviews.operation.datashader import rasterize
from PIL import Image, ImageChops

hv.extension('bokeh')


def trim(image_path, output_path=None, border_color=(255, 255, 255)):
    img = Image.open(image_path).convert("RGB")
    bg = Image.new("RGB", img.size, border_color)
    diff = ImageChops.difference(img, bg)
    bbox = diff.getbbox()
    if bbox:
        cropped = img.crop(bbox)
        if output_path:
            cropped.save(output_path)
        return cropped
    else:
        print("No borders found (image likely fully solid or transparent).")
        return img


# -- Load data
def load_lattice(path):
    with open(path) as f:
        return np.array([
            [np.float64(x) for x in line.split(";")]
            for line in f if line.strip()  # skip empty lines
        ])


# -- Load properties
def load_properties(path):
    p = 0
    q = 0
    i = 0
    steps = 0
    with open(path / "properties.dat") as f:
        data = f.readlines()
        for line in data:
            if line != "\n":
                aux = line.split("=")
                if aux[0] == "p":
                    p = int(float(aux[1]))
                elif aux[0] == "q":
                    q = int(float(aux[1]))
                elif aux[0] == "i":
                    i = int(float(aux[1]))
                elif aux[0] == "steps":
                    steps = int(float(aux[1]))
    return p, q, i, steps


base_data_path = Path("data")
dataset_dirs = sorted([d for d in base_data_path.iterdir() if d.is_dir()])
n_datasets = len(dataset_dirs)
k = n_datasets-1
print(dataset_dirs[k])

a = 2.46
p, q, j, steps = load_properties(dataset_dirs[k])
angle_i = np.acos((3.0*(q**2) - (p**2))/(3.0*(q**2) + (p**2)))
angle_f = np.acos((3.0*((q-1)**2) - (p**2))/(3.0*((q-1)**2) + (p**2)))
angle = angle_i + (j/steps)*(angle_f - angle_i)
print("Angle in radians: ", angle)
print("Angle in degrees: ", (angle * 180) / np.pi)
moire_period = a/(2*np.sin(angle/2))
print("D = ", moire_period)

latA1 = load_lattice(dataset_dirs[k] / "latticeA1.dat")
latB1 = load_lattice(dataset_dirs[k] / "latticeB1.dat")
latA2 = load_lattice(dataset_dirs[k] / "latticeA2.dat")
latB2 = load_lattice(dataset_dirs[k] / "latticeB2.dat")

latAA = load_lattice(dataset_dirs[k] / "latticeAA.dat")
latBA = load_lattice(dataset_dirs[k] / "latticeBA.dat")
latAB = load_lattice(dataset_dirs[k] / "latticeAB.dat")
latBB = load_lattice(dataset_dirs[k] / "latticeBB.dat")

# -- Combine data
r = moire_period
theta = np.linspace(0, 2 * np.pi, 200)
x = r * np.cos(theta)
y = r * np.sin(theta)
circle_data = np.stack([x, y], axis=1)

r = moire_period*2.0
theta = np.linspace(0, 2 * np.pi, 200)
x = r * np.cos(theta)
y = r * np.sin(theta)
circle_data2 = np.stack([x, y], axis=1)


### holoviews plot ###
all_pts = np.concatenate((latA1, latB1, latA2, latB2))
points_total = hv.Points(all_pts)

raster1 = rasterize(points_total, width=800, height=800)
raster1 = raster1.opts(
    cmap='bmy',
    colorbar=True,
    frame_width=800,
    title="Moire Pattern Visualization",
    aspect="equal"
)

pointsAB = hv.Points(latAB).opts(size=10)
pointsBA = hv.Points(latBA).opts(size=10)
circ1 = hv.Curve(circle_data)
circ2 = hv.Curve(circle_data2)

overlap = pointsAB * pointsBA * circ1 * circ2

overlay = raster1 * overlap
hv.save(overlay, 'my_plot.png', fmt='png', backend='bokeh')
trim('my_plot.png', 'my_plot.png')

import numpy as np
from pathlib import Path
import holoviews as hv
from holoviews.operation.datashader import rasterize
from PIL import Image, ImageChops
import logging
logging.getLogger('bokeh').setLevel(logging.FATAL)


hv.extension("bokeh")


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


# load data
def load_lattice(path):
    with open(path) as f:
        return np.array([
            [np.float64(x) for x in line.split(";")]
            for line in f if line.strip()
        ])


# load properties
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


angles = []
AA = []
BA = []
AB = []
BB = []
with open("AAstack_separations") as f:
    data = f.readlines()
    for line in data:
        if line[0:2] == "0.":
            angles.append(float(line))
        elif line[0:3] == "AA:":
            aux = line.split(":")
            AA.append(float(aux[1]))
        elif line[0:3] == "BA:":
            aux = line.split(":")
            BA.append(float(aux[1]))
        elif line[0:3] == "AB:":
            aux = line.split(":")
            AB.append(float(aux[1]))
        elif line[0:3] == "BB:":
            aux = line.split(":")
            BB.append(float(aux[1]))


min_index = np.argmin(AB)
print("Angle:      ", angles[min_index])
print("Separation: ", AB[min_index])


# load available data folders
base_data_path = Path("data")
dataset_dirs = sorted([d for d in base_data_path.iterdir() if d.is_dir()])

print(dataset_dirs[min_index+1])



k = min_index+1
# obtain lattice properties
a = 2.46
p, q, j, steps = load_properties(dataset_dirs[k])
angle_i = np.acos((3.0*(q**2) - (p**2))/(3.0*(q**2) + (p**2)))
angle_f = np.acos((3.0*((q-1)**2) - (p**2))/(3.0*((q-1)**2) + (p**2)))
angle = angle_i + (j/steps)*(angle_f - angle_i)
print("Angle in radians: ", angle)
print("Angle in degrees: ", (angle * 180) / np.pi)
moire_period = a/(2*np.sin(angle/2))
print("D = ", moire_period)

# load data from the lattices
latA1 = load_lattice(dataset_dirs[k] / "latticeA1_dist.dat")
latB1 = load_lattice(dataset_dirs[k] / "latticeB1_dist.dat")
latA2 = load_lattice(dataset_dirs[k] / "latticeA2.dat")
latB2 = load_lattice(dataset_dirs[k] / "latticeB2.dat")

latAA = load_lattice(dataset_dirs[k] / "latticeAA_dist.dat")
latBA = load_lattice(dataset_dirs[k] / "latticeBA_dist.dat")
latAB = load_lattice(dataset_dirs[k] / "latticeAB_dist.dat")
latBB = load_lattice(dataset_dirs[k] / "latticeBB_dist.dat")

# moire period circle visualization
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


# holoviews plot
folder_name = str(dataset_dirs[k])
first_underline = folder_name.find("_")
second_underline = folder_name.find("_", first_underline+1)
angle_name = folder_name[first_underline+1:second_underline]

all_pts = np.concatenate((latA1, latB1, latA2, latB2))
points_total = hv.Points(all_pts)

min_x = -850
max_x = 850
min_y = -850
max_y = 850
f_width = max_x-min_x
f_height = max_y-min_y
moire = rasterize(points_total, width=f_width+100, height=f_height+100)
moire = moire.opts(
    cmap="blues",
    # colorbar=True,
    yaxis=None,
    frame_width=f_width,
    frame_height=f_height,
    xlim=(min_x, max_x),
    ylim=(min_y, max_y),
    title=angle_name,
    aspect="equal"
)

overlap = hv.Overlay([])
if (latAB.size > 0):
    pointsAB = hv.Points(latAB, label="AB").opts(size=10, color="purple")
    overlap *= pointsAB
if (latBA.size > 0):
    pointsBA = hv.Points(latBA, label="BA").opts(size=10, color="orange")
    overlap *= pointsBA
if (latAA.size > 0):
    pointsAA = hv.Points(latAA, label="AA").opts(size=10, color="green")
    overlap *= pointsAA
if (latBB.size > 0):
    pointsBB = hv.Points(latBB, label="BB").opts(size=10, color="magenta")
    overlap *= pointsBB

circ1 = hv.Curve(circle_data).opts(color="gray")
circ2 = hv.Curve(circle_data2).opts(color="gray")

overlap = overlap * circ1 * circ2
overlap = overlap.opts(show_legend=True)

full_image = moire * overlap
hv.save(full_image, angle_name+".png", fmt="png", backend="bokeh")
# trim(angle_name+".png", angle_name+".png")

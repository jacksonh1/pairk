import os
from pathlib import Path

src_dir = Path(os.path.dirname(__file__)).parent

colormap_dir = src_dir / 'local_data/color_maps/'

COLOR_MAP_FILES = {
    'clustal': colormap_dir / 'clustal_hex_map.json',
}
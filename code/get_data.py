"""Download the three public datasets into the layout expected by c1_data.py (no data are shipped with this repository).

Usage: python get_data.py <dataset_dir>        then set THS_DATASETS=<dataset_dir> before running the experiments.

Layout produced:
  <dataset_dir>/UCI_AirQuality/AirQualityUCI.csv     UCI #360 Air Quality (CC BY 4.0)
  <dataset_dir>/UCI_GasDrift/conc/batch1..10.dat      UCI #270 Gas Sensor Array Drift at Different Concentrations (CC BY 4.0)
  <dataset_dir>/IntelLab/data.txt.gz, mote_locs.txt   Intel Berkeley Research Lab (no stated licence; do not redistribute)
Licences: read the dataset pages before use. The UCI pages state that the data "can be used exclusively for research purposes".
"""
import io, os, sys, urllib.request, zipfile

UA = {"User-Agent": "aiot-h35-drift/1.0 (research reproduction script)"}
SRC = {
    "air": "https://archive.ics.uci.edu/static/public/360/air+quality.zip",
    "gas": "https://archive.ics.uci.edu/static/public/270/gas+sensor+array+drift+dataset+at+different+concentrations.zip",
    "intel": ["http://db.csail.mit.edu/labdata/data.txt.gz", "http://db.csail.mit.edu/labdata/mote_locs.txt"],
}


def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r:
        return r.read()


def extract(blob, want, out_dir):
    """Extract members whose base name matches `want` (callable) from a zip, descending into nested zips."""
    os.makedirs(out_dir, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        for m in z.namelist():
            base = os.path.basename(m)
            if base.lower().endswith(".zip"):
                extract(z.read(m), want, out_dir)
            elif base and want(base):
                with open(os.path.join(out_dir, base), "wb") as f:
                    f.write(z.read(m))


def main(root):
    extract(get(SRC["air"]), lambda b: b == "AirQualityUCI.csv", os.path.join(root, "UCI_AirQuality"))
    extract(get(SRC["gas"]), lambda b: b.startswith("batch") and b.endswith(".dat"), os.path.join(root, "UCI_GasDrift", "conc"))
    d = os.path.join(root, "IntelLab"); os.makedirs(d, exist_ok=True)
    for u in SRC["intel"]:
        with open(os.path.join(d, u.rsplit("/", 1)[1]), "wb") as f:
            f.write(get(u))
    print("done:", root)


if __name__ == "__main__":
    main(sys.argv[1])

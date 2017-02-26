import glob
from pyPanair.preprocess import read_wgs


if __name__ == '__main__':
    print("converting LaWGS to stl")
    wgs_list = glob.glob("*.wgs")
    for w in wgs_list:
        wgs = read_wgs(w)
        wgs.create_stl(w.replace(".wgs", ".stl"))
    print("success!")
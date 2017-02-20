from pyPanair.preprocess import read_wgs

if __name__ == '__main__':
    print("converting LaWGS to stl")
    wgs = read_wgs("ADODG3.wgs", boun_cond=(1,1,18))
    wgs.create_stl("ADODG3.stl")
    print("success!")
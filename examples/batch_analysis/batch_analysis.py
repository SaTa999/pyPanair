#!/usr/bin/env python
# -*- coding: utf-8 -*-
from glob import glob
from logging import getLogger, StreamHandler, FileHandler, Formatter, INFO
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool
from os import remove, makedirs
from os.path import isdir, join
from shutil import copy2, rmtree
from subprocess import Popen, PIPE
from numpy import genfromtxt
import makewgs_aux_model3


def panin(folder, aux):
    process = Popen("./panin", stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=folder)
    aux = str.encode(aux)
    stdout, stderr = process.communicate(aux)
    if stderr:
        logger.error("\n".join((folder, stderr)))
    else:
        logger.debug("\n".join((folder, stdout)))


def panair(folder):
    process = Popen("./panair", stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=folder)
    stdout, stderr = process.communicate(b"a502.in")
    if "PanAir is stopping." in open(join(folder, "panair.err")).read():
        logger.critical("\n".join((folder, "fatal error with PanAir")))
    if stderr:
        logger.error("\n".join((folder, stderr)))
    else:
        logger.debug("\n".join((folder, stdout)))


def safe_makedirs(path):
    try:
        makedirs(path)
    except OSError:
        if not isdir(path):
            logger.error("could not create directory at {}".format(path))


def copy2folder(files, oldfolder, newfolder):
    safe_makedirs(newfolder)
    if type(files) is str:
        files = (files)
    for f in files:
        try:
            copy2("{0}/{1}".format(oldfolder, f), "{0}/{1}".format(newfolder, f))
        except FileNotFoundError as e:
            logger.error(e)

def clean(root, casename):
    # 前回使用した中間ファイルを削除
    folders = [f for f in glob(join(root, "panair[0-9]")) if isdir(f)]
    for f in folders:
        try:
            rmtree(f)
        except WindowsError:
            logger.warning("failed to delete {}".format(f))

    # auxを削除
    subcases = glob(join(root, "{}_[0-9].aux".format(casename)))
    for f in subcases:
        try:
            remove(f)
        except WindowsError:
            logger.warning("failed to delete {}".format(f))


if __name__ == '__main__':
    # 初期設定
    casename = "model3"
    wgsname = casename + ".wgs"
    enable_ramdisk = False # ramdiskを使用する場合はTrue
    if enable_ramdisk:
        root = "Z:/" # ramdiskのpath
    else:
        root = "./"

    # ロガーの設定
    logger = getLogger(__name__)
    shandler = StreamHandler()
    shandler.setFormatter(Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s"))
    fhandler = FileHandler(filename="panelmethod.log")
    fhandler.setFormatter(Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s"))
    shandler.setLevel(INFO)
    fhandler.setLevel(INFO)
    logger.setLevel(INFO)
    logger.addHandler(shandler)
    logger.addHandler(fhandler)
    logger.info("started program")

    # caselistを読み込む
    caselist = genfromtxt("caselist.csv", delimiter=',', names=True)
    for case in caselist:
        dvs = {"xfin":case["x"], "yfin":case["y"], "delta":case["delta"],
               "span_fin":case["b"], "rchord_fin":case["cr"],
               "taper":case["taper"], "sweep":case["sweep"]}
        try:
            dvs["zfin"] = case["z"]
        except ValueError:
            logger.debug("zfin is not specified")
        if case["flag"] == 1:
            logger.debug("skipping case{}".format(int(case["casenum"])))
            continue
        print(dvs)

        clean(root, casename)

        logger.info("calculating case{}".format(int(case["casenum"])))
        makewgs_aux_model3.main(**dvs)
        subcases = glob("{}_[0-9].aux".format(casename))
        subcase_folders = [join(root, "panair{}".format(i+1)) for (i, _) in enumerate(subcases)]
        # paninを起動
        num = min(len(subcases),cpu_count())  # 1度に使用するスレッド数（Noneなら上限まで使用）
        tp = ThreadPool(num)
        for (subcase, folder) in zip(subcases, subcase_folders):
            copy2folder((subcase, wgsname), ".", folder)
            tp.apply_async(panin, (folder, subcase, ))
        tp.close()
        tp.join()

        # panairを起動
        tp = ThreadPool(num)
        for folder in subcase_folders:
            tp.apply_async(panair, (folder,))
        tp.close()
        tp.join()

        # 結果を保存
        files_to_save = ("panair.out", "ffmf", "agps", "panair.err", "panin.dbg", wgsname)
        for folder in subcase_folders:
            newfolder = join("results", "case{0}".format(int(case["casenum"])), folder.lstrip(root))
            copy2folder(files_to_save, folder, newfolder)


    logger.info("finished program")

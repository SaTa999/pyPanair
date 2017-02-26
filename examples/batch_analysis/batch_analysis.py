#!/usr/bin/env python
from concurrent import futures
from logging import getLogger, StreamHandler, FileHandler, Formatter, INFO, DEBUG
import multiprocessing
import os
from shutil import copy2, rmtree
from subprocess import Popen, PIPE
import pandas as pd
from make_wgs_aux import main as mkwgsaux


# initialize logger
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
logger.info("start batch analysis")


def run_panin(directory, aux):
    process = Popen("./panin", stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=directory)
    aux = str.encode(aux)
    stdout, stderr = process.communicate(aux)
    if stderr:
        logger.error("\n".join((directory, str(stderr))))
    else:
        logger.debug("\n".join((directory, str(stdout))))


def run_panair(directory):
    process = Popen("./panair", stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = process.communicate(b"a502.in")
    if "PanAir is stopping." in open(os.path.join(directory, "panair.err")).read():
        logger.critical("\n".join((directory, "fatal error with PanAir")))
    if stderr:
        logger.error("\n".join((directory, str(stderr))))
    else:
        logger.debug("\n".join((directory, str(stdout))))


def safe_makedirs(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            logger.error("could not create directory at {}".format(path))


def copy2dir(files, olddir, newdir):
    safe_makedirs(newdir)
    if type(files) is str:
        files = (files, )
    for f in files:
        try:
            copy2("{0}/{1}".format(olddir, f), "{0}/{1}".format(newdir, f))
        except FileNotFoundError as e:
            logger.error(e)


def run_analysis(casenum, auxname, analysis_dir, params):
    logger.info("calculating case{}".format(casenum))
    # create directory to run panin and panair
    procid = int(multiprocessing.current_process().pid)
    target_dir = os.path.join("panair{}".format(procid))
    safe_makedirs(target_dir)

    # run panin and panair
    mkwgsaux(target_dir=target_dir, **params)
    run_panin(target_dir, auxname)
    run_panair(target_dir)

    # save results
    files_to_save = ("panair.out", "ffmf", "agps", "panair.err", "panin.dbg", "a502.in")
    newdir = os.path.join("results", "case{}".format(casenum))
    copy2dir(files_to_save, target_dir, newdir)

    # delete intermediate files
    rmtree(target_dir)


if __name__ == '__main__':
    # initialize
    N_PROCS = 3
    wgsname = "ADODG_case3.wgs"
    auxname = "ADODG_case3.aux"
    analysis_dir = "" # directory to run analysis (intermediate files will be stored here)

    # read caselist
    caselist = pd.read_csv("caselist.csv")
    with futures.ProcessPoolExecutor(max_workers=N_PROCS) as executor:
        # submit jobs
        fs = list()
        for _, case in caselist.iterrows():
            casenum = int(case["casenum"])
            params = case[2:].to_dict()
            if case["run"] == 1:
                fs.append(executor.submit(run_analysis, casenum, auxname, analysis_dir, params))
            else:
                logger.debug("skipping case{}".format(casenum))
                continue

    # run analysis
    for future in futures.as_completed(fs):
        pass

    logger.info("finish batch analysis")
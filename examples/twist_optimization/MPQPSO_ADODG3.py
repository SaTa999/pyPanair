#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
import seaborn as sns
from deap import base, benchmarks, creator, tools
from concurrent import futures
import multiprocessing
from logging import getLogger, StreamHandler, FileHandler, Formatter, INFO
import os
import shutil

from ADODG_case3 import main as adc3
import clean_rwms
from pyPanair.utilities import bspline
from pyPanair.postprocess import read_ffmf


def latin_hypercube_sampling(size, n):
    """standard latin hypercube sampling
    size: number of dimensions
    n: number of samples"""
    cut = np.linspace(0, 1, n + 1)
    u = np.random.rand(n, size)
    a = cut[:n, np.newaxis]
    b = cut[1:n + 1, np.newaxis]
    rdpoints = u * (b - a) + a
    # Make the random pairings
    h = np.zeros_like(rdpoints)
    for j in range(size):
        order = np.random.permutation(list(range(n)))
        h[:, j] = rdpoints[order, j]
    return h


def lhs_swarm(generator, size, n, pmin, pmax):
    """generate a swarm that is coposed of particles distributed in a latin hypercube
    generator: generator inheriting deap.base.Fitness
    size: number of dimensions
    n: number of samples
    pmin, pmax: min & max particle position"""
    h = latin_hypercube_sampling(size, n)
    h *= (pmax - pmin)
    h += pmin
    population = [generator(position=x) for x in h]
    return population


def generate(creator, size, pmin, pmax, position=None):
    """generate a particle
    creator: creator inheriting deap.base.Fitness
    size: number of dimensions
    pmin, pmax: lower & upper bound of the position of the particle
    position: initial position of the particle"""
    if position is None:
        part = creator(np.random.uniform(pmin, pmax, size))
    else:
        part = creator(position)
    part.pmin = pmin
    part.pmax = pmax
    return part


def updateParticle(part, best, meanbest, beta):
    """update the particle's position (type1: quantum behaved PSO by Sun 2004)"""
    phi = np.random.random(len(part))
    p = phi * part.best + (1. - phi) * best
    k = np.random.randint(0, 2, len(part)) * 2 - 1  # random array of -1 & 1
    u = np.random.random(len(part))
    #     part[:] = p + k * beta * np.abs(meanbest - part) * np.log(1./u)
    part[:] = p + np.einsum("i,...,i,i->i", k, beta, np.abs(meanbest - part), np.log(1. / u))
    part[:] = np.clip(part, part.pmin, part.pmax)
    return part


def updateParticle2(part, best, meanbest, beta):
    """update the particle's position (type2: gaussian QPSO by Sun 2011)
    increases diversity of particle position (lower chance of getting trapped in local optimum)"""
    phi = np.random.random(len(part))
    p = phi * part.best + (1. - phi) * best
    k = np.random.randint(0, 2, len(part)) * 2 - 1  # random array of -1 & 1
    u = np.random.random(len(part))
    p += np.random.normal(0, 1, len(part)) * np.sqrt(np.abs(meanbest - part.best))
    #     part[:] = p + k * beta * np.abs(meanbest - part) * np.log(1./u)
    part[:] = p + np.einsum("i,...,i,i->i", k, beta, np.abs(meanbest - part), np.log(1. / u))
    if part.pmin or part.pmax:
        part[:] = np.clip(part, part.pmin, part.pmax)
    return part


def sort_cp(part):
    col = (part.shape[0] - 1) // 2
    tmp = part[:-1].reshape(2, col)
    tmp[:] = tmp[:, tmp[0, :].argsort()]
    part[:-1] = tmp.ravel()


def run_panin(directory):
    process = Popen("panin", stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=directory)
    (stdout, stderr) = process.communicate(b"ADODG_case3.aux")
    if stderr:
        logger.error(stderr)
    else:
        logger.debug(stdout)


def run_panair(directory):
    process = Popen("panair", stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=directory)
    (stdout, stderr) = process.communicate(b"a502.in")
    if stderr:
        logger.error(stderr)
    else:
        logger.debug(stdout)


def evaluate(part):
    procid = int(multiprocessing.current_process().pid)
    target_dir = "panair{}".format(procid)
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)
    aoa_low = 0
    aoa_high = 10
    target_cl = 37.5
    auxpath = "{}/ADODG_case3.aux".format(target_dir)
    ffmfpath = "{}/ffmf".format(target_dir)

    nan = (float("nan"),)

    # define twist func and create aux & wgs
    n_cp = (part.shape[0] - 1) // 2
    cv = np.zeros((2, n_cp + 2))
    cv[:, 1:-1] = part[:-1].reshape(2, n_cp)
    cv[0, -1] = 1
    cv[1, -1] = part[-1]
    twist_func = bspline(cv.T, degree=3, periodic=False)
    adc3(twist_func, aoas=aoa_low, target_dir=target_dir)

    # run panair
    run_panin(target_dir)
    run_panair(target_dir)

    # read ffmf for aoa_low
    try:
        ffmf = read_ffmf(ffmfpath)
        cllow = float(ffmf.cl[0])
    except FileNotFoundError:
        return nan

    # modify aux file
    with open(auxpath, "r") as f:
        aux = f.readlines()
        aux[3] = "ALPHA {}\n".format(aoa_high)
        aux[4] = "ALPC {}\n".format(aoa_high)
        aux = "".join(aux)
    with open(auxpath, "w") as f:
        f.write(aux)

    # clean directory
    clean_rwms.main(target_dir)
    os.remove(ffmfpath)

    # run panair
    run_panin(target_dir)
    run_panair(target_dir)

    # read ffmf for aoa_low
    try:
        ffmf = read_ffmf(ffmfpath)
        clhigh = float(ffmf.cl[0])
    except FileNotFoundError:
        return nan

    # calc target aoa (CL=0.375)
    cla = (clhigh - cllow) / (aoa_high - aoa_low)
    cl0 = cllow - aoa_low * cla
    target_alpha = (target_cl - cl0) / cla

    # modify aux file
    with open(auxpath, "r") as f:
        aux = f.readlines()
        aux[3] = "ALPHA {}\n".format(target_alpha)
        aux[4] = "ALPC {}\n".format(target_alpha)
        aux = "".join(aux)
    with open(auxpath, "w") as f:
        f.write(aux)

    # clean directory
    clean_rwms.main(target_dir)
    os.remove(ffmfpath)

    # run panair
    run_panin(target_dir)
    run_panair(target_dir)

    try:
        ffmf = read_ffmf(ffmfpath)
    except FileNotFoundError:
        return nan
    cl = float(ffmf.cl)
    cdi = float(ffmf.cdi)

    shutil.rmtree(target_dir)

    if cl < 37.5:
        return nan
    else:
        return cdi,


# set logger
logger = getLogger(__name__)
shandler = StreamHandler()
shandler.setFormatter(Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s"))
fhandler = FileHandler(filename="QPSO.log")
fhandler.setFormatter(Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s"))
shandler.setLevel(INFO)
fhandler.setLevel(INFO)
logger.setLevel(INFO)
logger.addHandler(shandler)
logger.addHandler(fhandler)

if __name__ == '__main__':
    logger.info("started program")
    N_PROCS = 10
    NCP = 4  # number of control points
    NDV = NCP * 2 + 1  # number of design variables
    NPOP = 100  # number of particles (more than 3 to 4 times the number of design variables)
    GEN = 200  # number of maximum generations
    LOWER = np.array((0.,) * NCP + (-10.,) * (NCP + 1))  # lower bound of particle position x1 x2 x3 x4 y1 y2 y3 y4 y5
    UPPER = np.array((1.,) * NCP + (5.,) * (NCP + 1))  # upper bound of particle position
    FITNESS_WEIGHTS = (-1.0,)  # maximize:(1.0,), minimize:(-1.0,), don't forget the comma!
    BETA_INIT = 1.0  # initial value of contraction-expansion coefficient
    BETA_FIN = 0.5  # final value of contraction-expansion coefficient
    """TIPS: large beta -> global search, small beta -> local search
             beta_init: 0.8 to 1.2
             beta_fin : less than 0.6
             beta must be below e^gamma=1.781 to guarantee convergence of the particle"""
    LOG = True  # record logbook & graph

    # register creators
    creator.create("FitnessMax", base.Fitness, weights=FITNESS_WEIGHTS)
    creator.create("Particle", np.ndarray, fitness=creator.FitnessMax,
                   pmin=None, pmax=None, best=None)

    # register functions to the toolbox
    toolbox = base.Toolbox()
    toolbox.register("particle", generate, creator=creator.Particle, size=NDV, pmin=LOWER, pmax=UPPER)
    # toolbox.register("population", tools.initRepeat, list, toolbox.particle)
    toolbox.register("population", lhs_swarm, generator=toolbox.particle, size=NDV, pmin=LOWER, pmax=UPPER)
    toolbox.register("update", updateParticle)

    # register the functions used to calculate stats
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("mean", np.nanmean)
    stats.register("std", np.nanstd)
    stats.register("min", np.nanmin)
    stats.register("max", np.nanmax)

    # initialize the logbook
    logbook = tools.Logbook()
    logbook.header = ["gen", "evals", "gbest", "invalid"] + stats.fields

    pop = toolbox.population(n=NPOP)
    best = None
    betas = np.linspace(BETA_INIT, BETA_FIN, GEN)
    meanbest = np.empty(len(pop[0]))  # initialize the meanbest
    fitness_all = np.empty((NPOP,))

    with futures.ProcessPoolExecutor(max_workers=N_PROCS) as executor:

        for g, beta in zip(list(range(GEN)), betas):
            meanbest[:] = 0.  # reinitialize the meanbest
            INVALID = 0
            # sort the x-coordinates of the contril points
            for part in pop:
                sort_cp(part)

            # MP region
            mappings = {executor.submit(evaluate, np.array(part)): i for (i, part) in enumerate(pop)}
            for future in futures.as_completed(mappings):
                target = mappings[future]
                fitness_all[target] = future.result()[0]
            # MP region end

            for i, part in enumerate(pop):
                part.fitness.values = (fitness_all[i],)
                if part.fitness.values[0] is float("nan"):
                    INVALID += 1
                # reinitialize the particle if its initial fitness value is nan
                if part.best is None:
                    while np.isnan(part.fitness.values):
                        part[:] = toolbox.particle()
                        part.fitness.values = evaluate(part)
                # update the particle's best position
                if part.best is None or part.best.fitness < part.fitness:
                    part.best = creator.Particle(part)
                    part.best.fitness.values = part.fitness.values
                # update the global best position
                if best is None or best.fitness < part.fitness:
                    best = creator.Particle(part)
                    best.fitness.values = part.fitness.values
                meanbest += part.best
            meanbest /= NPOP

            # update the particle's position
            for part in pop:
                toolbox.update(part, best, meanbest, beta)

            # Gather all the fitnesses in one list and print the stats
            logbook.record(gen=g, evals=len(pop) * (g + 1), gbest=best.fitness.values[0],
                           invalid=INVALID, **stats.compile(pop))
            if LOG:
                print(logbook.stream)
            logger.info("end of generation {}".format(g))

    print("end of global optimization")
    print("best particle value is {}".format(best.fitness.values[0]))
    print("best particle location is {}".format(best))
    print()
    np.savetxt("res.txt", best)
    with open("resval.txt", "w") as f:
        f.write(str(best.fitness.values[0]))

    if LOG:
        # save logbook
        dframe = pd.DataFrame.from_dict(logbook)
        dframe.to_csv("logbook.csv", index=False)
        # plot the stats
        sns.set_context("paper", font_scale=1.8)
        sns.set_style("ticks", rc={"legend.frameon": True})
        # sns.set_palette("Set1", n_colors=9, desat=0.8)
        gen = logbook.select("gen")
        fit_mins = logbook.select("min")
        fit_avgs = logbook.select("mean")
        fit_stds = logbook.select("std")

        fig, ax1 = plt.subplots()
        line1 = ax1.plot(gen, fit_mins, "b-", label="Minimum Fitness")
        line2 = ax1.plot(gen, fit_avgs, "r-", label="Mean Fitness")
        ax1.set_xlabel("Generation")
        ax1.set_ylabel("Fitness (Min & Mean)")

        ax2 = ax1.twinx()
        ax2.set_ylabel("Fitness (Std)")
        line3 = ax2.plot(gen, fit_stds, "g-", label="Std Fitness")

        lns = line1 + line2 + line3
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, loc="center right")
        plt.tight_layout(pad=0.4)

        #     plt.show()
        plt.savefig("fitnees_log")

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
    b = cut[1:n+1, np.newaxis] 
    rdpoints = u * (b - a) + a
    # Make the random pairings
    H = np.zeros_like(rdpoints)
    for j in range(size):
        order = np.random.permutation(list(range(n)))
        H[:, j] = rdpoints[order, j]
    return H


def lhs_swarm(generator, size, n, pmin, pmax):
    """generate a swarm that is coposed of particles distributed in a latin hypercube
    generator: generator inheriting deap.base.Fitness
    size: number of dimensions
    n: number of samples
    pmin, pmax: min & max particle position"""
    H = latin_hypercube_sampling(size, n)
    H *= (pmax - pmin)
    H += pmin
    population = [generator(position=x) for x in H]
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
    part.pmin=pmin
    part.pmax=pmax 
    return part


def updateParticle(part, best, meanbest, beta):
    """update the particle's position (type1: quantum behaved PSO by Sun 2004)"""
    phi = np.random.random(len(part))
    p = phi * part.best + (1. - phi) * best
    k = np.random.randint(0,2,len(part)) * 2 - 1 # random array of -1 & 1
    u = np.random.random(len(part))
#     part[:] = p + k * beta * np.abs(meanbest - part) * np.log(1./u)
    part[:] = p + np.einsum("i,...,i,i->i",k,beta,np.abs(meanbest - part),np.log(1./u))
    part[:] = np.clip(part, part.pmin, part.pmax)
    return part


def updateParticle2(part, best, meanbest, beta):
    """update the particle's position (type2: gaussian QPSO by Sun 2011)
    increases diversity of particle position (lower chance of getting trapped in local optimum)"""
    phi = np.random.random(len(part))
    p = phi * part.best + (1. - phi) * best
    k = np.random.randint(0,2,len(part)) * 2 - 1 # random array of -1 & 1
    u = np.random.random(len(part))
    p += np.random.normal(0, 1, len(part)) * np.sqrt(np.abs(meanbest - part.best))
#     part[:] = p + k * beta * np.abs(meanbest - part) * np.log(1./u)
    part[:] = p + np.einsum("i,...,i,i->i",k,beta,np.abs(meanbest - part),np.log(1./u))
    if part.pmin or part.pmax:
        part[:] = np.clip(part, part.pmin, part.pmax)
    return part


def evaluate(part):
    procid = int(multiprocessing.current_process().pid)
    target_dir = "panair{}".format(procid)
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)

    # define twist func
    cv = np.zeros((4, 2))
    cv[1, :] = part[:2]
    cv[2, :] = part[2:4]
    cv[-1,0] = 1
    cv[-1,1] = part[-1]
    twist_func = bspline(cv, degree=3, periodic=False)
    aoa_low = 2
    aoa_high = 8
    adc3(twist_func, aoas=(aoa_low, aoa_high), target_dir=target_dir)
    # run panin
    process = Popen("panin", stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=target_dir)
    (stdout, stderr) = process.communicate(b"ADODG_case3.aux")
    if stderr:
        logger.error(stderr)
    else:
        logger.debug(stdout)

    # run panair
    process = Popen("panair", stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=target_dir)
    (stdout, stderr) = process.communicate(b"a502.in")
    if stderr:
        logger.error(stderr)
    else:
        logger.debug(stdout)

    # calc target aoa (CL=0.375)
    ffmf = read_ffmf("{}/ffmf".format(target_dir))
    cllow = float(ffmf.cl[0])
    clhigh = float(ffmf.cl[1])
    cla = (clhigh - cllow) / (aoa_high - aoa_low)
    cl0 = clhigh - aoa_high * cla
    target_alpha = (37.5 - cl0) / cla
    clean_rwms.main(target_dir)

    # modify aux file
    with open("{}/ADODG_case3.aux".format(target_dir), "r") as f:
        aux = f.readlines()
        aux[3] = "ALPHA {}\n".format(target_alpha)
        aux = "".join(aux)
    with open("{}/ADODG_case3.aux".format(target_dir), "w") as f:
        f.write(aux)

    # run panin & panair for target_alpha
    process = Popen("panin", stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=target_dir)
    (stdout, stderr) = process.communicate(b"ADODG_case3.aux")
    if stderr:
        logger.error(stderr)
    else:
        logger.debug(stdout)
    process = Popen("panair", stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=target_dir)
    (stdout, stderr) = process.communicate(b"a502.in")
    if stderr:
        logger.error(stderr)
    else:
        logger.debug(stdout)
    clean_rwms.main(target_dir)

    # evaluate cl & cdi at target_alpha
    ffmf = read_ffmf("{}/ffmf".format(target_dir))
    cl = float(ffmf.cl)
    cdi = float(ffmf.cdi)
    if cl < 37.5:
        return (float("nan"), )
    else:
        return (cdi, )

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
    NDV = 5 # number of design variables
    NPOP = 50 # number of particles (more than 3 to 4 times the number of design variables)
    GEN = 200 # number of maximum generations
    LOWER = np.array((0, -5, 0.5, -5, -5))  # lower bound of particle position x1 y1 x2 y2 y3
    UPPER = np.array((0.5, 5, 1., 5, 5))  # upper bound of particle position
    FITNESS_WEIGHTS = (-1.0,) # maximize:(1.0,), minimize:(-1.0,), don't forget the comma!
    BETA_INIT = 1.0 # initial value of contraction-expansion coefficient
    BETA_FIN = 0.5 # final value of contraction-expansion coefficient
    """TIPS: large beta -> global search, small beta -> local search
             beta_init: 0.8 to 1.2
             beta_fin : less than 0.6
             beta must be below e^gamma=1.781 to guarantee convergence of the particle"""
    LOG = True # record logbook & graph

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
    meanbest = np.empty(len(pop[0])) # initialize the meanbest
    fitness_all = np.empty((NPOP,))

    with futures.ProcessPoolExecutor(max_workers=N_PROCS) as executor:

        for g, beta in zip(list(range(GEN)), betas):
            meanbest[:] = 0. # reinitialize the meanbest
            INVALID = 0

            # MP region
            mappings = {executor.submit(evaluate, np.array(part)): i for (i, part) in enumerate(pop)}
            for future in futures.as_completed(mappings):
                target = mappings[future]
                fitness_all[target] = future.result()[0]
            # MP region end

            for i, part in enumerate(pop):
                part.fitness.values = (fitness_all[i], )
                if part.fitness.values[0] is float("nan"):
                    INVALID += 1
                # reinitialize the particle if its initial fitness value is nan
                if part.best is None:
                    while np.isnan(part.fitness.values):
                        part[:] = toolbox.particle()
                        part.fitness.values = toolbox.evaluate(part)
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
            logbook.record(gen=g, evals=len(pop)*(g+1), gbest=best.fitness.values[0], invalid=INVALID,  **stats.compile(pop))
            if LOG:
                print((logbook.stream))
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
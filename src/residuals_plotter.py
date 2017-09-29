#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import csv
import time as tt


def plot(logfile):

    plt.ion()

    print "logfile", logfile

    fig, (resp, contp) = plt.subplots(nrows=2)
    iterp = resp.twinx()

    resp.set_xlabel("Time")
    resp.set_ylabel("Vel Residuals")
    resp.set_yscale("log")

    iterp.set_ylabel("NumIters")

    contp.set_ylabel("Pres Cont Residuals")
    contp.set_yscale("log")

    i = 0

    while (True):

        tt = []
        uresInit = []
        uresFin = []
        uresInit = []
        uresFin = []
        vresInit = []
        vresFin = []
        contresInit = []
        contresFin = []
        presresInit = []
        presresFin = []
        iters = []

        with open(logfile, 'rb') as csvfile:

            reader = csv.DictReader(csvfile)

            for row in reader:
                tt.append(row["Time"])
                uresInit.append(row["U-Res-Initial"])
                uresFin.append(row["U-Res-Final"])
                vresInit.append(row["V-Res-Initial"])
                vresFin.append(row["V-Res-Final"])
                contresInit.append(row["Cont-Res-Initial"])
                contresFin.append(row["Cont-Res-Final"])
                presresInit.append(row["P-Res-Initial"])
                presresFin.append(row["P-Res-Final"])
                iters.append(row["NumIters"])

        resp.plot(
            tt,
            uresInit,
            label='U-Res-Initial',
            marker='*',
            ls="--",
            color='r')
        resp.plot(
            tt, uresFin, label='U-Res-Final', marker='o', ls="-.", color='r')
        resp.plot(
            tt,
            vresInit,
            label='V-Res-Initial',
            marker='*',
            ls="--",
            color='b')
        resp.plot(
            tt, vresFin, label='V-Res-Final', marker='o', ls="-.", color='b')

        iterp.plot(
            tt, iters, '-', marker='.', ls='-.', label='NumIters', color='b')

        contp.plot(
            tt,
            contresInit,
            marker='*',
            label='Cont-Res-Initial',
            ls="--",
            color='r')
        contp.plot(
            tt,
            contresFin,
            marker='o',
            label='Cont-Res-Final',
            ls="-.",
            color='r')
        contp.plot(
            tt,
            presresInit,
            marker='*',
            label='Pres-Res-Initial',
            ls="--",
            color='b')
        contp.plot(
            tt,
            presresFin,
            marker='o',
            label='Pres-Res-Final',
            ls="-.",
            color='b')

        if i == 0:
            resp.legend(loc=2)
            contp.legend()
            iterp.legend(loc=1)
            i = 1

        # print "Printing..."
        plt.pause(10)


if __name__ == '__main__':
    #plot("/Users/calebbuahin/Documents/Projects/HydroCouple/FVHMComponent/examples/test.csv")
    plot(sys.argv[1])
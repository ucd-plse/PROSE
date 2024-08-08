import sys
import os
import time
import numpy as np
from numpy import linalg as LA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def get_time(log_path):

    with open(os.path.join(log_path, "outlog.txt"), "r") as f:
        lines = f.readlines()

    total_runtime = 0
    number_clocks = 0
    number_runs = 0

    # gather total runtime
    i = -1
    while i + 1 < len(lines):
        i += 1

        if lines[i].strip().startswith("clock runtime"):
            total_runtime += float(lines[i].split()[2])
            number_clocks += 1

        elif lines[i].strip().startswith("MPI terminated with Status"):
            number_runs += 1

    if total_runtime == 0:
        total_runtime = np.nan
    else:
        total_runtime = number_runs * total_runtime/number_clocks

    return total_runtime


class ADCrun():
    def __init__(self, casepath):
        self._casepath = casepath
        self._read_14()

    def _read_14(self):
        f14 = open(os.path.join(self._casepath,"fort.14"))
        header = f14.readline()
        line2 = f14.readline()
        self.ne = int(line2.split()[0])
        self.nn = int(line2.split()[1])

        # read nodes
        self.x = np.zeros(self.nn)
        self.y = np.zeros(self.nn)
        self.d = np.zeros(self.nn)
        for i in range(self.nn):
            sline = f14.readline().split()
            self.x[i] = float(sline[1])
            self.y[i] = float(sline[2])
            self.d[i] = float(sline[3])

        # read elements (converted to zero-based indices)
        self.elem = [[None,None,None] for i in range(self.ne)]
        for i in range(self.ne):
            sline = f14.readline().split()
            self.elem[i][0] = int(sline[2])-1
            self.elem[i][1] = int(sline[3])-1
            self.elem[i][2] = int(sline[4])-1

    def read_maxele(self, fname="maxele.63"):
        m63 = open(os.path.join(self._casepath,fname))
        line = m63.readline() # header
        sline = m63.readline().split()
        assert int(sline[1]) == self.nn, "incorrect number of nodes in maxele"
        line = m63.readline()
        self.maxele = np.zeros(self.nn)
        for i in range(self.nn):
            sline = m63.readline().split()
            self.maxele[i] = float(sline[1])
            if self.maxele[i] == -99999.0: # dry node correction
                self.maxele[i] = np.NaN

    def read_maxvel(self, fname="maxvel.63"):
        mv63 = open(os.path.join(self._casepath,fname))
        line = mv63.readline() # header
        sline = mv63.readline().split()
        assert int(sline[1]) == self.nn, "incorrect number of nodes in maxele"
        line = mv63.readline()
        self.maxvel = np.zeros(self.nn)
        for i in range(self.nn):
            sline = mv63.readline().split()
            self.maxvel[i] = float(sline[1])
            if self.maxvel[i] == -99999.0: # dry node correction
                self.maxvel[i] = np.NaN

    def diff_norm(self, q1, q2, fillval=None, ord=None):
        if fillval is None:
            fillval = -self.d
        q1_masked = np.ma.masked_invalid(q1)
        q1_masked = q1_masked.filled(fill_value=fillval)
        q2_masked = np.ma.masked_invalid(q2)
        q2_masked = q2_masked.filled(fill_value=fillval)
        return LA.norm(q1_masked-q2_masked, ord)


    def plot_quantity(self, q, filename, xlim=None, ylim=None, qlim=None, cmap=None, nbins=20, title=None, symmetric_q=False):
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri

        # preprocess the quantity by masking the invalid values:
        qprime = np.ma.masked_invalid(q)
        qmin, qmax = qprime.min(), qprime.max()
        qprime = qprime.filled(fill_value=-99999.0)
        if qlim:
            qmin, qmax = qlim[0], qlim[1]
        if symmetric_q:
            qmax = max(abs(qmin), abs(qmax))
            qmin = - qmax
        levels = np.linspace(qmin, qmax, nbins)

        triang = mtri.Triangulation(self.x, self.y, self.elem)

        if not (xlim or ylim):
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16,9))
            try:
                z = ax.tricontourf(triang, qprime,  levels=levels, cmap=cmap)
            except ValueError:
                return
            ax.triplot(triang, color='0.9', alpha=0.3)
        else:
            fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(16,6))
            try:
                z = axs[0].tricontourf(triang, qprime, levels=levels, cmap=cmap)
                z = axs[1].tricontourf(triang, qprime, levels=levels, cmap=cmap)
            except ValueError:
                return
            axs[0].triplot(triang, color='0.99', alpha=0.3)
            axs[1].triplot(triang, color='0.99', alpha=0.4)

            if xlim:
                axs[1].set_xlim(xlim)
            if ylim:
                axs[1].set_ylim(ylim)

            axs[0].set_xlabel('Longitude [$^\circ$E]')
            axs[0].set_ylabel('Latitude [$^\circ$N]')
            axs[1].set_xlabel('Longitude [$^\circ$E]')

        if title:
            fig.suptitle(title)

        fig.tight_layout()
        fig.colorbar(z)
        #plt.show()

        fig.savefig(filename)


if __name__ == "__main__":

    log_path = sys.argv[1]
    log_path_for_baseline = sys.argv[2]
    ELE_L2_NORM_CUTOFF = 1.00e-1
    ELE_MAX_NORM_CUTOFF = 1.00e-1
    VEL_L2_NORM_CUTOFF = 3.00e-1
    VEL_MAX_NORM_CUTOFF = 3.00e-1

    timestr = time.strftime("%Y%m%d-%H%M%S")

    case = ADCrun('./')
    case.read_maxele()
    case.read_maxvel()

    # copy maxele and maxvel files
    import shutil
    shutil.copy("./maxele.63", log_path)
    shutil.copy("./maxvel.63", log_path)

    total_runtime = get_time(log_path)
    if np.isnan(total_runtime):
        total_runtime = 0 # error code recognized by Prose for a failing configuration that terminated gracefully
    elif "prose_logs/0000" not in log_path:
        case0 = ADCrun('./')
        case0.read_maxele(os.path.join(log_path_for_baseline, "maxele.63"))
        case0.read_maxvel(os.path.join(log_path_for_baseline, "maxvel.63"))

        ele_err_l2_norm = case.diff_norm(case0.maxele, case.maxele)
        ele_err_max_norm = case.diff_norm(case0.maxele, case.maxele, ord=np.inf)

        vel_err_l2_norm = case.diff_norm(case0.maxvel, case.maxvel, fillval=0.0)
        vel_err_max_norm = case.diff_norm(case0.maxvel, case.maxvel, fillval=0.0, ord=np.inf)

        feasible = ele_err_l2_norm < ELE_L2_NORM_CUTOFF

        error_metrics_file = open(os.path.join(log_path,"error_metrics.txt"), "a")
        error_metrics_file.write("{}  {:10.4f}  {:2}  {:8.6f}  {:8.6f}  {:8.6f}  {:8.6f}\n".format(log_path[-3:], total_runtime,
            feasible, ele_err_l2_norm, ele_err_max_norm,
            vel_err_l2_norm, vel_err_max_norm))

        # plot maximum elevation error
        case.plot_quantity(case0.maxele - case.maxele, filename=os.path.join(log_path, 'maxele_{}.png'.format(timestr)),
            xlim=(-80,-72), ylim=(32,40), qlim=(-ELE_MAX_NORM_CUTOFF, ELE_MAX_NORM_CUTOFF), cmap='turbo',nbins=40,
            title="Max elevation error (m) - l2-norm: {:e}".format(ele_err_l2_norm), symmetric_q=True)
        case.plot_quantity(case0.maxvel - case.maxvel, filename=os.path.join(log_path, 'maxvel_{}.png'.format(timestr)),
            xlim=(-80,-72), ylim=(32,40), qlim=(-VEL_MAX_NORM_CUTOFF, VEL_MAX_NORM_CUTOFF), cmap='turbo',nbins=40,
            title="Max velocity error (m/s) - l2-norm: {:e}".format(vel_err_l2_norm), symmetric_q=True)

        if (not feasible):
            total_runtime = -1 * abs(total_runtime)

    print(total_runtime)
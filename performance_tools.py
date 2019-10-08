import numpy as np
import time
import os


class PerformanceTimer:

    def __init__(self, allow_profiling = True):
        self.performance = {}
        self.num_calls   = {}
        self._start_time = {}
        self.allow_profiling = allow_profiling

        self.ignore_in_sum = {}

        return

    def timing_on(self):
        self.allow_profiling = True
        return

    def timing_off(self):
        self.allow_profiling = False
        return

    def get_timers(self):
        return self.performance.keys()

    def start_timer(self, timername, subprocess=False):
        if not self.allow_profiling:
            return

        if timername in self._start_time.keys():
            if self._start_time[timername] > 0:
                print("Error. Timer already started: ", self._start_time[timername], timername)
                raise RuntimeError

        self.ignore_in_sum[timername] = subprocess

        self._start_time[timername] = time.process_time()

        if not timername in self.performance.keys():
            self.performance[timername] = 0.0
            self.num_calls[timername]   = 0

        return

    def end_timer(self, timername):
        if not self.allow_profiling:
            return

        if not timername in self._start_time.keys():
            print("Trying to end a performance timer does not exist ", timername)
            raise RuntimeError

        self.performance[timername] += (time.process_time() - self._start_time[timername])
        self.num_calls[timername]   += 1
        self._start_time[timername] = 0.0

        return

    def print_performance(self):

        if not self.allow_profiling:
            return

        total_time = np.sum([ self.performance[k] for k in self.performance.keys() if (not self.ignore_in_sum[k])])

        for k in self.performance.keys():
            n = self.num_calls[k]
            t = self.performance[k]
            if self.ignore_in_sum[k]:
                frac = 0.0
            else:
                frac = t / total_time
            e
            print("%20s  %00005i  %5.5E   %5.5E"%(k,n,t,frac))

        return

    def write_performance(self, outstr = None, outname = 'performance.out'):
        if not self.allow_profiling:
            return

        write_header=False
        if not os.path.isfile(outname):
            write_header = True

        f = open(outname, 'a')

        if write_header:
            f.write("name     count     time    frac\n")

        if not (outstr is None):
            f.write(outstr)

        total_time = np.sum([ self.performance[k] for k in self.performance.keys() if (not self.ignore_in_sum[k])])

        for k in self.performance.keys():
            n = self.num_calls[k]
            t = self.performance[k]
            if self.ignore_in_sum[k]:
                frac = 0.0
            else:
                frac = t / total_time

            f.write("%20s  %00005i  %5.5E   %5.5E\n"%(k,n,t,frac))

        f.close()

        return

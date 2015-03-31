#!/usr/bin/env python
from __future__ import print_function

from multiprocessing import Queue, Process
import subprocess as sp
import shlex
import psutil
import time

class get_runtime_data(object):
    def __init__(self, command, sleeptime=2, refresh_freq=0.1, analyze=[]):
        self._command = command
        self.sleeptime = sleeptime
        self.refresh_freq = refresh_freq
        if len(analyze) > 0:
            self.analyze = analyze
        else:
            self.analyze = ['cpu_user_time', 'cpu_sys_time', 'max_mex', 'real_time']

        self.analfunc = {
            'max_mem': self._max_mem,
            'real_time': self._pass,
            'cpu_user_time': self._cpu_user_time,
            'cpu_sys_time': self._cpu_sys_time,
        }
        self.anres = {}

    def run(self):
        p = Process(self._wrapper())
        p.start()
        print(self.anres)
        return self.anres

    def _wrapper(self):
        resultq = Queue()
        p = Process(target=self._subproc, args=(resultq, self._command,))
        p.start()
        ps = psutil.Process(p.pid)
        time.sleep(self.sleeptime)
        self._analyzing(ps)
        self.anres['real_time'] = resultq.get()
                
    def _subproc(self, q, command):
        ti = time.time()
        proc = sp.Popen(shlex.split(command))
        proc.communicate()
        tf = time.time()
        q.put(tf-ti)

    def _analyzing(self,ps):
        while len(ps.get_children()) > 0:
            time.sleep(self.refresh_freq)
            for an in self.analyze:
                # try:
                    self.anres[an] = self.analfunc[an](ps)
                # except:
                #     pass
            

    def _pass(self,ps):
        pass

    def _cpu_user_time(self,ps):
        return ps.get_children()[0].get_cpu_times()
    
    def _cpu_sys_time(self,ps):
        return ps.get_children()[0].get_cpu_times()[1]
        
    def _get_thread_num(self,ps):
        return ps.get_children()[0].get_num_threads()

    def _max_mem(self,ps):
        mem = ps.get_memory_info()[0]
        try:
            if self.membuffer < mem:
                self.membuffer = mem
        except:
            self.membuffer = mem
        return self.membuffer

if __name__ == '__main__':
    import os

    command = "./op_f.x"
    thread_test = [1,2,4,8]

    for thread in thread_test:
        os.environ['OMP_NUM_THREADS'] = '{:d}'.format(thread)
        asd = get_runtime_data(command)
        asd.run()
        del(asd)


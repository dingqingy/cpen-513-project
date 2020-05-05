#! /usr/bin/env python3

# Copyright (c) 2019, NVIDIA CORPORATION. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of NVIDIA CORPORATION nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import inspect
import os
import pprint
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np
import libconf
import yaml
import json

this_file_path = os.path.abspath(inspect.getfile(inspect.currentframe()))
this_directory = os.path.dirname(this_file_path)
root_dir = os.path.join(os.path.dirname(this_file_path), '..')

sys.path.append(os.path.join(root_dir, 'scripts'))
from cnn_layers import *
import timeloop
import parse_timeloop_output

config_abspath = os.path.join(root_dir, 'configs/mapper/VGGN.yaml')

# Just test that path points to a valid config file.
with open(config_abspath, 'r') as f:
    config = yaml.load(f)
    #config = libconf.load(f)

edp = {}
result = {}
for name in names:
    edp[name] = sys.float_info.max

# for t in [1e24]:  # [1e24, 1e12, 1e6, 1e3]
t = 1e24
for cooling in [10, 100]:
    for max_iter in [1000]:
        for beta in [0.8, 0.9]:
            for name in names:
                problem = layers[name]

                print("Preparing to run timeloop for problem conv ", i)

                dirname = f'run/{name}/{t}_{cooling}_{max_iter}_{beta}/'
                subprocess.check_call(['mkdir', '-p', dirname])

                timeloop.run_timeloop(dirname,
                                      configfile=config_abspath,
                                      workload_bounds=problem,
                                      t=t,
                                      cooling=cooling,
                                      max_iter=max_iter,
                                      beta=beta)

                stats = parse_timeloop_output.parse_timeloop_stats(dirname)
                if stats == {}:
                    print("Timeloop couldn't find a mapping for this problem within the search parameters, please check the log for more details.")
                else:
                    print("Run successful, see log for text stats, or use the Python parser to parse the XML stats.")
                    # print("Stats from run:")
                    # pprint.pprint(stats)

                energy_delay_product = stats['energy_pJ'] * stats['cycles']
                if energy_delay_product < edp[name]:
                    edp[name] = energy_delay_product
                    result[name] = stats

            # print("DONE.")

with open(f'edp_{t}.json', 'w') as f:
    json.dump(edp, f)

with open(f'result_{t}.json', 'w') as f:
    json.dump(result, f)

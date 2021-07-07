"""The main program that runs gSpan."""
# -*- coding=utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys

from .config import parser
from .gspan import gSpan
from .gspan import flat_trans
import codecs
import collections
import copy
import itertools
import time
graph_cnt=0

def read_graphs(FLAGS=None):
    global graph_cnt
    if FLAGS is None:
        FLAGS, _ = parser.parse_known_args(args=sys.argv[1:])

    database_file_name=FLAGS.database_file_name
    with codecs.open(database_file_name, 'r', 'utf-8') as f:
        lines = [line.strip() for line in f.readlines()]
        for i, line in enumerate(lines):
            cols = line.split(' ')
            if cols[0] == 't':
                graph_cnt += 1
    return graph_cnt

def main(FLAGS=None):
    """Run gSpan."""

    if FLAGS is None:
        FLAGS, _ = parser.parse_known_args(args=sys.argv[1:])

    if not os.path.exists(FLAGS.database_file_name):
        print('{} does not exist.'.format(FLAGS.database_file_name))
        sys.exit()

    gs = gSpan(
        database_file_name=FLAGS.database_file_name,
        min_support=FLAGS.min_support,
        min_num_vertices=FLAGS.lower_bound_of_num_vertices,
        max_num_vertices=FLAGS.upper_bound_of_num_vertices,
        max_ngraphs=FLAGS.num_graphs,
        is_undirected=(not FLAGS.directed),
        verbose=FLAGS.verbose,
        visualize=FLAGS.plot,
        where=FLAGS.where
    )
    graph_cnt=read_graphs()
    gs.run()
    gs.time_stats()
    min_sup=sys.argv[2]
    s=sys.argv[5]
    p=sorted(flat_trans)
    # print("aa",graph_cnt)
    #for i in range(0,graph_cnt):
        # print(len(flat_trans[i]))
        #if(len(flat_trans[i])==0):
        #    flat_trans[i].append(-1)
        # print(i,flat_trans[i])
    f=open(str(s)+"_Flat_tra.txt",'w')
    for i in range(0,graph_cnt):
        sarr=[str(a) for a in flat_trans[i]]
        outstr=str(" ".join(sarr))+"\n"
        f.write(outstr)
    f.close()
    return gs
    

if __name__ == '__main__':
    main()

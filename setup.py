#!/usr/bin/env python

from setuptools import setup

setup(
    name='cytospace',
    version='1.0.0',
    packages=['cytospace',
              'cytospace.common',
              'cytospace.linear_assignment_solvers',
              'cytospace.post_processing'],
    package_dir={'cytospace': 'cytospace'},
    entry_points={'console_scripts':
                      ['cytospace=cytospace.cytospace:run_cytospace']#,
                       #'cytospace-plot=cytospace.post_processing.plot:run_cytospace_plot']
                  }
)

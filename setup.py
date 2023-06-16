#!/usr/bin/env python

from setuptools import setup

setup(
    name='cytospace',
    version='1.0.5',
    packages=['cytospace',
              'cytospace.common',
              'cytospace.linear_assignment_solvers',
              'cytospace.post_processing'],
    package_dir={'cytospace': 'cytospace'},
    package_data={'cytospace': ['get_cellfracs_seuratv3.R']},
    entry_points={'console_scripts':
                      ['cytospace=cytospace.cytospace:run_cytospace']
                  },
    include_package_data=True
)

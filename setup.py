#!/usr/bin/env python

# Instructions:
# Replace ?????? with package (i.e. folder name)
#

from setuptools import setup, find_packages
from setuptools.command.develop import develop

class run_post_develop(develop):
    def run(self):
        develop.run(self)
        
        from fsc.plugins import install_plugin
        install_plugin.main()
        

setup(name='fsc',
      version='0.1',
      description='calculates fsc/frc',
      author='Kenny Chung',
      author_email='kenny.chung@yale.edu',
      url='',
      packages=find_packages(),
      # package_data={
      #       # include all svg and html files, otherwise conda will miss them
      #       '': ['*.svg', '*.html'],
      # }
      cmdclass = {
              'develop': run_post_develop,
              },
     )

from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='kinact',
      version='0.1',
      description='Toolbox for kinase activity estimation',
      url='https://github.com/saezlab/kinact',
      author='Jakob Wirbel',
      author_email='jakob.wirbel@gmail.com',
      license='GNU GPLv3 License',
      packages=['kinact'],
      long_description=readme(),
      keywords='Kinase activity',
      zip_safe=False,
      install_requires=[
          'pandas',
          'numpy',
          'scipy',
          'statsmodels'
      ],
      include_package_data=True
      )

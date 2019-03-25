from setuptools import setup

setup(name='aneuploidy_survival',
      version='0.1',
      description='Python scripts to support analysis of whole chromosome and arm-length aneuploidy.',
      url='http://github.com/joansmith/aneuploidy-survival',
      author='Joan Smith',
      author_email='joans@alum.mit.edu',
      license='MIT',
      packages=['aneuploidy_survival'],
      zip_safe=False,
      install_requires=[
        'biomarker_survival'])

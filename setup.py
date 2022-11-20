from setuptools import setup, find_packages

setup(
    name='AIMDprobe',
    version='0.0.1',
    url='https://github.com/tjunewson/AIMDprobe.git',
    author='Liu et al.',
    author_email='newson@tju.edu.cn',
    description='A simple toolkit to analyze ab initio molecular dynamics (AIMD) simulation results',
    packages=find_packages(),    
    install_requires=['numpy >= 1.11.1', 'matplotlib >= 1.5.1'],
)
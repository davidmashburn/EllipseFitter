from distutils.core import setup

# Read the version number
with open("EllipseFitter/_version.py") as f:
    exec(f.read())

setup(
    name='EllipseFitter',
    version=__version__, # use the same version that's in _version.py
    author='David N. Mashburn',
    author_email='david.n.mashburn@gmail.com',
    packages=['EllipseFitter'],
    scripts=[],
    url='http://pypi.python.org/pypi/EllipseFitter/',
    license='LICENSE.txt',
    description='find best fit ellipses of binary images',
    long_description=open('README.rst').read(),
    install_requires=[
                      'numpy>=1.0'
                     ],
)

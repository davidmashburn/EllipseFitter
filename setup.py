from distutils.core import setup

setup(
    name='EllipseFitter',
    version='0.1.2',
    author='David N. Mashburn',
    author_email='david.n.mashburn@gmail.com',
    packages=['EllipseFitter'],
    scripts=[],
    url='http://pypi.python.org/pypi/EllipseFitter/',
    license='LICENSE.txt',
    description='',
    long_description=open('README.rst').read(),
    install_requires=[
                      'numpy>=1.0'
                     ],
)

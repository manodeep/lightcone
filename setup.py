#!/usr/bin/env python

import os
import sys
import subprocess
import codecs


if sys.version_info[:2] < (2, 7) or (3, 0) <= sys.version_info[:2] < (3, 4):
    raise RuntimeError("Python version 2.7 or >= 3.4 required.")

if sys.version_info[0] >= 3:
    import builtins
else:
    import __builtin__ as builtins
            

MAJOR               = 0
MINOR               = 0
MICRO               = 1
ISRELEASED          = False
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        
        for k in ['SYSTEMROOT', 'PATH', 'HOME']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v

        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out
    
    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"
        
    return GIT_REVISION


# Most of the following are directly taken from numpy/setup.py: MS 30/08/2017

# BEFORE importing setuptools, remove MANIFEST. Otherwise it may not be
# properly updated when the contents of directories change (true for distutils,
# not sure about setuptools).
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

# This is a bit hackish: we are setting a global variable so that the main
# __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet. 
builtins.__LIGHTCONE_SETUP__ = True


def get_version_info():
    # Adding the git rev number needs to be done inside write_version_py(),
    # otherwise the import of lightcone.version messes up the build under
    # Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('lightcone/version.py'):
        # must be a source distribution, use existing version file
        try:
            from lightcone.version import git_revision as GIT_REVISION
        except ImportError:
            raise ImportError("Unable to import git_revision. Try removing " \
                              "lightcone/version.py and the build directory " \
                              "before building.")
    else:
        GIT_REVISION = "Unknown"
        
    if not ISRELEASED:
        FULLVERSION += '.dev0+' + GIT_REVISION[:7]
                
    return FULLVERSION, GIT_REVISION


def write_version_py(filename='lightcone/version.py'):
    cnt = """
# THIS FILE IS AUTO-GENERATED FROM LIGHTCONE SETUP.PY
#
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION, GIT_REVISION = get_version_info()
    
    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version' : FULLVERSION,
                       'git_revision' : GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()



curdir = os.path.dirname(os.path.abspath(__file__))
readme_file = os.path.join(curdir, 'README.rst')
LONG_DESCRIPTION = codecs.open(readme_file, 'r', 'utf-8').read()

CLASSIFIERS="""\
Development Status :: 3 - Alpha
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: MIT License
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering
"""

NAME = "lightcone"
PACKAGES = [NAME]
PACKAGE_DATA = {NAME: [os.path.join('data', '*')]}
INSTALL_REQUIRES = ["setuptools"]
RUN_REQUIRES = ["numpy"]

def setup_package():
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)
    
    # Rewrite the version file everytime
    write_version_py()
                            
    metadata = dict(name=NAME,
                    maintainer="Manodeep Sinha",
                    maintainer_email="manodeep@gmail.com",
                    description="lightcone: Create lightcones from simulated "\
                    "galaxies and halos",
                    long_description=LONG_DESCRIPTION,
                    url="https://github.com/manodeep/lightcone",
                    license="MIT",
                    classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
                    author="Manodeep Sinha",
                    author_email="manodeep@gmail.com",
                    platforms="OS Independent",
                    version=VERSION,
                    packages=PACKAGES,
                    package_data=PACKAGE_DATA,
                    install_requires=INSTALL_REQUIRES,
                    requires=RUN_REQUIRES)
    
    from setuptools import setup    
    try:
        setup(**metadata)
    finally:
        del sys.path[0]
        os.chdir(old_path)
        
    return


if __name__ == '__main__':
    setup_package()
    del builtins.__LIGHTCONE_SETUP__

from setuptools import setup, find_packages, find_namespace_packages

# To install in linux, type in terminal "sudo python3 setup.py install"

# For developer install, use "sudo python3 setup.py develop" to install 
# and "sudo python3 setup.py develop -u" to remove.

setup(
    name="mission_tools",
    version="1.0",
    # packages=find_packages(),
    packages=find_namespace_packages(),
    #scripts=["say_hello.py"],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=["Cartopy >= 0.16.0",
                    "matplotlib == 2.2.3", # Newer versions of matplotlib break Cartopy
                    "netCDF4 >= 1.5.3",
                    "numpy >= 1.18.1",
                    "pandas >= 1.0.0",
                    "sgp4 >= 1.4",
                    "spacepy >= 0.1.6"],

    # package_data={
    #     # If any package contains *.txt or *.rst files, include them:
    #     "": ["*.txt", "*.rst"],
    #     # And include any *.msg files found in the "hello" package, too:
    #     "hello": ["*.msg"],
    # },

    # metadata to display on PyPI
    author="Mykhaylo Shumko",
    author_email="msshumko@gmail.com",
    description="Data processing tools for a handful of magnetospheric missions",
    keywords="FIREBIRD, AC6, RBSP, GOES, SGP4",
    url="https://github.com/mshumko/mission_tools",   # project home page, if any
    classifiers=[
        "License :: OSI Approved :: Python Software Foundation License"
    ]
)
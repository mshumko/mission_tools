import setuptools

setuptools.setup(
    name="mission_tools",
    version="1.0",
    packages=setuptools.find_packages(),
    #packages=find_namespace_packages(include=["ac6.*", 'firebird.*']),
    #scripts=["say_hello.py"],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    # install_requires=["Cartopy >= 0.16.0",
    #                 "matplotlib == 3.1.2", # Newer versions of matplotlib break Cartopy
    #                 "netCDF4 >= 1.5.3",
    #                 "numpy >= 1.18.1",
    #                 "pandas >= 1.0.0",
    #                 "sgp4 >= 1.4",
    #                 "spacepy >= 0.1.6"],

    # metadata to display on PyPI
    author="Mykhaylo Shumko",
    author_email="msshumko@gmail.com",
    description="Data processing tools for a handful of magnetospheric missions",
    keywords="FIREBIRD, AC6, RBSP, GOES, SGP4, gemagetic indicies, SAMPEX",
    url="https://github.com/mshumko/mission_tools",   # project home page, if any
)
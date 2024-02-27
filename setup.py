from setuptools import setup, find_packages

setup(
    name="vast_mw",
    version="0.1.0",
    description="Multi-wavelength searches for VAST sources",
    author="David Kaplan",
    author_email="kaplan@uwm.edu",
    url="",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "check_gaia=vast_mw.scripts.check_gaia:main",
            "check_pulsarscraper=vast_mw.scripts.check_pulsarscraper:main",
            "check_simbad=vast_mw.scripts.check_simbad:main",
        ]
    },
    python_requires=">=3.7",
    install_requires=["astropy", "astroquery", "loguru", "numpy"],
    package_data={},
    include_package_data=False,
    zip_safe=False,
)

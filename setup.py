from setuptools import setup, find_packages

setup(
    name="schto",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "PyYAML"
    ],
    entry_points={
        "console_scripts": [
            "schto=src.main:main"
        ]
    }
)
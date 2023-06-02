from setuptools import setup
from Cython.Build import cythonize


# setup(
#     name='formating app',
#     ext_modules=cythonize("formating.pyx"),
#     zip_safe=False,
# )

setup(
    name='formating app',
    ext_modules=cythonize("formating.pyx"),
    zip_safe=False,
)
# thermostat
A toolbox for creating statistics for thermodynamic models.

All properties are calculated with https://github.com/usnistgov/teqp provided by Ian Bell.



# Build (cmake based)

To build the thermostat tool type in console in root if thermostat:

mkdir build
cd build
cmake .. -DTEQP_NO_PYTHON=ON -DTEQP_NO_TESTS=ON
cmake --build . --config Release

(So far only tested on Windows using Visual C++ compiler)

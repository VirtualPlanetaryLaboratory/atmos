CLIMA/SETUP

This subdirectory has the subroutines to start the climate code.

- choose_star.f
Reads the flux from stars that are not the Sun

- grid.f
Sets the pressure grid

- interp.f, interpozone.f, irexpsums.f
They read data used by the old IR code

- profile.f
Constructs a new temperature and water profiles from specific surface and 
stratospheric temperatures.

- readsol.f
Reads all the data needed by the solar subroutine

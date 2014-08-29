#include "brg/global.h"

#include "brg/physics/units/units.h"

#include "position_grid_cache.h"

// Initialisation for brgastro::grid_cache
#if (1)

int brgastro::grid_cache::_ra_grid_change_num_ = 0;
int brgastro::grid_cache::_dec_grid_change_num_ = 0;
int brgastro::grid_cache::_z_grid_change_num_ = 0;
BRG_ANGLE brgastro::grid_cache::_ra_grid_min_ = -pi;
BRG_ANGLE brgastro::grid_cache::_ra_grid_max_ = pi;
BRG_ANGLE brgastro::grid_cache::_ra_grid_step_ = pi / 8;
BRG_ANGLE brgastro::grid_cache::_dec_grid_min_ = -pi / 2;
BRG_ANGLE brgastro::grid_cache::dec_grid_max_val = pi / 2;
BRG_ANGLE brgastro::grid_cache::_dec_grid_step_ = pi / 8;
double brgastro::grid_cache::_z_grid_min_ = 0;
double brgastro::grid_cache::_z_grid_max_ = 2;
double brgastro::grid_cache::_z_grid_step_ = 0.1;

#endif // End initialisation for brgastro::grid_cache


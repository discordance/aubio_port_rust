use ::libc;
extern "C" {
    #[no_mangle]
    fn fabsf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
}
/*
  Copyright (C) 2003-2015 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Definition of data types used in aubio

*/
/* * defined to 1 if aubio is compiled in double precision */
/* * short sample format (32 or 64 bits) */
pub type smpl_t = libc::c_float;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = libc::c_uint;
pub type C2RustUnnamed = libc::c_uint;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_parameter_t {
    pub current_value: smpl_t,
    pub target_value: smpl_t,
    pub increment: smpl_t,
    pub max_value: smpl_t,
    pub min_value: smpl_t,
    pub steps: uint_t,
}
/*
  Copyright (C) 2003-2013 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Parameter with linear interpolation

  This object manages a parameter, with minimum and maximum values, and a
  number of steps to compute linear interpolation between two values.

  \example utils/test-parameter.c

*/
/* * parameter object */
pub type aubio_parameter_t = _aubio_parameter_t;
/* * create new parameter object

  \param min_value the minimum value of the new parameter
  \param max_value the maximum value of the new parameter
  \param steps the number of steps to interpolate from the old value to the target value

  \return the newly created ::aubio_parameter_t

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_parameter(mut min_value: smpl_t,
                                             mut max_value: smpl_t,
                                             mut steps: uint_t)
 -> *mut aubio_parameter_t {
    let mut param: *mut aubio_parameter_t =
        calloc(::std::mem::size_of::<aubio_parameter_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_parameter_t;
    (*param).min_value = min_value;
    (*param).max_value = max_value;
    (*param).steps = steps;
    (*param).current_value = (*param).min_value;
    (*param).target_value = (*param).current_value;
    (*param).increment = 0.0f64 as smpl_t;
    return param;
}
/* * set target value of the parameter

  \param param parameter, created by ::new_aubio_parameter
  \param value new target value

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_parameter_set_target_value(mut param:
                                                              *mut aubio_parameter_t,
                                                          mut value: smpl_t)
 -> uint_t {
    let mut err: uint_t = AUBIO_OK as libc::c_int as uint_t;
    if value < (*param).min_value {
        (*param).target_value = (*param).min_value;
        err = AUBIO_FAIL as libc::c_int as uint_t
    } else if value > (*param).max_value {
        (*param).target_value = (*param).max_value;
        err = AUBIO_FAIL as libc::c_int as uint_t
    } else { (*param).target_value = value }
    (*param).increment =
        ((*param).target_value - (*param).current_value) /
            (*param).steps as libc::c_float;
    return err;
}
/* * set current parameter value, skipping interpolation

  \param param parameter, created by ::new_aubio_parameter
  \param value new parameter value

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_parameter_set_current_value(mut param:
                                                               *mut aubio_parameter_t,
                                                           mut value: smpl_t)
 -> uint_t {
    let mut err: uint_t = AUBIO_OK as libc::c_int as uint_t;
    if value < (*param).min_value {
        (*param).current_value = (*param).min_value;
        err = AUBIO_FAIL as libc::c_int as uint_t
    } else if value > (*param).max_value {
        (*param).current_value = (*param).max_value;
        err = AUBIO_FAIL as libc::c_int as uint_t
    } else { (*param).current_value = value }
    (*param).target_value = (*param).current_value;
    (*param).increment = 0 as libc::c_int as smpl_t;
    return err;
}
/* * get current parameter value, without interpolation

  \param param parameter, created by ::new_aubio_parameter

  \return current value

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_parameter_get_current_value(mut s:
                                                               *const aubio_parameter_t)
 -> smpl_t {
    return (*s).current_value;
}
/* * get next parameter

  \param param parameter, created by ::new_aubio_parameter

  \return new interpolated parameter value

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_parameter_get_next_value(mut s:
                                                            *mut aubio_parameter_t)
 -> smpl_t {
    if fabsf((*s).current_value - (*s).target_value) > fabsf((*s).increment) {
        (*s).current_value += (*s).increment
    } else { (*s).current_value = (*s).target_value }
    return (*s).current_value;
}
/* * set number of steps used for interpolation

  \param param parameter, created by ::new_aubio_parameter
  \param steps new number of steps

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_parameter_set_steps(mut param:
                                                       *mut aubio_parameter_t,
                                                   mut steps: uint_t)
 -> uint_t {
    if steps < 1 as libc::c_int as libc::c_uint ||
           steps > 2000 as libc::c_int as libc::c_uint {
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    (*param).steps = steps;
    (*param).increment =
        ((*param).target_value - (*param).current_value) /
            (*param).steps as libc::c_float;
    return AUBIO_OK as libc::c_int as uint_t;
}
/* * get number of steps of this parameter

  \param param parameter, created by ::new_aubio_parameter

  \return number of steps

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_parameter_get_steps(mut param:
                                                       *const aubio_parameter_t)
 -> uint_t {
    return (*param).steps;
}
/* * set minimum value of this parameter

  \param param parameter, created by ::new_aubio_parameter
  \param min_value new minimum value

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_parameter_set_min_value(mut param:
                                                           *mut aubio_parameter_t,
                                                       mut min_value: smpl_t)
 -> uint_t {
    (*param).min_value = min_value;
    return AUBIO_OK as libc::c_int as uint_t;
}
/* * get minimum value of this parameter

  \param param parameter, created by ::new_aubio_parameter

  \return minimum value

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_parameter_get_min_value(mut param:
                                                           *const aubio_parameter_t)
 -> smpl_t {
    return (*param).min_value;
}
/* * set maximum value of this parameter

  \param param parameter, created by ::new_aubio_parameter
  \param max_value new maximum value

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_parameter_set_max_value(mut param:
                                                           *mut aubio_parameter_t,
                                                       mut max_value: smpl_t)
 -> uint_t {
    (*param).max_value = max_value;
    return AUBIO_OK as libc::c_int as uint_t;
}
/* * get maximum value of this parameter

  \param param parameter, created by ::new_aubio_parameter

  \return maximum value

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_parameter_get_max_value(mut param:
                                                           *const aubio_parameter_t)
 -> smpl_t {
    return (*param).max_value;
}
/* * destroy ::aubio_parameter_t object

  \param param parameter, created by ::new_aubio_parameter

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_parameter(mut param:
                                                 *mut aubio_parameter_t) {
    free(param as *mut libc::c_void);
}

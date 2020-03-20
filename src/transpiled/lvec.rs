extern "C" {
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    /* Logging */
    /* * internal logging function, defined in utils/log.c */
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
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
pub type smpl_t = f32;
/* * print format for sample in single precision */
/* * long sample format (64 bits or more) */
pub type lsmp_t = f64;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = u32;
/* * signed integer */
pub type sint_t = i32;
/* * character */
pub type char_t = i8;
/*
  Copyright (C) 2016 Paul Brossier <piem@aubio.org>

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

  Logging features

  This file specifies ::aubio_log_set_function and
  ::aubio_log_set_level_function, which let you define one or several custom
  logging functions to redirect warnings and errors from aubio to your
  application. The custom function should have the prototype defined in
  ::aubio_log_function_t.

  After a call to ::aubio_log_set_level_function, ::aubio_log_reset can be used
  to reset each logging functions to the default ones.

  \example utils/test-log.c

*/
/* * list of logging levels */
pub type aubio_log_level = u32;
/* *< number of valid levels */
/* *< warnings */
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
/* *< debug messages */
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
/* *< general messages */
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
/* *< infos */
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
/* *< critical errors */
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct lvec_t {
    pub length: uint_t,
    pub data: *mut lsmp_t,
}
/* * lvec_t buffer creation function

  \param length the length of the buffer to create

*/
/*
  Copyright (C) 2003-2009 Paul Brossier <piem@aubio.org>

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
#[no_mangle]
pub unsafe extern "C" fn new_lvec(length: uint_t) -> *mut lvec_t {
    let mut s: *mut lvec_t = 0 as *mut lvec_t;
    if length as sint_t <= 0 as i32 { return 0 as *mut lvec_t }
    s =
        calloc(::std::mem::size_of::<lvec_t>() as u64,
               1 as i32 as u64) as *mut lvec_t;
    (*s).length = length;
    (*s).data =
        calloc(((*s).length as
                    u64).wrapping_mul(::std::mem::size_of::<lsmp_t>()
                                                    as u64),
               1 as i32 as u64) as *mut lsmp_t;
    return s;
}
/* * lvec_t buffer deletion function

  \param s buffer to delete as returned by new_lvec()

*/
#[no_mangle]
pub unsafe extern "C" fn del_lvec(s: *mut lvec_t) {
    free((*s).data as *mut core::ffi::c_void);
    free(s as *mut core::ffi::c_void);
}
/* * write sample value in a buffer

  \param s vector to write to
  \param data value to write in s->data[position]
  \param position sample position to write to

*/
#[no_mangle]
pub unsafe extern "C" fn lvec_set_sample(s: *mut lvec_t, data: lsmp_t,
                                         position: uint_t) {
    *(*s).data.offset(position as isize) = data;
}
/* * read sample value in a buffer

  \param s vector to read from
  \param position sample position to read from

*/
#[no_mangle]
pub unsafe extern "C" fn lvec_get_sample(s: *mut lvec_t,
                                         position: uint_t) -> lsmp_t {
    return *(*s).data.offset(position as isize);
}
/* * read data from a buffer

  \param s vector to read from

*/
#[no_mangle]
pub unsafe extern "C" fn lvec_get_data(s: *const lvec_t) -> *mut lsmp_t {
    return (*s).data;
}
/* * print out lvec data

  \param s vector to print out

*/
/* helper functions */
#[no_mangle]
pub unsafe extern "C" fn lvec_print(s: *const lvec_t) {
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        aubio_log(AUBIO_LOG_MSG as i32,
                  b"%lf \x00" as *const u8 as *const i8,
                  *(*s).data.offset(j as isize));
        j = j.wrapping_add(1)
    }
    aubio_log(AUBIO_LOG_MSG as i32,
              b"\n\x00" as *const u8 as *const i8);
}
/* * set all elements to a given value

  \param s vector to modify
  \param val value to set elements to

*/
#[no_mangle]
pub unsafe extern "C" fn lvec_set_all(s: *mut lvec_t, val: smpl_t) {
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        *(*s).data.offset(j as isize) = val as lsmp_t;
        j = j.wrapping_add(1)
    };
}
/* * set all elements to zero

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn lvec_zeros(s: *mut lvec_t) {
    lvec_set_all(s, 0.0f64 as smpl_t);
}
/* * set all elements to ones

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn lvec_ones(s: *mut lvec_t) {
    lvec_set_all(s, 1.0f64 as smpl_t);
}

extern "C" {
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    #[no_mangle]
    fn logf(_: f32) -> f32;
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
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = u32;
/* * signed integer */
pub type sint_t = i32;
/* * character */
pub type char_t = i8;
pub type aubio_log_level = u32;
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct cvec_t {
    pub length: uint_t,
    pub norm: *mut smpl_t,
    pub phas: *mut smpl_t,
}
/* * cvec_t buffer creation function

  This function creates a cvec_t structure holding two arrays of size
  [length/2+1], corresponding to the norm and phase values of the
  spectral frame. The length stored in the structure is the actual size of both
  arrays, not the length of the complex and symmetrical vector, specified as
  creation argument.

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
pub unsafe extern "C" fn new_cvec(length: uint_t) -> *mut cvec_t {
    let mut s: *mut cvec_t = 0 as *mut cvec_t;
    if length as sint_t <= 0 as i32 { return 0 as *mut cvec_t }
    s =
        calloc(::std::mem::size_of::<cvec_t>() as u64,
               1 as i32 as u64) as *mut cvec_t;
    (*s).length =
        length.wrapping_div(2 as i32 as
                                u32).wrapping_add(1 as i32 as
                                                               u32);
    (*s).norm =
        calloc(((*s).length as
                    u64).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                    as u64),
               1 as i32 as u64) as *mut smpl_t;
    (*s).phas =
        calloc(((*s).length as
                    u64).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                    as u64),
               1 as i32 as u64) as *mut smpl_t;
    return s;
}
/* * cvec_t buffer deletion function

  \param s buffer to delete as returned by new_cvec()

*/
#[no_mangle]
pub unsafe extern "C" fn del_cvec(s: *mut cvec_t) {
    free((*s).norm as *mut core::ffi::c_void);
    free((*s).phas as *mut core::ffi::c_void);
    free(s as *mut core::ffi::c_void);
}
/* * write norm value in a complex buffer

  This is equivalent to:
  \code
  s->norm[position] = val;
  \endcode

  \param s vector to write to
  \param val norm value to write in s->norm[position]
  \param position sample position to write to

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_norm_set_sample(s: *mut cvec_t,
                                              data: smpl_t,
                                              position: uint_t) {
    *(*s).norm.offset(position as isize) = data;
}
/* * write phase value in a complex buffer

  This is equivalent to:
  \code
  s->phas[position] = val;
  \endcode

  \param s vector to write to
  \param val phase value to write in s->phas[position]
  \param position sample position to write to

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_phas_set_sample(s: *mut cvec_t,
                                              data: smpl_t,
                                              position: uint_t) {
    *(*s).phas.offset(position as isize) = data;
}
/* * read norm value from a complex buffer

  This is equivalent to:
  \code
  smpl_t foo = s->norm[position];
  \endcode

  \param s vector to read from
  \param position sample position to read from

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_norm_get_sample(s: *mut cvec_t,
                                              position: uint_t)
 -> smpl_t {
    return *(*s).norm.offset(position as isize);
}
/* * read phase value from a complex buffer

  This is equivalent to:
  \code
  smpl_t foo = s->phas[position];
  \endcode

  \param s vector to read from
  \param position sample position to read from
  \returns the value of the sample at position

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_phas_get_sample(s: *mut cvec_t,
                                              position: uint_t)
 -> smpl_t {
    return *(*s).phas.offset(position as isize);
}
/* * read norm data from a complex buffer

  \code
  smpl_t *data = s->norm;
  \endcode

  \param s vector to read from

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_norm_get_data(s: *const cvec_t)
 -> *mut smpl_t {
    return (*s).norm;
}
/* * read phase data from a complex buffer

  This is equivalent to:
  \code
  smpl_t *data = s->phas;
  \endcode

  \param s vector to read from

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_phas_get_data(s: *const cvec_t)
 -> *mut smpl_t {
    return (*s).phas;
}
/* * print out cvec data

  \param s vector to print out

*/
/* helper functions */
#[no_mangle]
pub unsafe extern "C" fn cvec_print(s: *const cvec_t) {
    let mut j: uint_t = 0;
    aubio_log(AUBIO_LOG_MSG as i32,
              b"norm: \x00" as *const u8 as *const i8);
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        aubio_log(AUBIO_LOG_MSG as i32,
                  b"%f \x00" as *const u8 as *const i8,
                  *(*s).norm.offset(j as isize) as f64);
        j = j.wrapping_add(1)
    }
    aubio_log(AUBIO_LOG_MSG as i32,
              b"\n\x00" as *const u8 as *const i8);
    aubio_log(AUBIO_LOG_MSG as i32,
              b"phas: \x00" as *const u8 as *const i8);
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        aubio_log(AUBIO_LOG_MSG as i32,
                  b"%f \x00" as *const u8 as *const i8,
                  *(*s).phas.offset(j as isize) as f64);
        j = j.wrapping_add(1)
    }
    aubio_log(AUBIO_LOG_MSG as i32,
              b"\n\x00" as *const u8 as *const i8);
}
/* * make a copy of a vector

  \param s source vector
  \param t vector to copy to

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_copy(s: *const cvec_t, t: *mut cvec_t) {
    if (*s).length != (*t).length {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: trying to copy %d elements to %d elements \n\x00"
                      as *const u8 as *const i8, (*s).length,
                  (*t).length);
        return
    }
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*t).length {
        *(*t).norm.offset(j as isize) = *(*s).norm.offset(j as isize);
        *(*t).phas.offset(j as isize) = *(*s).phas.offset(j as isize);
        j = j.wrapping_add(1)
    };
}
/* * set all norm elements to a given value

  \param s vector to modify
  \param val value to set elements to

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_norm_set_all(s: *mut cvec_t,
                                           val: smpl_t) {
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        *(*s).norm.offset(j as isize) = val;
        j = j.wrapping_add(1)
    };
}
/* * set all norm elements to zero

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_norm_zeros(s: *mut cvec_t) {
    cvec_norm_set_all(s, 0.0f64 as smpl_t);
}
/* * set all norm elements to one

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_norm_ones(s: *mut cvec_t) {
    cvec_norm_set_all(s, 1.0f64 as smpl_t);
}
/* * set all phase elements to a given value

  \param s vector to modify
  \param val value to set elements to

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_phas_set_all(s: *mut cvec_t,
                                           val: smpl_t) {
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        *(*s).phas.offset(j as isize) = val;
        j = j.wrapping_add(1)
    };
}
/* * set all phase elements to zero

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_phas_zeros(s: *mut cvec_t) {
    cvec_phas_set_all(s, 0.0f64 as smpl_t);
}
/* * set all phase elements to one

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_phas_ones(s: *mut cvec_t) {
    cvec_phas_set_all(s, 1.0f64 as smpl_t);
}
/* * set all norm and phas elements to zero

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_zeros(s: *mut cvec_t) {
    cvec_norm_zeros(s);
    cvec_phas_zeros(s);
}
/* * take logarithmic magnitude

  \param s input cvec to compress
  \param lambda value to use for normalisation

  \f$ S_k = log( \lambda * S_k + 1 ) \f$

*/
#[no_mangle]
pub unsafe extern "C" fn cvec_logmag(s: *mut cvec_t, lambda: smpl_t) {
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        *(*s).norm.offset(j as isize) =
            logf(lambda * *(*s).norm.offset(j as isize) +
                     1 as i32 as f32);
        j = j.wrapping_add(1)
    };
}

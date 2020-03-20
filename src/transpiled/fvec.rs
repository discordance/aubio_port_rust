extern "C" {
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    #[no_mangle]
    fn floorf(_: f32) -> f32;
    /* file interface */
    /* strings */
    /* Error reporting */
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
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = u32;
/* * signed integer */
pub type sint_t = i32;
/* * character */
pub type char_t = i8;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct fvec_t {
    pub length: uint_t,
    pub data: *mut smpl_t,
}
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
pub type aubio_log_level = u32;
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
pub const AUBIO_LOG_INF: aubio_log_level = 1;
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
pub unsafe extern "C" fn new_fvec(length: uint_t) -> *mut fvec_t {
    let mut s: *mut fvec_t = 0 as *mut fvec_t;
    if length as sint_t <= 0 as i32 { return 0 as *mut fvec_t }
    s =
        calloc(::std::mem::size_of::<fvec_t>() as u64,
               1 as i32 as u64) as *mut fvec_t;
    (*s).length = length;
    (*s).data =
        calloc(((*s).length as
                    u64).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                    as u64),
               1 as i32 as u64) as *mut smpl_t;
    return s;
}
#[no_mangle]
pub unsafe extern "C" fn del_fvec(s: *mut fvec_t) {
    free((*s).data as *mut core::ffi::c_void);
    free(s as *mut core::ffi::c_void);
}
#[no_mangle]
pub unsafe extern "C" fn fvec_set_sample(s: *mut fvec_t, data: smpl_t,
                                         position: uint_t) {
    *(*s).data.offset(position as isize) = data;
}
#[no_mangle]
pub unsafe extern "C" fn fvec_get_sample(s: *const fvec_t,
                                         position: uint_t) -> smpl_t {
    return *(*s).data.offset(position as isize);
}
#[no_mangle]
pub unsafe extern "C" fn fvec_get_data(s: *const fvec_t) -> *mut smpl_t {
    return (*s).data;
}
/* helper functions */
#[no_mangle]
pub unsafe extern "C" fn fvec_print(s: *const fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        aubio_log(AUBIO_LOG_MSG as i32,
                  b"%f \x00" as *const u8 as *const i8,
                  *(*s).data.offset(j as isize) as f64);
        j = j.wrapping_add(1)
    }
    aubio_log(AUBIO_LOG_MSG as i32,
              b"\n\x00" as *const u8 as *const i8);
}
#[no_mangle]
pub unsafe extern "C" fn fvec_set_all(s: *mut fvec_t, val: smpl_t) {
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        *(*s).data.offset(j as isize) = val;
        j = j.wrapping_add(1)
    };
}
#[no_mangle]
pub unsafe extern "C" fn fvec_zeros(s: *mut fvec_t) {
    fvec_set_all(s, 0.0f64 as smpl_t);
}
#[no_mangle]
pub unsafe extern "C" fn fvec_ones(s: *mut fvec_t) {
    fvec_set_all(s, 1.0f64 as smpl_t);
}
#[no_mangle]
pub unsafe extern "C" fn fvec_rev(s: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while (j as f32) <
              floorf((*s).length as smpl_t /
                         2 as i32 as f32) {
        let t: smpl_t = *(*s).data.offset(j as isize);
        *(*s).data.offset(j as isize) =
            *(*s).data.offset((*s).length.wrapping_sub(1 as i32 as
                                                           u32).wrapping_sub(j)
                                  as isize);
        *(*s).data.offset((*s).length.wrapping_sub(1 as i32 as
                                                       u32).wrapping_sub(j)
                              as isize) = t;
        j = j.wrapping_add(1)
    };
}
#[no_mangle]
pub unsafe extern "C" fn fvec_weight(s: *mut fvec_t,
                                     weight: *const fvec_t) {
    let length: uint_t =
        if (*s).length < (*weight).length {
            (*s).length
        } else { (*weight).length };
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < length {
        let ref mut fresh0 = *(*s).data.offset(j as isize);
        *fresh0 *= *(*weight).data.offset(j as isize);
        j = j.wrapping_add(1)
    };
    /* HAVE_ACCELERATE */
}
/* * make a copy of a vector, applying weights to each element

  \param in input vector
  \param weight weights vector
  \param out output vector

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_weighted_copy(in_0: *const fvec_t,
                                            weight: *const fvec_t,
                                            out: *mut fvec_t) {
    let length: uint_t =
        if (*in_0).length <
               (if (*out).length < (*weight).length {
                    (*out).length
                } else { (*weight).length }) {
            (*in_0).length
        } else if (*out).length < (*weight).length {
            (*out).length
        } else { (*weight).length };
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < length {
        *(*out).data.offset(j as isize) =
            *(*in_0).data.offset(j as isize) *
                *(*weight).data.offset(j as isize);
        j = j.wrapping_add(1)
    };
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

  Vector of real-valued data

  This file specifies the ::fvec_t buffer type, which is used throughout aubio
  to store vector of real-valued ::smpl_t.

  \example test-fvec.c

*/
/* * Buffer for real data

  Vector of real-valued data

  ::fvec_t is is the structure used to store vector of real-valued data, ::smpl_t .

  \code

  uint_t buffer_size = 1024;

  // create a vector of 512 values
  fvec_t * input = new_fvec (buffer_size);

  // set some values of the vector
  input->data[23] = 2.;
  // ..

  // compute the mean of the vector
  mean = fvec_mean(a_vector);

  // destroy the vector
  del_fvec(a_vector);

  \endcode

  See `examples/` and `tests/src` directories for more examples.

 */
/* *< length of buffer */
/* *< data vector of length ::fvec_t.length */
/* * fvec_t buffer creation function

  \param length the length of the buffer to create

*/
/* * fvec_t buffer deletion function

  \param s buffer to delete as returned by new_fvec()

*/
/* * read sample value in a buffer

  \param s vector to read from
  \param position sample position to read from

*/
/* * write sample value in a buffer

  \param s vector to write to
  \param data value to write in s->data[position]
  \param position sample position to write to

*/
/* * read data from a buffer

  \param s vector to read from

*/
/* * print out fvec data

  \param s vector to print out

*/
/* * set all elements to a given value

  \param s vector to modify
  \param val value to set elements to

*/
/* * set all elements to zero

  \param s vector to modify

*/
/* * set all elements to ones

  \param s vector to modify

*/
/* * revert order of vector elements

  \param s vector to revert

*/
/* * apply weight to vector

  If the weight vector is longer than s, only the first elements are used. If
  the weight vector is shorter than s, the last elements of s are not weighted.

  \param s vector to weight
  \param weight weighting coefficients

*/
/* * make a copy of a vector

  \param s source vector
  \param t vector to copy to

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_copy(s: *const fvec_t, t: *mut fvec_t) {
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
        *(*t).data.offset(j as isize) = *(*s).data.offset(j as isize);
        j = j.wrapping_add(1)
    };
}

use ::libc;
extern "C" {
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn floorf(_: libc::c_float) -> libc::c_float;
    /* * set all elements to zero

  \param s vector to modify

*/
    #[no_mangle]
    fn fvec_zeros(s: *mut fvec_t);
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
pub type smpl_t = libc::c_float;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = libc::c_uint;
/* * signed integer */
pub type sint_t = libc::c_int;
/* * character */
pub type char_t = libc::c_char;
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
#[derive(Copy, Clone)]
#[repr(C)]
pub struct fvec_t {
    pub length: uint_t,
    pub data: *mut smpl_t,
}
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
pub type aubio_log_level = libc::c_uint;
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
pub struct fmat_t {
    pub length: uint_t,
    pub height: uint_t,
    pub data: *mut *mut smpl_t,
}
/* * fmat_t buffer creation function

  \param length the length of the matrix to create
  \param height the height of the matrix to create

*/
/*
  Copyright (C) 2009 Paul Brossier <piem@aubio.org>

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
pub unsafe extern "C" fn new_fmat(mut height: uint_t, mut length: uint_t)
 -> *mut fmat_t {
    let mut s: *mut fmat_t = 0 as *mut fmat_t;
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    if length as sint_t <= 0 as libc::c_int ||
           height as sint_t <= 0 as libc::c_int {
        return 0 as *mut fmat_t
    }
    s =
        calloc(::std::mem::size_of::<fmat_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut fmat_t;
    (*s).height = height;
    (*s).length = length;
    (*s).data =
        calloc(((*s).height as
                    libc::c_ulong).wrapping_mul(::std::mem::size_of::<*mut smpl_t>()
                                                    as libc::c_ulong),
               1 as libc::c_int as libc::c_ulong) as *mut *mut smpl_t;
    i = 0 as libc::c_int as uint_t;
    while i < (*s).height {
        let ref mut fresh0 = *(*s).data.offset(i as isize);
        *fresh0 =
            calloc(((*s).length as
                        libc::c_ulong).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                        as libc::c_ulong),
                   1 as libc::c_int as libc::c_ulong) as *mut smpl_t;
        j = 0 as libc::c_int as uint_t;
        while j < (*s).length {
            *(*(*s).data.offset(i as isize)).offset(j as isize) =
                0.0f64 as smpl_t;
            j = j.wrapping_add(1)
        }
        i = i.wrapping_add(1)
    }
    return s;
}
/* * fmat_t buffer deletion function

  \param s buffer to delete as returned by new_fmat()

*/
#[no_mangle]
pub unsafe extern "C" fn del_fmat(mut s: *mut fmat_t) {
    let mut i: uint_t = 0;
    i = 0 as libc::c_int as uint_t;
    while i < (*s).height {
        free(*(*s).data.offset(i as isize) as *mut libc::c_void);
        i = i.wrapping_add(1)
    }
    free((*s).data as *mut libc::c_void);
    free(s as *mut libc::c_void);
}
/* * write sample value in a buffer

  \param s vector to write to
  \param data value to write in s->data[channel][position]
  \param channel channel to write to
  \param position sample position to write to

*/
#[no_mangle]
pub unsafe extern "C" fn fmat_set_sample(mut s: *mut fmat_t, mut data: smpl_t,
                                         mut channel: uint_t,
                                         mut position: uint_t) {
    *(*(*s).data.offset(channel as isize)).offset(position as isize) = data;
}
/* * read sample value in a buffer

  \param s vector to read from
  \param channel channel to read from
  \param position sample position to read from

*/
#[no_mangle]
pub unsafe extern "C" fn fmat_get_sample(mut s: *const fmat_t,
                                         mut channel: uint_t,
                                         mut position: uint_t) -> smpl_t {
    return *(*(*s).data.offset(channel as isize)).offset(position as isize);
}
/* * read channel vector from a buffer

  \param s vector to read from
  \param channel channel to read from
  \param output ::fvec_t to output to

*/
#[no_mangle]
pub unsafe extern "C" fn fmat_get_channel(mut s: *const fmat_t,
                                          mut channel: uint_t,
                                          mut output: *mut fvec_t) {
    (*output).data = *(*s).data.offset(channel as isize);
    (*output).length = (*s).length;
}
/* * get vector buffer from an fmat data

  \param s vector to read from
  \param channel channel to read from

*/
#[no_mangle]
pub unsafe extern "C" fn fmat_get_channel_data(mut s: *const fmat_t,
                                               mut channel: uint_t)
 -> *mut smpl_t {
    return *(*s).data.offset(channel as isize);
}
/* * read data from a buffer

  \param s vector to read from

*/
#[no_mangle]
pub unsafe extern "C" fn fmat_get_data(mut s: *const fmat_t)
 -> *mut *mut smpl_t {
    return (*s).data;
}
/* * print out fmat data

  \param s vector to print out

*/
/* helper functions */
#[no_mangle]
pub unsafe extern "C" fn fmat_print(mut s: *const fmat_t) {
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    i = 0 as libc::c_int as uint_t;
    while i < (*s).height {
        j = 0 as libc::c_int as uint_t;
        while j < (*s).length {
            aubio_log(AUBIO_LOG_MSG as libc::c_int,
                      b"%f \x00" as *const u8 as *const libc::c_char,
                      *(*(*s).data.offset(i as isize)).offset(j as isize) as
                          libc::c_double);
            j = j.wrapping_add(1)
        }
        aubio_log(AUBIO_LOG_MSG as libc::c_int,
                  b"\n\x00" as *const u8 as *const libc::c_char);
        i = i.wrapping_add(1)
    };
}
/* * set all elements to a given value

  \param s vector to modify
  \param val value to set elements to

*/
#[no_mangle]
pub unsafe extern "C" fn fmat_set(mut s: *mut fmat_t, mut val: smpl_t) {
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    i = 0 as libc::c_int as uint_t;
    while i < (*s).height {
        j = 0 as libc::c_int as uint_t;
        while j < (*s).length {
            *(*(*s).data.offset(i as isize)).offset(j as isize) = val;
            j = j.wrapping_add(1)
        }
        i = i.wrapping_add(1)
    };
}
/* * set all elements to zero

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn fmat_zeros(mut s: *mut fmat_t) {
    /* HAVE_MEMCPY_HACKS */
    fmat_set(s, 0.0f64 as smpl_t);
    /* HAVE_MEMCPY_HACKS */
}
/* * set all elements to ones

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn fmat_ones(mut s: *mut fmat_t) {
    fmat_set(s, 1.0f64 as smpl_t);
}
/* * revert order of vector elements

  \param s vector to revert

*/
#[no_mangle]
pub unsafe extern "C" fn fmat_rev(mut s: *mut fmat_t) {
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    i = 0 as libc::c_int as uint_t;
    while i < (*s).height {
        j = 0 as libc::c_int as uint_t;
        while (j as libc::c_float) <
                  floorf((*s).length as smpl_t /
                             2 as libc::c_int as libc::c_float) {
            let mut t: smpl_t =
                *(*(*s).data.offset(i as isize)).offset(j as isize);
            *(*(*s).data.offset(i as isize)).offset(j as isize) =
                *(*(*s).data.offset(i as
                                        isize)).offset((*s).length.wrapping_sub(1
                                                                                    as
                                                                                    libc::c_int
                                                                                    as
                                                                                    libc::c_uint).wrapping_sub(j)
                                                           as isize);
            *(*(*s).data.offset(i as
                                    isize)).offset((*s).length.wrapping_sub(1
                                                                                as
                                                                                libc::c_int
                                                                                as
                                                                                libc::c_uint).wrapping_sub(j)
                                                       as isize) = t;
            j = j.wrapping_add(1)
        }
        i = i.wrapping_add(1)
    };
}
/* * apply weight to vector

  If the weight vector is longer than s, only the first elements are used. If
  the weight vector is shorter than s, the last elements of s are not weighted.

  \param s vector to weight
  \param weight weighting coefficients

*/
#[no_mangle]
pub unsafe extern "C" fn fmat_weight(mut s: *mut fmat_t,
                                     mut weight: *const fmat_t) {
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    let mut length: uint_t =
        if (*s).length < (*weight).length {
            (*s).length
        } else { (*weight).length };
    i = 0 as libc::c_int as uint_t;
    while i < (*s).height {
        j = 0 as libc::c_int as uint_t;
        while j < length {
            let ref mut fresh1 =
                *(*(*s).data.offset(i as isize)).offset(j as isize);
            *fresh1 *=
                *(*(*weight).data.offset(0 as libc::c_int as
                                             isize)).offset(j as isize);
            j = j.wrapping_add(1)
        }
        i = i.wrapping_add(1)
    };
}
/* * make a copy of a matrix

  \param s source vector
  \param t vector to copy to

*/
#[no_mangle]
pub unsafe extern "C" fn fmat_copy(mut s: *const fmat_t, mut t: *mut fmat_t) {
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    /* HAVE_MEMCPY_HACKS */
    if (*s).height != (*t).height {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: trying to copy %d rows to %d rows \n\x00" as
                      *const u8 as *const libc::c_char, (*s).height,
                  (*t).height);
        return
    }
    if (*s).length != (*t).length {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: trying to copy %d columns to %d columns\n\x00"
                      as *const u8 as *const libc::c_char, (*s).length,
                  (*t).length);
        return
    }
    /* HAVE_MEMCPY_HACKS */
    i = 0 as libc::c_int as uint_t;
    while i < (*t).height {
        j = 0 as libc::c_int as uint_t;
        while j < (*t).length {
            *(*(*t).data.offset(i as isize)).offset(j as isize) =
                *(*(*s).data.offset(i as isize)).offset(j as isize);
            j = j.wrapping_add(1)
        }
        i = i.wrapping_add(1)
    };
    /* HAVE_MEMCPY_HACKS */
}
/* * compute the product of a matrix by a vector

   \param s matrix to compute product with
   \param scale vector to compute product with
   \param output vector to store restults in

*/
#[no_mangle]
pub unsafe extern "C" fn fmat_vecmul(mut s: *const fmat_t,
                                     mut scale: *const fvec_t,
                                     mut output: *mut fvec_t) {
    let mut k: uint_t = 0;
    let mut j: uint_t = 0;
    fvec_zeros(output);
    j = 0 as libc::c_int as uint_t;
    while j < (*s).length {
        k = 0 as libc::c_int as uint_t;
        while k < (*s).height {
            let ref mut fresh2 = *(*output).data.offset(k as isize);
            *fresh2 +=
                *(*scale).data.offset(j as isize) *
                    *(*(*s).data.offset(k as isize)).offset(j as isize);
            k = k.wrapping_add(1)
        }
        j = j.wrapping_add(1)
    };
}

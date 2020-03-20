use ::libc;
extern "C" {
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    #[no_mangle]
    fn logf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn sqrtf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn floorf(_: libc::c_float) -> libc::c_float;
    /* * fvec_t buffer creation function

  \param length the length of the buffer to create

*/
    #[no_mangle]
    fn new_fvec(length: uint_t) -> *mut fvec_t;
    /* * fvec_t buffer deletion function

  \param s buffer to delete as returned by new_fvec()

*/
    #[no_mangle]
    fn del_fvec(s: *mut fvec_t);
    /* * make a copy of a vector

  \param s source vector
  \param t vector to copy to

*/
    #[no_mangle]
    fn fvec_copy(s: *const fvec_t, t: *mut fvec_t);
    /* * return 1 if a is a power of 2, 0 otherwise */
    #[no_mangle]
    fn aubio_is_power_of_two(a: uint_t) -> uint_t;
    #[no_mangle]
    fn aubio_ooura_ddct(_: libc::c_int, _: libc::c_int, _: *mut smpl_t,
                        _: *mut libc::c_int, _: *mut smpl_t);
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
pub type aubio_log_level = libc::c_uint;
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_dct_ooura_t {
    pub size: uint_t,
    pub input: *mut fvec_t,
    pub w: *mut smpl_t,
    pub ip: *mut libc::c_int,
    pub scalers: [smpl_t; 5],
}
/*
  Copyright (C) 2017 Paul Brossier <piem@aubio.org>

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
pub type aubio_dct_ooura_t = _aubio_dct_ooura_t;
#[no_mangle]
pub unsafe extern "C" fn new_aubio_dct_ooura(mut size: uint_t)
 -> *mut aubio_dct_ooura_t {
    let mut s: *mut aubio_dct_ooura_t =
        calloc(::std::mem::size_of::<aubio_dct_ooura_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_dct_ooura_t;
    if aubio_is_power_of_two(size) != 1 as libc::c_int as libc::c_uint ||
           size as sint_t <= 0 as libc::c_int {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: dct_ooura: can only create with sizes power of two, requested %d\n\x00"
                      as *const u8 as *const libc::c_char, size);
        free(s as *mut libc::c_void);
        return 0 as *mut aubio_dct_ooura_t
    } else {
        (*s).size = size;
        (*s).input = new_fvec((*s).size);
        (*s).w =
            calloc(((*s).size.wrapping_mul(5 as libc::c_int as
                                               libc::c_uint).wrapping_div(4 as
                                                                              libc::c_int
                                                                              as
                                                                              libc::c_uint)
                        as
                        libc::c_ulong).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                        as libc::c_ulong),
                   1 as libc::c_int as libc::c_ulong) as *mut smpl_t;
        (*s).ip =
            calloc(((3 as libc::c_int +
                         ((1 as libc::c_int) <<
                              floorf(logf((*s).size.wrapping_div(2 as
                                                                     libc::c_int
                                                                     as
                                                                     libc::c_uint)
                                              as libc::c_float) /
                                         logf(2 as libc::c_int as
                                                  libc::c_float)) as
                                  libc::c_int) / 2 as libc::c_int) as
                        libc::c_ulong).wrapping_mul(::std::mem::size_of::<libc::c_int>()
                                                        as libc::c_ulong),
                   1 as libc::c_int as libc::c_ulong) as *mut libc::c_int;
        *(*s).ip.offset(0 as libc::c_int as isize) = 0 as libc::c_int;
        (*s).scalers[0 as libc::c_int as usize] =
            (2.0f64 *
                 sqrtf((1.0f64 / (4.0f64 * (*s).size as libc::c_double)) as
                           libc::c_float) as libc::c_double) as smpl_t;
        (*s).scalers[1 as libc::c_int as usize] =
            (2.0f64 *
                 sqrtf((1.0f64 / (2.0f64 * (*s).size as libc::c_double)) as
                           libc::c_float) as libc::c_double) as smpl_t;
        (*s).scalers[2 as libc::c_int as usize] =
            (1.0f64 /
                 (*s).scalers[0 as libc::c_int as usize] as libc::c_double) as
                smpl_t;
        (*s).scalers[3 as libc::c_int as usize] =
            (1.0f64 /
                 (*s).scalers[1 as libc::c_int as usize] as libc::c_double) as
                smpl_t;
        (*s).scalers[4 as libc::c_int as usize] =
            (2.0f64 / (*s).size as libc::c_double) as smpl_t;
        return s
    };
}
#[no_mangle]
pub unsafe extern "C" fn del_aubio_dct_ooura(mut s: *mut aubio_dct_ooura_t) {
    del_fvec((*s).input);
    free((*s).ip as *mut libc::c_void);
    free((*s).w as *mut libc::c_void);
    free(s as *mut libc::c_void);
}
#[no_mangle]
pub unsafe extern "C" fn aubio_dct_ooura_do(mut s: *mut aubio_dct_ooura_t,
                                            mut input: *const fvec_t,
                                            mut output: *mut fvec_t) {
    let mut i: uint_t = 0 as libc::c_int as uint_t;
    fvec_copy(input, (*s).input);
    aubio_ooura_ddct((*s).size as libc::c_int, -(1 as libc::c_int),
                     (*(*s).input).data, (*s).ip, (*s).w);
    // apply orthonormal scaling
    let ref mut fresh0 =
        *(*(*s).input).data.offset(0 as libc::c_int as isize);
    *fresh0 *= (*s).scalers[0 as libc::c_int as usize];
    i = 1 as libc::c_int as uint_t;
    while i < (*(*s).input).length {
        let ref mut fresh1 = *(*(*s).input).data.offset(i as isize);
        *fresh1 *= (*s).scalers[1 as libc::c_int as usize];
        i = i.wrapping_add(1)
    }
    fvec_copy((*s).input, output);
}
#[no_mangle]
pub unsafe extern "C" fn aubio_dct_ooura_rdo(mut s: *mut aubio_dct_ooura_t,
                                             mut input: *const fvec_t,
                                             mut output: *mut fvec_t) {
    let mut i: uint_t = 0 as libc::c_int as uint_t;
    fvec_copy(input, (*s).input);
    let ref mut fresh2 =
        *(*(*s).input).data.offset(0 as libc::c_int as isize);
    *fresh2 *= (*s).scalers[2 as libc::c_int as usize];
    i = 1 as libc::c_int as uint_t;
    while i < (*(*s).input).length {
        let ref mut fresh3 = *(*(*s).input).data.offset(i as isize);
        *fresh3 *= (*s).scalers[3 as libc::c_int as usize];
        i = i.wrapping_add(1)
    }
    let ref mut fresh4 =
        *(*(*s).input).data.offset(0 as libc::c_int as isize);
    *fresh4 = (*fresh4 as libc::c_double * 0.5f64) as smpl_t;
    aubio_ooura_ddct((*s).size as libc::c_int, 1 as libc::c_int,
                     (*(*s).input).data, (*s).ip, (*s).w);
    i = 0 as libc::c_int as uint_t;
    while i < (*(*s).input).length {
        let ref mut fresh5 = *(*(*s).input).data.offset(i as isize);
        *fresh5 *= (*s).scalers[4 as libc::c_int as usize];
        i = i.wrapping_add(1)
    }
    fvec_copy((*s).input, output);
}
// !defined(HAVE_ACCELERATE) && !defined(HAVE_FFTW3)

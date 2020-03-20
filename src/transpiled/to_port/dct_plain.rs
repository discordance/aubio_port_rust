use ::libc;
extern "C" {
    #[no_mangle]
    fn cosf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn sqrtf(_: libc::c_float) -> libc::c_float;
    /* * return 1 if a is a power of 2, 0 otherwise */
    #[no_mangle]
    fn aubio_is_power_of_two(a: uint_t) -> uint_t;
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    /* * fmat_t buffer creation function

  \param length the length of the matrix to create
  \param height the height of the matrix to create

*/
    #[no_mangle]
    fn new_fmat(height: uint_t, length: uint_t) -> *mut fmat_t;
    /* * fmat_t buffer deletion function

  \param s buffer to delete as returned by new_fmat()

*/
    #[no_mangle]
    fn del_fmat(s: *mut fmat_t);
    /* * compute the product of a matrix by a vector

   \param s matrix to compute product with
   \param scale vector to compute product with
   \param output vector to store restults in

*/
    #[no_mangle]
    fn fmat_vecmul(s: *const fmat_t, scale: *const fvec_t,
                   output: *mut fvec_t);
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
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_dct_plain_t {
    pub size: uint_t,
    pub dct_coeffs: *mut fmat_t,
    pub idct_coeffs: *mut fmat_t,
}
/*
  Copyright (C) 2018 Paul Brossier <piem@aubio.org>

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
pub type aubio_dct_plain_t = _aubio_dct_plain_t;
#[no_mangle]
pub unsafe extern "C" fn new_aubio_dct_plain(mut size: uint_t)
 -> *mut aubio_dct_plain_t {
    let mut s: *mut aubio_dct_plain_t =
        calloc(::std::mem::size_of::<aubio_dct_plain_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_dct_plain_t;
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    let mut scaling: smpl_t = 0.;
    if aubio_is_power_of_two(size) == 1 as libc::c_int as libc::c_uint &&
           size > 16 as libc::c_int as libc::c_uint {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: dct_plain: using plain dct but size %d is a power of two\n\x00"
                      as *const u8 as *const libc::c_char, size);
    }
    if size as sint_t <= 0 as libc::c_int {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: dct_plain: can only create with size > 0, requested %d\n\x00"
                      as *const u8 as *const libc::c_char, size);
        del_aubio_dct_plain(s);
        return 0 as *mut aubio_dct_plain_t
    } else {
        (*s).size = size;
        (*s).dct_coeffs = new_fmat(size, size);
        (*s).idct_coeffs = new_fmat(size, size);
        /* compute DCT type-II transformation matrix
     dct_coeffs[j][i] = cos ( j * (i+.5) * PI / n_filters )
  */
        scaling = sqrtf((2.0f64 / size as libc::c_double) as libc::c_float);
        i = 0 as libc::c_int as uint_t;
        while i < size {
            j = 1 as libc::c_int as uint_t;
            while j < size {
                *(*(*(*s).dct_coeffs).data.offset(j as
                                                      isize)).offset(i as
                                                                         isize)
                    =
                    scaling *
                        cosf((j as libc::c_double *
                                  (i as libc::c_double + 0.5f64) *
                                  3.14159265358979323846264338327950288f64 /
                                  size as libc::c_double) as libc::c_float);
                j = j.wrapping_add(1)
            }
            *(*(*(*s).dct_coeffs).data.offset(0 as libc::c_int as
                                                  isize)).offset(i as isize) =
                (1.0f64 / sqrtf(size as libc::c_float) as libc::c_double) as
                    smpl_t;
            i = i.wrapping_add(1)
        }
        /* compute DCT type-III transformation matrix
     idct_coeffs[j][i] = cos ( i * (j+.5) * PI / n_filters )
  */
        scaling = sqrtf((2.0f64 / size as libc::c_double) as libc::c_float);
        j = 0 as libc::c_int as uint_t;
        while j < size {
            i = 1 as libc::c_int as uint_t;
            while i < size {
                *(*(*(*s).idct_coeffs).data.offset(j as
                                                       isize)).offset(i as
                                                                          isize)
                    =
                    scaling *
                        cosf((i as libc::c_double *
                                  (j as libc::c_double + 0.5f64) *
                                  3.14159265358979323846264338327950288f64 /
                                  size as libc::c_double) as libc::c_float);
                i = i.wrapping_add(1)
            }
            *(*(*(*s).idct_coeffs).data.offset(j as
                                                   isize)).offset(0 as
                                                                      libc::c_int
                                                                      as
                                                                      isize) =
                (1.0f64 / sqrtf(size as libc::c_float) as libc::c_double) as
                    smpl_t;
            j = j.wrapping_add(1)
        }
        return s
    };
}
/* * DCT type III orthonormal transform, size * size */
#[no_mangle]
pub unsafe extern "C" fn del_aubio_dct_plain(mut s: *mut aubio_dct_plain_t) {
    if !(*s).dct_coeffs.is_null() { del_fmat((*s).dct_coeffs); }
    if !(*s).idct_coeffs.is_null() { del_fmat((*s).idct_coeffs); }
    free(s as *mut libc::c_void);
}
#[no_mangle]
pub unsafe extern "C" fn aubio_dct_plain_do(mut s: *mut aubio_dct_plain_t,
                                            mut input: *const fvec_t,
                                            mut output: *mut fvec_t) {
    if (*input).length != (*output).length || (*input).length != (*s).size {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: dct_plain: using input length %d, but output length = %d and size = %d\n\x00"
                      as *const u8 as *const libc::c_char, (*input).length,
                  (*output).length, (*s).size);
    }
    fmat_vecmul((*s).dct_coeffs, input, output);
}
#[no_mangle]
pub unsafe extern "C" fn aubio_dct_plain_rdo(mut s: *mut aubio_dct_plain_t,
                                             mut input: *const fvec_t,
                                             mut output: *mut fvec_t) {
    if (*input).length != (*output).length || (*input).length != (*s).size {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: dct_plain: using input length %d, but output length = %d and size = %d\n\x00"
                      as *const u8 as *const libc::c_char, (*input).length,
                  (*output).length, (*s).size);
    }
    fmat_vecmul((*s).idct_coeffs, input, output);
}

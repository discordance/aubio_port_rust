use ::libc;
extern "C" {
    #[no_mangle]
    fn cosf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn sinf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn expf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn logf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn log10f(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn fabsf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn powf(_: libc::c_float, _: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn sqrtf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn ceilf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn floorf(_: libc::c_float) -> libc::c_float;
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
  Copyright (C) 2009-2015 Paul Brossier <piem@aubio.org>

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

  Utility functions for ::fvec_t

 */
/* * compute \f$e^x\f$ of each vector elements

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_exp(mut o: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < (*o).length {
        *(*o).data.offset(j as isize) = expf(*(*o).data.offset(j as isize));
        j = j.wrapping_add(1)
    };
}
/* * compute \f$cos(x)\f$ of each vector elements

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_cos(mut o: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < (*o).length {
        *(*o).data.offset(j as isize) = cosf(*(*o).data.offset(j as isize));
        j = j.wrapping_add(1)
    };
}
/* * compute \f$sin(x)\f$ of each vector elements

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_sin(mut o: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < (*o).length {
        *(*o).data.offset(j as isize) = sinf(*(*o).data.offset(j as isize));
        j = j.wrapping_add(1)
    };
}
/* * compute the \f$abs(x)\f$ of each vector elements

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_abs(mut o: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < (*o).length {
        *(*o).data.offset(j as isize) = fabsf(*(*o).data.offset(j as isize));
        j = j.wrapping_add(1)
    };
}
/* * compute the \f$sqrt(x)\f$ of each vector elements

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_sqrt(mut o: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < (*o).length {
        *(*o).data.offset(j as isize) = sqrtf(*(*o).data.offset(j as isize));
        j = j.wrapping_add(1)
    };
}
/* * compute the \f$log10(x)\f$ of each vector elements

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_log10(mut o: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < (*o).length {
        *(*o).data.offset(j as isize) =
            log10f(if (*(*o).data.offset(j as isize) as libc::c_double) <
                          2.0e-42f64 {
                       2.0e-42f64
                   } else { *(*o).data.offset(j as isize) as libc::c_double }
                       as libc::c_float);
        j = j.wrapping_add(1)
    };
}
/* * compute the \f$log(x)\f$ of each vector elements

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_log(mut o: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < (*o).length {
        *(*o).data.offset(j as isize) =
            logf(if (*(*o).data.offset(j as isize) as libc::c_double) <
                        2.0e-42f64 {
                     2.0e-42f64
                 } else { *(*o).data.offset(j as isize) as libc::c_double } as
                     libc::c_float);
        j = j.wrapping_add(1)
    };
}
/* * compute the \f$floor(x)\f$ of each vector elements

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_floor(mut o: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < (*o).length {
        *(*o).data.offset(j as isize) = floorf(*(*o).data.offset(j as isize));
        j = j.wrapping_add(1)
    };
}
/* * compute the \f$ceil(x)\f$ of each vector elements

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_ceil(mut o: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < (*o).length {
        *(*o).data.offset(j as isize) = ceilf(*(*o).data.offset(j as isize));
        j = j.wrapping_add(1)
    };
}
/* * compute the \f$round(x)\f$ of each vector elements

  \param s vector to modify

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_round(mut o: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < (*o).length {
        *(*o).data.offset(j as isize) =
            floorf((*(*o).data.offset(j as isize) as libc::c_double + 0.5f64)
                       as libc::c_float);
        j = j.wrapping_add(1)
    };
}
/* * raise each vector elements to the power pow

  \param s vector to modify
  \param pow power to raise to

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_pow(mut s: *mut fvec_t, mut power: smpl_t) {
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < (*s).length {
        *(*s).data.offset(j as isize) =
            powf(*(*s).data.offset(j as isize), power);
        j = j.wrapping_add(1)
    };
}

use ::libc;
extern "C" {
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    /* Logging */
    /* * internal logging function, defined in utils/log.c */
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
    /* * make a copy of a matrix

  \param s source vector
  \param t vector to copy to

*/
    #[no_mangle]
    fn fmat_copy(s: *const fmat_t, t: *mut fmat_t);
    /* * compute the product of a matrix by a vector

   \param s matrix to compute product with
   \param scale vector to compute product with
   \param output vector to store restults in

*/
    #[no_mangle]
    fn fmat_vecmul(s: *const fmat_t, scale: *const fvec_t,
                   output: *mut fvec_t);
    /* * raise each vector elements to the power pow

  \param s vector to modify
  \param pow power to raise to

*/
    #[no_mangle]
    fn fvec_pow(s: *mut fvec_t, pow: smpl_t);
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
/* file interface */
/* strings */
/* Error reporting */
pub type C2RustUnnamed = libc::c_uint;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
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
pub struct cvec_t {
    pub length: uint_t,
    pub norm: *mut smpl_t,
    pub phas: *mut smpl_t,
}
/*
  Copyright (C) 2007-2009 Paul Brossier <piem@aubio.org>
                      and Amaury Hazan <ahazan@iua.upf.edu>

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
/* * \brief A structure to store a set of n_filters filters of lenghts win_s */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_filterbank_t {
    pub win_s: uint_t,
    pub n_filters: uint_t,
    pub filters: *mut fmat_t,
    pub norm: smpl_t,
    pub power: smpl_t,
}
/*
  Copyright (C) 2007-2013 Paul Brossier <piem@aubio.org>
                      and Amaury Hazan <ahazan@iua.upf.edu>

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

  Filterbank object

  General-purpose spectral filterbank object.

  \example spectral/test-filterbank.c

*/
/* * filterbank object

  This object stores a matrix of spectral filter coefficients.

 */
pub type aubio_filterbank_t = _aubio_filterbank_t;
/* * create filterbank object

  \param n_filters number of filters to create
  \param win_s size of analysis buffer (and length the FFT transform)

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_filterbank(mut n_filters: uint_t,
                                              mut win_s: uint_t)
 -> *mut aubio_filterbank_t {
    /* allocate space for filterbank object */
    let mut fb: *mut aubio_filterbank_t =
        calloc(::std::mem::size_of::<aubio_filterbank_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_filterbank_t;
    if n_filters as sint_t <= 0 as libc::c_int {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: filterbank: n_filters should be > 0, got %d\n\x00"
                      as *const u8 as *const libc::c_char, n_filters);
    } else if win_s as sint_t <= 0 as libc::c_int {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: filterbank: win_s should be > 0, got %d\n\x00"
                      as *const u8 as *const libc::c_char, win_s);
    } else {
        (*fb).win_s = win_s;
        (*fb).n_filters = n_filters;
        /* allocate filter tables, a matrix of length win_s and of height n_filters */
        (*fb).filters =
            new_fmat(n_filters,
                     win_s.wrapping_div(2 as libc::c_int as
                                            libc::c_uint).wrapping_add(1 as
                                                                           libc::c_int
                                                                           as
                                                                           libc::c_uint));
        (*fb).norm = 1 as libc::c_int as smpl_t;
        (*fb).power = 1 as libc::c_int as smpl_t;
        return fb
    }
    free(fb as *mut libc::c_void);
    return 0 as *mut aubio_filterbank_t;
}
/* * destroy filterbank object

  \param f filterbank object, as returned by new_aubio_filterbank()

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_filterbank(mut fb:
                                                  *mut aubio_filterbank_t) {
    del_fmat((*fb).filters);
    free(fb as *mut libc::c_void);
}
/* * compute filterbank

  \param f filterbank object, as returned by new_aubio_filterbank()
  \param in input spectrum containing an input spectrum of length `win_s`
  \param out output vector containing the energy found in each band, `nfilt` output values

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_filterbank_do(mut f: *mut aubio_filterbank_t,
                                             mut in_0: *const cvec_t,
                                             mut out: *mut fvec_t) {
    /* apply filter to all input channel, provided out has enough channels */
  //uint_t max_filters = MIN (f->n_filters, out->length);
  //uint_t max_length = MIN (in->length, f->filters->length);
    // view cvec->norm as fvec->data
    let mut tmp: fvec_t = fvec_t{length: 0, data: 0 as *mut smpl_t,};
    tmp.length = (*in_0).length;
    tmp.data = (*in_0).norm;
    if (*f).power as libc::c_double != 1.0f64 {
        fvec_pow(&mut tmp, (*f).power);
    }
    fmat_vecmul((*f).filters, &mut tmp, out);
}
/* * return a pointer to the matrix object containing all filter coefficients

  \param f filterbank object, as returned by new_aubio_filterbank()

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_filterbank_get_coeffs(mut f:
                                                         *const aubio_filterbank_t)
 -> *mut fmat_t {
    return (*f).filters;
}
/* * copy filter coefficients to the filterbank

  \param f filterbank object, as returned by new_aubio_filterbank()
  \param filters filter bank coefficients to copy from

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_filterbank_set_coeffs(mut f:
                                                         *mut aubio_filterbank_t,
                                                     mut filter_coeffs:
                                                         *const fmat_t)
 -> uint_t {
    fmat_copy(filter_coeffs, (*f).filters);
    return 0 as libc::c_int as uint_t;
}
/* * set norm parameter

  \param f filterbank object, as returned by new_aubio_filterbank()
  \param norm `1` to norm the filters, `0` otherwise.

  If set to `0`, the filters will not be normalized. If set to `1`,
  each filter will be normalized to one. Defaults to `1`.

  This function should be called *before* setting the filters with one of
  aubio_filterbank_set_triangle_bands(), aubio_filterbank_set_mel_coeffs(),
  aubio_filterbank_set_mel_coeffs_htk(), or
  aubio_filterbank_set_mel_coeffs_slaney().

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_filterbank_set_norm(mut f:
                                                       *mut aubio_filterbank_t,
                                                   mut norm: smpl_t)
 -> uint_t {
    if norm != 0 as libc::c_int as libc::c_float &&
           norm != 1 as libc::c_int as libc::c_float {
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    (*f).norm = norm;
    return AUBIO_OK as libc::c_int as uint_t;
}
/* * get norm parameter

  \param f filterbank object, as returned by new_aubio_filterbank()
  \returns `1` if norm is set, `0` otherwise. Defaults to `1`.

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_filterbank_get_norm(mut f:
                                                       *mut aubio_filterbank_t)
 -> smpl_t {
    return (*f).norm;
}
/* * set power parameter

  \param f filterbank object, as returned by new_aubio_filterbank()
  \param power Raise norm of the input spectrum norm to this power before
  computing filterbank.  Defaults to `1`.

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_filterbank_set_power(mut f:
                                                        *mut aubio_filterbank_t,
                                                    mut power: smpl_t)
 -> uint_t {
    (*f).power = power;
    return AUBIO_OK as libc::c_int as uint_t;
}
/* * get power parameter

  \param f filterbank object, as returned by new_aubio_filterbank()
  \return current power parameter. Defaults to `1`.

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_filterbank_get_power(mut f:
                                                        *mut aubio_filterbank_t)
 -> smpl_t {
    return (*f).power;
}

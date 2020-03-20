use crate::transpiled::filter::aubio_filter_t;

extern "C" {
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
    /* * @file
 * Private include file
 *
 * This file is for inclusion from _within_ the library only.
 */
    /* ********************
 *
 * External includes
 *
 */
    /* must be included before fftw3.h */
    // for CHAR_BIT, in C99 standard
    // --enable-blas=true
    /* HAVE_ACCELERATE */
    /* HAVE_BLAS */
    /* ***
 *
 * SYSTEM INTERFACE
 *
 */
    /* Memory management */
    /* file interface */
    /* strings */
    /* Error reporting */
    /* Logging */
    /* * internal logging function, defined in utils/log.c */
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    /* * returns a pointer to feedback coefficients \f$ a_i \f$

  \param f filter object to get parameters from

  \return a pointer to the \f$ a_0 ... a_i ... a_P \f$ coefficients

*/
    #[no_mangle]
    fn aubio_filter_get_feedback(f: *const aubio_filter_t) -> *mut lvec_t;
    /* * returns a pointer to feedforward coefficients \f$ b_i \f$

  \param f filter object to get coefficients from

  \return a pointer to the \f$ b_0 ... b_i ... b_P \f$ coefficients

*/
    #[no_mangle]
    fn aubio_filter_get_feedforward(f: *const aubio_filter_t) -> *mut lvec_t;
    /* * get order of the filter

  \param f filter to get order from

  \return the order of the filter

*/
    #[no_mangle]
    fn aubio_filter_get_order(f: *const aubio_filter_t) -> uint_t;
    /* * create new filter object

  This function creates a new ::aubio_filter_t object, given the order of the
  filter.

  \param order order of the filter (number of coefficients)

  \return the newly created filter object

*/
    #[no_mangle]
    fn new_aubio_filter(order: uint_t) -> *mut aubio_filter_t;
}
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
pub type C2RustUnnamed = u32;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
pub type aubio_log_level = u32;
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct lvec_t {
    pub length: uint_t,
    pub data: *mut lsmp_t,
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

  Digital filter

  This object stores a digital filter of order \f$n\f$.
  It contains the following data:
    - \f$ n*1 b_i \f$ feedforward coefficients
    - \f$ n*1 a_i \f$ feedback coefficients
    - \f$ n*c x_i \f$ input signal
    - \f$ n*c y_i \f$ output signal

  For convenience, the samplerate of the input signal is also stored in the
  object.

  Feedforward and feedback parameters can be modified using
  aubio_filter_get_feedback() and aubio_filter_get_feedforward().

  The function aubio_filter_do_outplace() computes the following output signal
  \f$ y[n] \f$ from the input signal \f$ x[n] \f$:

  \f{eqnarray*}{
     y[n] = b_0 x[n] & + & b_1 x[n-1] + b_2 x[n-2] + ... + b_P x[n-P] \\
                     & - & a_1 y[n-1] - a_2 y[n-2] - ... - a_P y[n-P] \\
  \f}

  The function aubio_filter_do() executes the same computation but modifies
  directly the input signal (in-place).

  The function aubio_filter_do_filtfilt() version runs the filter twice, first
  forward then backward, to compensate with the phase shifting of the forward
  operation.

  Some convenience functions are provided:
    - new_aubio_filter_a_weighting() and aubio_filter_set_a_weighting(),
    - new_aubio_filter_c_weighting() and aubio_filter_set_c_weighting().
    - new_aubio_filter_biquad() and aubio_filter_set_biquad().

  \example temporal/test-filter.c

*/
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

  Second order Infinite Impulse Response filter

  This file implements a normalised biquad filter (second order IIR):

  \f$ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] - a_1 y[n-1] - a_2 y[n-2] \f$

  The filtfilt version runs the filter twice, forward and backward, to
  compensate the phase shifting of the forward operation.

  See also <a href="http://en.wikipedia.org/wiki/Digital_biquad_filter">Digital
  biquad filter</a> on wikipedia.

  \example temporal/test-biquad.c

*/
/* * set coefficients of a biquad filter

  \param f filter object as returned by new_aubio_filter()
  \param b0 forward filter coefficient
  \param b1 forward filter coefficient
  \param b2 forward filter coefficient
  \param a1 feedback filter coefficient
  \param a2 feedback filter coefficient

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
pub unsafe extern "C" fn aubio_filter_set_biquad(mut f: *mut aubio_filter_t,
                                                 mut b0: lsmp_t,
                                                 mut b1: lsmp_t,
                                                 mut b2: lsmp_t,
                                                 mut a1: lsmp_t,
                                                 mut a2: lsmp_t) -> uint_t {
    let mut order: uint_t = aubio_filter_get_order(f);
    let mut bs: *mut lvec_t = aubio_filter_get_feedforward(f);
    let mut as_0: *mut lvec_t = aubio_filter_get_feedback(f);
    if order != 3 as i32 as u32 {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: order of biquad filter must be 3, not %d\n\x00"
                      as *const u8 as *const i8, order);
        return AUBIO_FAIL as i32 as uint_t
    }
    *(*bs).data.offset(0 as i32 as isize) = b0;
    *(*bs).data.offset(1 as i32 as isize) = b1;
    *(*bs).data.offset(2 as i32 as isize) = b2;
    *(*as_0).data.offset(0 as i32 as isize) = 1.0f64;
    *(*as_0).data.offset(1 as i32 as isize) = a1;
    *(*as_0).data.offset(2 as i32 as isize) = a2;
    return AUBIO_OK as i32 as uint_t;
}
/* * create biquad filter with `b0`, `b1`, `b2`, `a1`, `a2` coeffs

  \param b0 forward filter coefficient
  \param b1 forward filter coefficient
  \param b2 forward filter coefficient
  \param a1 feedback filter coefficient
  \param a2 feedback filter coefficient

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_filter_biquad(mut b0: lsmp_t,
                                                 mut b1: lsmp_t,
                                                 mut b2: lsmp_t,
                                                 mut a1: lsmp_t,
                                                 mut a2: lsmp_t)
 -> *mut aubio_filter_t {
    let mut f: *mut aubio_filter_t =
        new_aubio_filter(3 as i32 as uint_t);
    aubio_filter_set_biquad(f, b0, b1, b2, a1, a2);
    return f;
}

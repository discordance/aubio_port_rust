extern "C" {
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    #[no_mangle]
    fn fabsf(_: f32) -> f32;
    /* * make a copy of a vector

  \param s source vector
  \param t vector to copy to

*/
    #[no_mangle]
    fn fvec_copy(s: *const fvec_t, t: *mut fvec_t);
    /* * lvec_t buffer creation function

  \param length the length of the buffer to create

*/
    #[no_mangle]
    fn new_lvec(length: uint_t) -> *mut lvec_t;
    /* * lvec_t buffer deletion function

  \param s buffer to delete as returned by new_lvec()

*/
    #[no_mangle]
    fn del_lvec(s: *mut lvec_t);
    /* * set all elements to zero

  \param s vector to modify

*/
    #[no_mangle]
    fn lvec_zeros(s: *mut lvec_t);
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
pub type C2RustUnnamed = u32;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct lvec_t {
    pub length: uint_t,
    pub data: *mut lsmp_t,
}
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
/* Requires lsmp_t to be long or double. float will NOT give reliable 
 * results */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_filter_t {
    pub order: uint_t,
    pub samplerate: uint_t,
    pub a: *mut lvec_t,
    pub b: *mut lvec_t,
    pub y: *mut lvec_t,
    pub x: *mut lvec_t,
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
/* * Digital filter

*/
pub type aubio_filter_t = _aubio_filter_t;
/* * filter input vector (out-of-place)

  \param f filter object as returned by new_aubio_filter()
  \param in input vector to filter
  \param out output vector to store filtered input

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_filter_do_outplace(f: *mut aubio_filter_t,
                                                  in_0: *const fvec_t,
                                                  out: *mut fvec_t) {
    fvec_copy(in_0, out);
    aubio_filter_do(f, out);
}
/* * filter input vector (in-place)

  \param f filter object as returned by new_aubio_filter()
  \param in input vector to filter

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_filter_do(f: *mut aubio_filter_t,
                                         in_0: *mut fvec_t) {
    let mut j: uint_t = 0;
    let mut l: uint_t = 0;
    let order: uint_t = (*f).order;
    let x: *mut lsmp_t = (*(*f).x).data;
    let y: *mut lsmp_t = (*(*f).y).data;
    let a: *mut lsmp_t = (*(*f).a).data;
    let b: *mut lsmp_t = (*(*f).b).data;
    j = 0 as i32 as uint_t;
    while j < (*in_0).length {
        /* new input */
        *x.offset(0 as i32 as isize) =
            if (fabsf(*(*in_0).data.offset(j as isize)) as f64) <
                   2.0e-42f64 {
                0.0f64
            } else { *(*in_0).data.offset(j as isize) as f64 };
        *y.offset(0 as i32 as isize) =
            *b.offset(0 as i32 as isize) *
                *x.offset(0 as i32 as isize);
        l = 1 as i32 as uint_t;
        while l < order {
            let ref mut fresh0 = *y.offset(0 as i32 as isize);
            *fresh0 += *b.offset(l as isize) * *x.offset(l as isize);
            let ref mut fresh1 = *y.offset(0 as i32 as isize);
            *fresh1 -= *a.offset(l as isize) * *y.offset(l as isize);
            l = l.wrapping_add(1)
        }
        /* new output */
        *(*in_0).data.offset(j as isize) =
            *y.offset(0 as i32 as isize) as smpl_t;
        /* store for next sample */
        l = order.wrapping_sub(1 as i32 as u32);
        while l > 0 as i32 as u32 {
            *x.offset(l as isize) =
                *x.offset(l.wrapping_sub(1 as i32 as u32) as
                              isize);
            *y.offset(l as isize) =
                *y.offset(l.wrapping_sub(1 as i32 as u32) as
                              isize);
            l = l.wrapping_sub(1)
        }
        j = j.wrapping_add(1)
    };
}
/* * filter input vector forward and backward

  \param f ::aubio_filter_t object as returned by new_aubio_filter()
  \param in ::fvec_t input vector to filter
  \param tmp memory space to use for computation

*/
/* The rough way: reset memory of filter between each run to avoid end effects. */
#[no_mangle]
pub unsafe extern "C" fn aubio_filter_do_filtfilt(f: *mut aubio_filter_t,
                                                  in_0: *mut fvec_t,
                                                  tmp: *mut fvec_t) {
    let mut j: uint_t = 0;
    let length: uint_t = (*in_0).length;
    /* apply filtering */
    aubio_filter_do(f, in_0);
    aubio_filter_do_reset(f);
    /* mirror */
    j = 0 as i32 as uint_t;
    while j < length {
        *(*tmp).data.offset(length.wrapping_sub(j).wrapping_sub(1 as
                                                                    i32
                                                                    as
                                                                    u32)
                                as isize) = *(*in_0).data.offset(j as isize);
        j = j.wrapping_add(1)
    }
    /* apply filtering on mirrored */
    aubio_filter_do(f, tmp);
    aubio_filter_do_reset(f);
    /* invert back */
    j = 0 as i32 as uint_t;
    while j < length {
        *(*in_0).data.offset(j as isize) =
            *(*tmp).data.offset(length.wrapping_sub(j).wrapping_sub(1 as
                                                                        i32
                                                                        as
                                                                        u32)
                                    as isize);
        j = j.wrapping_add(1)
    };
}
/* * returns a pointer to feedback coefficients \f$ a_i \f$

  \param f filter object to get parameters from

  \return a pointer to the \f$ a_0 ... a_i ... a_P \f$ coefficients

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_filter_get_feedback(f:
                                                       *const aubio_filter_t)
 -> *mut lvec_t {
    return (*f).a;
}
/* * returns a pointer to feedforward coefficients \f$ b_i \f$

  \param f filter object to get coefficients from

  \return a pointer to the \f$ b_0 ... b_i ... b_P \f$ coefficients

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_filter_get_feedforward(f:
                                                          *const aubio_filter_t)
 -> *mut lvec_t {
    return (*f).b;
}
/* * get order of the filter

  \param f filter to get order from

  \return the order of the filter

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_filter_get_order(f: *const aubio_filter_t)
 -> uint_t {
    return (*f).order;
}
/* * get sampling rate of the filter

  \param f filter to get sampling rate from

  \return the sampling rate of the filter, in Hz

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_filter_get_samplerate(f:
                                                         *const aubio_filter_t)
 -> uint_t {
    return (*f).samplerate;
}
/* * get sampling rate of the filter

  \param f filter to get sampling rate from
  \param samplerate sample rate to set the filter to

  \return the sampling rate of the filter, in Hz

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_filter_set_samplerate(mut f:
                                                         *mut aubio_filter_t,
                                                     samplerate: uint_t)
 -> uint_t {
    (*f).samplerate = samplerate;
    return AUBIO_OK as i32 as uint_t;
}
/* * reset filter memory

  \param f filter object as returned by new_aubio_filter()

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_filter_do_reset(f: *mut aubio_filter_t) {
    lvec_zeros((*f).x);
    lvec_zeros((*f).y);
}
/* * create new filter object

  This function creates a new ::aubio_filter_t object, given the order of the
  filter.

  \param order order of the filter (number of coefficients)

  \return the newly created filter object

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_filter(order: uint_t)
 -> *mut aubio_filter_t {
    let mut f: *mut aubio_filter_t =
        calloc(::std::mem::size_of::<aubio_filter_t>() as u64,
               1 as i32 as u64) as *mut aubio_filter_t;
    if (order as sint_t) < 1 as i32 {
        free(f as *mut core::ffi::c_void);
        return 0 as *mut aubio_filter_t
    }
    (*f).x = new_lvec(order);
    (*f).y = new_lvec(order);
    (*f).a = new_lvec(order);
    (*f).b = new_lvec(order);
    /* by default, samplerate is not set */
    (*f).samplerate = 0 as i32 as uint_t;
    (*f).order = order;
    /* set default to identity */
    *(*(*f).a).data.offset(0 as i32 as isize) = 1.0f64;
    *(*(*f).b).data.offset(0 as i32 as isize) = 1.0f64;
    return f;
}
/* * delete a filter object

  \param f filter object to delete

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_filter(f: *mut aubio_filter_t) {
    del_lvec((*f).a);
    del_lvec((*f).b);
    del_lvec((*f).x);
    del_lvec((*f).y);
    free(f as *mut core::ffi::c_void);
}
